"""Shared raw LDSC preparation and optional backend munging; no regressions."""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import warnings
from importlib.resources import files

import polars as pl

from .interfaces import read_manifest
from .restart import file_digest, write_checkpoint
from .utils import optional_prevalence
from .workers import run_parallel_jobs


def add_mode_options(parser):
    group = parser.add_argument_group('What to prepare')
    parser._action_groups.remove(group)
    parser._action_groups.insert(3, group)
    group.add_argument('--mode', choices=['gpca', 'ldsc', 'both'], default='gpca',
                       help='GPCA/GWAMA tables, shared LDSC tables plus munging, or both.')
    group.add_argument('--raw_only', action='store_true',
                       help='With ldsc/both, stop after writing shared raw LDSC tables; skip munging.')
    group.add_argument('--munge_backend', choices=['python', 'genomicsem'], default='python',
                       help='Program used for standalone preparation munging only. ldsc/pipeline use --ldsc_backend instead. Select genomicsem for later GenomicSEM LDSC reuse. Omit with --raw_only or gpca mode.')
    group = parser.add_argument_group('Munging filters (ldsc/both, without --raw_only)')
    group.add_argument('--info_filter', type=float, default=.9, help='INFO threshold for the selected munger; does not filter GPCA tables.')
    group.add_argument('--maf_filter', type=float, default=.01, help='MAF threshold; Python uses >, GenomicSEM uses >=. Does not filter GPCA tables.')
    group = parser.add_argument_group('Munging runtimes')
    group.add_argument('--conda_executable', default=os.environ.get('CONDA_EXE', 'conda'),
                       help='Python munging: Conda executable. Default: CONDA_EXE, otherwise conda on PATH.')
    target = group.add_mutually_exclusive_group()
    target.add_argument('--ldsc_env', help='Python munging environment name. Default: ldsc-cbiit unless a prefix is configured.')
    target.add_argument('--ldsc_env_prefix', help='Python munging environment path. Default: LDSC_GPCA_LDSC_PREFIX, otherwise named environment.')
    group.add_argument('--rscript', default='Rscript', help='GenomicSEM munging: Rscript executable with GenomicSEM installed.')
    parser._option_string_actions['--hm3'].help = 'Tab-separated HapMap reference; required for ldsc/both or legacy --write_munge_inputs. SNP header for raw export; SNP,A1,A2 for munging. Default: unset.'
    parser._option_string_actions['--write_munge_inputs'].help = 'Legacy GPCA-plus-raw export; original schema without INFO. Default: off. Do not combine with LDSC modes.'
    parser._option_string_actions['--n_cores'].help = 'Parallel workers for preparation and optional munging. Default: 4.'
    parser._option_string_actions['--outdir'].help = 'Output directory; creates only requested input folders, QC and optional munged files.'
    parser.description = 'Prepare GPCA/GWAMA inputs, shared raw LDSC tables and optional munged files from GWAS VCFs.'
    parser.epilog = '''INPUT FILE CONTRACT
  --input: CSV with unique traitname and vcf_files; one GWAS sample per VCF.
  Relative VCF paths resolve beside the manifest. Plain VCF and .vcf.gz accepted.
  Every mode needs coordinates, REF/ALT and FORMAT AF,ES,SE,LP. LP is -log10(P).
  --mode gpca (default): FORMAT/NEF is also required. Writes nine-column
  GPCA tables: SNPID,CHR,BP,EA,OA,EAF,N,Z,P. Z=ES/SE and N=NEF.
  GPCA split mode requires valid records on every chromosome 1-22 per trait;
  nosplit writes one autosomal file per trait. LDSC tables are always unsplit.

LDSC MODES
  --mode ldsc or --mode both: shared raw TSVs (SNP,CHR,BP,A1,A2,EAF,BETA,SE,P,N,INFO).
  Requires FORMAT/SI for INFO. N comes from NEF unless manifest N overrides it.
  Binary rows require explicit total N and both sample_prevalence and
  population_prevalence. Quantitative rows leave both prevalences blank.
  In both mode, GWAMA still uses NEF; the N override is LDSC-only.
  --raw_only stops before munging; otherwise Python munging is the default.
  --munge_backend genomicsem selects GenomicSEM R munging. Backend choices do not
  change the shared raw format. Munged core columns: SNP,N,Z,A1,A2.
  For later ldsc/pipeline --ldsc_backend genomicsem (or genomicsem ldsc) reuse, select
  --munge_backend genomicsem and the same --hm3 allele reference for all traits.
  Do not directly reuse Python-munged files in GenomicSEM: matching the
  same HapMap reference does not guarantee the required allele orientation.
  GenomicSEM preparation requires the pinned installed package. Keep the complete
  preparation directory for reuse: manifest, settings, sidecars, and copied allele
  reference. Reuse verifies completion, checksums, N/prevalences and allele order.
  New LDSC modes write Prepared_LDSC_Manifest.csv and LDSC_Input_QC_* reports.
  Raw tables remain available if munging fails; final munged/ is published only
  after all traits succeed. Failed jobs receive one retry (two total attempts).
  No LDSC regression, PCA or GWAMA is run. No LD scores are needed here.
  Existing --write_munge_inputs keeps its legacy GPCA-plus-raw behavior and schema;
  do not combine it with --mode ldsc or --mode both or --raw_only.
  Legacy raw schema: SNP,CHR,POS,A1,A2,eaf_A1,beta,se,N,p (space-separated).
  Invalid records are removed and audited. Original VCFs are never changed.
'''


def sample_metadata(manifest):
    """Use declared total N for binary data; never infer it from VCF NEF."""
    _, rows = read_manifest(manifest)
    result = {}
    for row in rows:
        name = row['traitname'].strip()
        prev = {key: optional_prevalence(row.get(key), f'{name}: {key}')
                for key in ('sample_prevalence', 'population_prevalence')}
        if (prev['sample_prevalence'] is None) != (prev['population_prevalence'] is None):
            raise ValueError(f'{name}: supply both prevalences or neither')
        value = (row.get('N') or '').strip()
        n = None if value.lower() in ('', '.', 'na', 'nan') else float(value)
        if n is not None and (not math.isfinite(n) or n <= 0):
            raise ValueError(f'{name}: manifest N must be finite and positive')
        if prev['population_prevalence'] is not None and n is None:
            raise ValueError(f'{name}: binary LDSC preparation requires explicit total N in the manifest and both prevalences')
        result[name] = {'N': n, **prev}
    return result


def write_ldsc_manifest(path, rows, metadata, outdir, *, munged=False):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', newline='', dir=Path(path).parent,
                                         prefix='.manifest-', delete=False) as stream:
            temporary = stream.name
            writer = csv.DictWriter(stream, fieldnames=['traitname', 'ref', 'sumstats_file', 'munged_file',
                                                       'N', 'sample_prevalence', 'population_prevalence'])
            writer.writeheader()
            for name, _ in rows:
                writer.writerow({'traitname': name, 'ref': 'yes', **metadata[name],
                    'sumstats_file': str(outdir/'munge_inputs'/f'{name}_munge_inputs.tsv'),
                    'munged_file': str(outdir/'munged'/f'{name}.sumstats.gz') if munged else ''})
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


def prepare_shared_trait(name, vcf, stage, executable, splitby_chr, gpca_id_source,
                         munge_id_source, hm3, p_min, mode, metadata):
    from .prepare import (QUERY, RAW_COLUMNS, extract_table, validate_and_transform,
                          write_original_issues, write_gpca, selected_id, validate_ids)
    query = QUERY.replace(r'\t%NEF', r'\t%NEF\t%SI')
    if mode == 'ldsc' and metadata['N'] is not None:
        query = query.replace('%NEF', format(metadata['N'], '.17g'))
    with tempfile.NamedTemporaryFile(dir=stage, suffix='.tsv') as temporary:
        raw = extract_table(vcf, executable, temporary.name, query=query, columns=[*RAW_COLUMNS, 'INFO'])
    gpca_rows = 0
    if mode == 'both':
        gpca, issues, summary = validate_and_transform(raw, p_min, gpca_id_source)
        write_original_issues(vcf, stage/'qc'/f'{name}.csv', name, issues)
        pl.DataFrame([{'traitname': name, **summary}]).write_csv(stage/'qc'/f'{name}.summary.csv')
        if gpca.is_empty():
            raise ValueError(f'{name}: no usable GPCA records; see GPCA_Input_QC_Issues.csv')
        gpca_rows = write_gpca(gpca, name, stage, splitby_chr, gpca_id_source)
    if metadata['N'] is not None:
        raw = raw.with_columns(pl.lit(metadata['N']).alias('N'))
    frame, issues, summary = validate_and_transform(raw, p_min, munge_id_source, require_info=True)
    frame = frame.with_columns(selected_id(munge_id_source).alias('SNP'))
    removed = frame.join(hm3, on='SNP', how='anti').select('_row').with_columns(
        pl.lit('SNP not in HapMap reference').alias('QC_reason'), pl.lit('removed').alias('QC_action'))
    issues = pl.concat([issues.join(removed.select('_row'), on='_row', how='anti'), removed]).sort('_row')
    frame = frame.join(hm3, on='SNP', how='semi')
    summary.update(retained_rows=frame.height, removed_rows=summary['input_rows']-frame.height,
                   hapmap_removed_rows=removed.height,
                   p_adjusted_rows=issues.filter(pl.col('QC_action') == 'p_adjusted').height)
    write_original_issues(vcf, stage/'ldsc_qc'/f'{name}.csv', name, issues)
    pl.DataFrame([{'traitname': name, **summary, 'N_source': 'manifest_N' if metadata['N'] is not None else 'FORMAT/NEF',
                   'N_override': metadata['N'], 'p_min': p_min}]).write_csv(stage/'ldsc_qc'/f'{name}.summary.csv')
    if frame.is_empty():
        raise ValueError(f'{name}: no usable LDSC records matched HapMap; see LDSC_Input_QC_Issues.csv')
    validate_ids(frame, 'SNP')
    frame.select('SNP', 'CHR', pl.col('POS').alias('BP'), 'A1', 'A2', pl.col('eaf_A1').alias('EAF'),
                 pl.col('beta').alias('BETA'), pl.col('se').alias('SE'), pl.col('p').alias('P'), 'N', 'INFO').write_csv(
        stage/'munge_inputs'/f'{name}_munge_inputs.tsv', separator='\t')
    return {'traitname': name, 'vcf_files': str(vcf), 'gpca_rows': gpca_rows,
            'excluded_non_autosomal_rows': summary['excluded_non_autosomal_rows'], 'munge_rows': frame.height,
            'p_underflow_rows': 0, 'success': True, 'error': ''}


def validate_mode(opts, argv):
    """Reject irrelevant options and missing runtimes before publishing inputs."""
    flags = {arg.split('=')[0] for arg in argv}
    munging = opts.mode != 'gpca' and not opts.raw_only
    python_flags = {'--conda_executable', '--ldsc_env', '--ldsc_env_prefix'}
    munging_flags = python_flags | {'--munge_backend', '--rscript', '--info_filter', '--maf_filter'}
    if opts.raw_only and opts.mode == 'gpca':
        raise ValueError('--raw_only requires --mode ldsc or both')
    if not munging and flags & munging_flags:
        raise ValueError('Munging options require --mode ldsc or --mode both without --raw_only: ' + ', '.join(sorted(flags & munging_flags)))
    if not munging:
        return None
    if not 0 <= opts.info_filter <= 1 or not 0 <= opts.maf_filter <= .5:
        raise ValueError('--info_filter must be in [0,1]; --maf_filter in [0,0.5]')
    if not opts.hapmap_file:
        raise ValueError('--hm3 is required for LDSC preparation')
    with open(opts.hapmap_file) as stream:
        if not {'SNP', 'A1', 'A2'}.issubset(stream.readline().split()):
            raise ValueError('Munging --hm3 needs SNP,A1,A2 headers')
    for name in ('munged', 'munge_logs', 'Preparation_Munging_Worker_Attempts.csv', 'Preparation_Munging_Status.csv'):
        if (Path(opts.outdir)/name).exists() or (Path(opts.outdir)/name).is_symlink():
            raise FileExistsError(f'Refusing to overwrite {Path(opts.outdir)/name}')
    if opts.munge_backend == 'python':
        if '--rscript' in flags:
            raise ValueError('--rscript applies only to --munge_backend genomicsem')
        from .ldsc_runtime import ldsc_munging_command
        command = ldsc_munging_command(conda=opts.conda_executable, environment=opts.ldsc_env or 'ldsc-cbiit',
            prefix=None if opts.ldsc_env else opts.ldsc_env_prefix or os.environ.get('LDSC_GPCA_LDSC_PREFIX'))
        probe = command + ['--help']
    else:
        if flags & python_flags:
            raise ValueError('Conda/LDSC environment options apply only to Python munging')
        command = [opts.rscript, str(files('ldsc_gpca').joinpath('r', 'genomicsem', 'munge_entry.R'))]
        probe = command + ['--check']
    check = subprocess.run(probe, capture_output=True, text=True)
    if check.returncode:
        raise ValueError(f'{opts.munge_backend} munging runtime unavailable: {check.stderr or check.stdout}')
    return command


def munge_prepared(opts, command):
    """Keep raw inputs; publish munged files together only after all traits succeed."""
    out = Path(opts.outdir).resolve()
    settings_path = out/'Preparation_Settings.json'
    settings = json.loads(settings_path.read_text())
    settings.update(munge_backend=opts.munge_backend, info_filter=opts.info_filter,
                    maf_filter=opts.maf_filter, munging_status='running')
    write_checkpoint(settings_path, settings)
    try:
        _munge_and_publish(opts, command, out)
        settings['munging_manifest_sha256'] = file_digest(out/'Prepared_LDSC_Manifest.csv')
        settings['munging_artifacts'] = {path.name: file_digest(path) for path in sorted((out/'munged').iterdir()) if path.is_file()}
        settings['munging_status'] = 'completed'
        write_checkpoint(settings_path, settings)
    except BaseException as error:
        settings.update(munging_status='failed', munging_error=str(error))
        try:
            write_checkpoint(settings_path, settings)
        except OSError as status_error:
            warnings.warn(f'Could not record failed munging status: {status_error}; preparation is not complete.')
        raise
    print(f'Munged inputs: {out/"munged"}')


def _munge_and_publish(opts, command, out):
    from .ldsc_munge import FLOAT_FORMAT
    from .genomicsem_inputs import (REFERENCE_FILE, read_allele_reference, validate_alleles,
                                   write_input_audit, InputCompatibilityError)
    _, rows = read_manifest(out/'Prepared_LDSC_Manifest.csv')
    logs = out/'munge_logs'
    logs.mkdir()
    with tempfile.TemporaryDirectory(prefix='.munge-', dir=out) as temporary:
        stage = Path(temporary)
        reference_path = Path(opts.hapmap_file).resolve()
        if opts.munge_backend == 'genomicsem':
            # Munge every trait against the same immutable snapshot, retained for reuse.
            reference_path = stage/REFERENCE_FILE
            shutil.copyfile(opts.hapmap_file, reference_path)
            reference_sha256 = file_digest(reference_path)
            reference = read_allele_reference(reference_path)
        def worker(row):
            name = row['traitname']
            with tempfile.TemporaryDirectory(prefix='trait-', dir=stage) as scratch:
                prefix = Path(scratch)/name
                if opts.munge_backend == 'python':
                    args = ['--sumstats', row['sumstats_file'], '--out', str(prefix), '--snp', 'SNP',
                        '--a1', 'A1', '--a2', 'A2', '--p', 'P', '--N-col', 'N', '--frq', 'EAF', '--info', 'INFO',
                        '--signed-sumstats', 'BETA,0', '--info-min', str(opts.info_filter), '--maf-min', str(opts.maf_filter),
                        '--merge-alleles', str(Path(opts.hapmap_file).resolve())]
                else:
                    args = [row['sumstats_file'], str(prefix), str(reference_path),
                            str(opts.info_filter), str(opts.maf_filter)]
                (logs/f'{name}.command.json').write_text(json.dumps(command+args, indent=2)+'\n')
                with (logs/f'{name}.console.log').open('a') as handle:
                    result = subprocess.run(command+args, stdout=handle, stderr=subprocess.STDOUT)
                for path in Path(scratch).glob('*.log'):
                    shutil.copyfile(path, logs/path.name)
                if result.returncode:
                    raise ValueError(f'{name}: {opts.munge_backend} munging exited {result.returncode}; see {logs/name}.console.log')
                path = prefix.with_name(name+'.sumstats.gz')
                frame = pl.read_csv(path, separator='\t', schema_overrides={'SNP': pl.String, 'A1': pl.String, 'A2': pl.String})
                if not {'SNP','A1','A2','N','Z'}.issubset(frame.columns):
                    raise ValueError(f'{name}: munged output lacks SNP,A1,A2,N,Z')
                valid = frame.drop_nulls(['SNP','A1','A2','N','Z'])
                if valid.is_empty() or not valid.select((pl.col('N').is_finite() & (pl.col('N') > 0) & pl.col('Z').is_finite()).all()).item():
                    raise ValueError(f'{name}: munged output has no usable rows or invalid N/Z')
                metadata = {'preparation_method': 'shared', 'munge_backend': opts.munge_backend,
                    'sample_size_column': 'N_TOTAL' if row['population_prevalence'] else 'NEF',
                    'sample_prevalence': float(row['sample_prevalence']) if row['sample_prevalence'] else None,
                    'sha256': file_digest(path), 'info_filter': opts.info_filter, 'maf_filter': opts.maf_filter,
                    'N_source': 'manifest_N' if row['N'] else 'FORMAT/NEF',
                    'float_format': FLOAT_FORMAT if opts.munge_backend == 'python' else 'GenomicSEM'}
                metadata.update(schema_version=1, traitname=name,
                    population_prevalence=float(row['population_prevalence']) if row['population_prevalence'] else None,
                    N_override=float(row['N']) if row['N'] else None,
                    N_convention='total' if row['population_prevalence'] else 'manifest_N' if row['N'] else 'NEF')
                if opts.munge_backend == 'genomicsem':
                    with prefix.with_name(name+'.runtime.tsv').open(newline='') as stream:
                        runtimes = list(csv.DictReader(stream, delimiter='\t'))
                    if len(runtimes) != 1:
                        raise ValueError(f'{name}: missing GenomicSEM runtime identity')
                    metadata.update(reference_file=REFERENCE_FILE, reference_sha256=reference_sha256, runtime=runtimes[0])
                    try:
                        audit = validate_alleles(path, name, reference, n_override=metadata['N_override'])
                    except InputCompatibilityError as error:
                        write_input_audit(logs/f'{name}.input_qc.csv', [error.audit])
                        raise
                    audit.update(sumstats_sha256=metadata['sha256'], reference_sha256=reference_sha256)
                    write_input_audit(logs/f'{name}.input_qc.csv', [audit])
                path.rename(stage/path.name)
                (stage/f'{name}.prevalence.json').write_text(json.dumps(metadata, indent=2)+'\n')
                return {'traitname': name, 'output_rows': frame.height, 'usable_rows': valid.height,
                        'missing_rows': frame.height-valid.height, 'success': True, 'error': ''}
        values = run_parallel_jobs([(r['traitname'], worker, (r,)) for r in rows], opts.prepare_workers,
            stage='Preparation munging', status_path=out/'Preparation_Munging_Worker_Attempts.csv', collect_failures=True)
        statuses = [value if not isinstance(value, Exception) else {'traitname': row['traitname'],
                    'success': False, 'error': str(value)} for row, value in zip(rows, values)]
        pl.DataFrame(statuses, infer_schema_length=None).write_csv(out/'Preparation_Munging_Status.csv')
        failed = [r for r in statuses if not r['success']]
        if failed:
            raise ValueError('Munging failed after 2 attempts; raw inputs retained, no munged/ published.\n' +
                             '\n'.join(r['error'] for r in failed))
        stage.rename(out/'munged')
    metadata = sample_metadata(out/'Prepared_LDSC_Manifest.csv')
    write_ldsc_manifest(out/'Prepared_LDSC_Manifest.csv', [(r['traitname'], None) for r in rows], metadata, out, munged=True)
