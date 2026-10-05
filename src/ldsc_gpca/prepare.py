"""Prepare GPCA tables and optional HapMap-filtered munging inputs from VCFs."""
import math
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

import polars as pl
from .postprocess import filename_component
from .helptext import HelpParser, PREPARE_INPUT_HELP
from .interfaces import read_manifest
from .workers import run_parallel_jobs
from .vcf_common import (ID_CHOICES, RAW_COLUMNS, QUERY, extract_table,
                         validate_and_transform, write_original_issues, merge_issue_reports,
                         selected_id, validate_ids)


def add_prepare_options(parser, *, standalone=False):
    group = parser.add_argument_group('VCF input preparation')
    group.add_argument('--p_min', type=float, default=1e-300,
                       help='Floor for P calculated from valid LP; smaller values are retained and adjusted. Default: 1e-300.')
    group.add_argument('--gpca_id_source', choices=ID_CHOICES, default='chr_pos_ref_alt',
                       help='GPCA SNPID values. Default: chr_pos_ref_alt.')
    group.add_argument('--write_munge_inputs', action='store_true',
                       help='Also write HapMap-filtered LDSC munging inputs. Default: GPCA files only; does not run munging.')
    group.add_argument('--hm3', dest='hapmap_file', metavar='REFERENCE.tsv', help='Tab-delimited file with a SNP header. Default: unset; required with --write_munge_inputs.')
    group.add_argument('--munge_id_source', choices=ID_CHOICES, default='vcf_id',
                       help='Munging SNP values; must match HapMap identifiers. Default: vcf_id (usually rsIDs).')
    worker_flag = '--n_cores' if standalone else '--prepare_workers'
    group.add_argument(worker_flag, dest='prepare_workers', type=int, default=4,
                       help='Parallel VCF preparation workers. Default: 4.')
    group.add_argument('--bcftools', default='bcftools', help='Local bcftools executable or path. Default: bcftools on PATH.')


def manifest_inputs(path):
    path = Path(path).resolve()
    columns, manifest_rows = read_manifest(path)
    if not {'traitname', 'vcf_files'}.issubset(columns) or not manifest_rows:
        raise ValueError('Preparation requires a non-empty manifest with traitname and vcf_files columns')
    rows, seen = [], set()
    for row in manifest_rows:
        name, vcf = row['traitname'], row['vcf_files']
        name = filename_component((name or '').strip())
        if name in seen:
            raise ValueError(f'Duplicate traitname: {name}')
        seen.add(name)
        if not vcf or not vcf.strip():
            raise ValueError(f'{name}: missing vcf_files path')
        source = Path(vcf.strip())
        source = (path.parent / source).resolve() if not source.is_absolute() else source.resolve()
        if not source.is_file():
            raise ValueError(f'{name}: VCF file not found: {source}')
        rows.append((name, source))
    return rows


def write_gpca(frame, name, stage, splitby_chr, gpca_id_source):
    gpca = frame.select([
        selected_id(gpca_id_source).alias('SNPID'), 'CHR', pl.col('POS').alias('BP'),
        pl.col('A1').alias('EA'), pl.col('A2').alias('OA'), pl.col('eaf_A1').alias('EAF'),
        'N', 'Z', pl.col('p').alias('P')])
    validate_ids(gpca, 'SNPID')
    if splitby_chr == 'split':
        missing = sorted(set(range(1,23)) - set(gpca['CHR'].unique().to_list()))
        if missing:
            raise ValueError(f'Missing chromosomes {missing}; split GWAMA requires files for 1–22')
        for chrom in gpca.partition_by('CHR', include_key=True):
            chrom.write_csv(stage/'gpca_inputs'/f'{name}_chr{chrom["CHR"][0]}_GenomicPCA_inputs.tsv', separator='\t')
    else:
        gpca.write_csv(stage/'gpca_inputs'/f'{name}_GenomicPCA_inputs.tsv', separator='\t')
    return gpca.height


def prepare_trait(name, vcf, stage, executable, splitby_chr, gpca_id_source, munge_id_source, hm3, p_min=1e-300):
    with tempfile.NamedTemporaryFile(dir=stage, suffix='.tsv') as temporary:
        raw = extract_table(vcf, executable, temporary.name)
        frame, issues, summary = validate_and_transform(raw, p_min, gpca_id_source)
    write_original_issues(vcf, stage/'qc'/f'{name}.csv', name, issues)
    pl.DataFrame([{'traitname':name, **summary}]).write_csv(stage/'qc'/f'{name}.summary.csv')
    if frame.is_empty():
        raise ValueError('No usable variants remain; see GPCA_Input_QC_Issues.csv for removal reasons')
    gpca_rows = write_gpca(frame, name, stage, splitby_chr, gpca_id_source)
    munge_rows = 0
    if hm3 is not None:
        munge = frame.with_columns(selected_id(munge_id_source).alias('SNP')).join(hm3, on='SNP', how='semi')
        if munge.is_empty():
            raise ValueError('No identifiers matched the HapMap SNP list; check --munge_id_source and --hm3')
        validate_ids(munge, 'SNP')
        munge.select(['SNP','CHR','POS','A1','A2','eaf_A1','beta','se','N','p']).write_csv(
            stage/'munge_inputs'/f'{name}_munge_inputs.txt', separator=' ')
        munge_rows = munge.height
    return {'traitname': name, 'vcf_files': str(vcf), 'gpca_rows': gpca_rows,
            'excluded_non_autosomal_rows': summary['excluded_non_autosomal_rows'], 'munge_rows': munge_rows,
            'p_underflow_rows': int((frame['p'] == 0).sum()), 'success': True, 'error': ''}


def prepare_shared_trait(name, vcf, stage, executable, splitby_chr, gpca_id_source,
                         munge_id_source, hm3, p_min, mode, metadata):
    from .prepare_ldsc import write_ldsc_trait
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
    result = write_ldsc_trait(raw, name, vcf, stage, munge_id_source, hm3, p_min, metadata)
    return {**result, 'gpca_rows': gpca_rows}


def prepare_inputs(manifest, outdir, splitby_chr='split', gpca_id_source='chr_pos_ref_alt',
                   write_munge_inputs=False, hapmap_file=None, munge_id_source='vcf_id',
                   prepare_workers=4, bcftools='bcftools', p_min=1e-300, *, mode='gpca'):
    if not math.isfinite(p_min) or not 0 < p_min < 1:
        raise ValueError('--p_min must be finite and strictly between 0 and 1')
    if splitby_chr not in ('split','nosplit') or gpca_id_source not in ID_CHOICES or munge_id_source not in ID_CHOICES:
        raise ValueError('Invalid layout or identifier source')
    if not isinstance(prepare_workers, int) or prepare_workers < 1:
        raise ValueError('VCF preparation workers (--n_cores for prepare; --prepare_workers within gpca) must be a positive integer')
    if mode not in ('gpca', 'ldsc', 'both') or (mode != 'gpca' and write_munge_inputs):
        raise ValueError('Use --mode ldsc or --mode both without legacy --write_munge_inputs')
    shared = mode != 'gpca'
    if (write_munge_inputs or shared) != bool(hapmap_file):
        raise ValueError('--hm3 is required for LDSC preparation or --write_munge_inputs; omit it for GPCA-only preparation')
    executable = shutil.which(bcftools)
    if not executable:
        raise ValueError(f'Local bcftools executable not found: {bcftools}')
    rows = manifest_inputs(manifest)
    hm3 = None
    if write_munge_inputs or shared:
        with open(hapmap_file) as stream:
            if 'SNP' not in stream.readline().rstrip('\r\n').split('\t'):
                raise ValueError('--hm3 must be tab-separated with a SNP header')
        hm3 = pl.read_csv(hapmap_file, separator='\t', columns=['SNP'], schema_overrides={'SNP': pl.String}).drop_nulls().unique()
        if hm3.is_empty():
            raise ValueError('HapMap SNP list is empty')
    if shared:
        from .prepare_ldsc import sample_metadata
        metadata = sample_metadata(manifest)
    outdir = Path(outdir).resolve()
    report_sets = [('qc', 'GPCA_Input')] if mode != 'ldsc' else []
    if shared:
        report_sets.append(('ldsc_qc', 'LDSC_Input'))
    reports = [f'{prefix}_QC_{kind}.csv' for _, prefix in report_sets for kind in ('Issues', 'Summary')]
    outputs = ['Preparation_Status.csv', 'Preparation_Settings.json', *reports]
    if mode != 'ldsc':
        outputs.insert(0, 'gpca_inputs')
    if write_munge_inputs or shared:
        outputs.append('munge_inputs')
    if shared:
        outputs.append('Prepared_LDSC_Manifest.csv')
    for name in outputs:
        path = outdir/name
        if path.exists() or path.is_symlink():
            raise FileExistsError(f'Refusing to overwrite {path}; use a fresh preparation outdir or pass existing --gpca_input_folder to gpca')
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=outdir) as temporary:
        stage = Path(temporary)
        for folder in ['gpca_inputs'] * (mode != 'ldsc') + ['munge_inputs'] * bool(write_munge_inputs or shared) + [folder for folder, _ in report_sets]:
            (stage/folder).mkdir()
        statuses = {}
        jobs = [(name, prepare_shared_trait if shared else prepare_trait,
                 (name, vcf, stage, executable, splitby_chr, gpca_id_source, munge_id_source, hm3, p_min) +
                 ((mode, metadata[name]) if shared else ())) for name, vcf in rows]
        values = run_parallel_jobs(jobs,
            prepare_workers, stage='VCF preparation', collect_failures=True,
            status_path=outdir/'Preparation_Worker_Attempts.csv')
        for (name, vcf), value in zip(rows, values):
            if isinstance(value, Exception):
                detail = getattr(value, 'stderr', '') or str(value)
                statuses[name] = {'traitname':name, 'vcf_files':str(vcf), 'gpca_rows':0,
                                  'excluded_non_autosomal_rows':0, 'munge_rows':0, 'p_underflow_rows':0,
                                  'success':False, 'error':str(detail)}
            else:
                statuses[name] = value
                print(f'Prepared {name}: {value["gpca_rows"]} GPCA variants, {value["munge_rows"]} raw LDSC variants', flush=True)
        ordered = [statuses[name] for name,_ in rows]
        for folder, prefix in report_sets:
            summaries = []
            for status in ordered:
                path = stage/folder/f'{status["traitname"]}.summary.csv'
                summary = pl.read_csv(path, schema_overrides={'traitname':pl.String}).to_dicts()[0] if path.exists() else {
                    'traitname':status['traitname'], 'input_rows':None, 'retained_rows':None,
                    'removed_rows':None, 'p_adjusted_rows':None, 'excluded_non_autosomal_rows':None}
                summaries.append({**summary, 'success':status['success'], 'error':status['error']})
            pl.DataFrame(summaries, infer_schema_length=None).write_csv(stage/f'{prefix}_QC_Summary.csv')
            merge_issue_reports([stage/folder/f'{name}.csv' for name,_ in rows
                                 if (stage/folder/f'{name}.csv').exists()], stage/f'{prefix}_QC_Issues.csv')
        failed = [s for s in ordered if not s['success']]
        if failed:
            pl.DataFrame(ordered).write_csv(outdir/'Preparation_Status.csv')
            for report in reports:
                os.rename(stage/report, outdir/report)
            raise ValueError('Analysis stopped: VCF preparation failed after 2 attempts; no prepared tables published. '
                             'See Preparation_Worker_Attempts.csv.\n' + '\n'.join(f'{s["traitname"]}: {s["error"]}' for s in failed))
        pl.DataFrame(ordered).write_csv(stage/'Preparation_Status.csv')
        settings = {'manifest':str(Path(manifest).resolve()), 'splitby_chr':splitby_chr,
                    'gpca_id_source':gpca_id_source, 'munge_id_source':munge_id_source,
                    'write_munge_inputs':write_munge_inputs, 'hapmap_file':str(Path(hapmap_file).resolve()) if hapmap_file else None,
                    'N_source':'FORMAT/NEF', 'Z_source':'FORMAT/ES / FORMAT/SE',
                    'P_source':'10 ** (-FORMAT/LP), floored at p_min for valid LP', 'p_min':p_min,
                    'duplicate_policy':'Largest valid LP, first original row on ties',
                    'QC_issues_source':'Original VCF fields, unchanged; union of headers across samples',
                    'bcftools':executable}
        if shared:
            from .prepare_ldsc import write_ldsc_manifest
            write_ldsc_manifest(stage/'Prepared_LDSC_Manifest.csv', rows, metadata, outdir)
            settings.update(mode=mode, LDSC_N_source='manifest N override, otherwise FORMAT/NEF',
                            LDSC_INFO_source='FORMAT/SI', LDSC_effect_source='FORMAT/ES for ALT',
                            LDSC_QC='basic QC and HapMap IDs; INFO/MAF thresholds applied only during munging')
        (stage/'Preparation_Settings.json').write_text(json.dumps(settings,indent=2)+'\n')
        published = []
        try:
            for name in outputs:
                if (outdir/name).exists():
                    raise FileExistsError(f'Output appeared during preparation: {outdir/name}')
                os.rename(stage/name, outdir/name)
                published.append(name)
        except Exception:
            for name in reversed(published):
                os.rename(outdir/name, stage/name)
            raise
    folder = outdir/('munge_inputs' if mode == 'ldsc' else 'gpca_inputs')
    print(f'Prepared inputs: {folder}')
    return folder


def preparation_kwargs(opts):
    return {name: getattr(opts,name) for name in ['gpca_id_source','write_munge_inputs','hapmap_file',
                                                'munge_id_source','prepare_workers','bcftools','p_min']}


def build_parser():
    parser = HelpParser(prog='ldsc-gpca prepare', description=__doc__,
                        usage='%(prog)s --input MANIFEST.csv --outdir DIRECTORY [options]', epilog=PREPARE_INPUT_HELP)
    required = parser.add_argument_group('Required inputs')
    layout = parser.add_argument_group('Output layout')
    required.add_argument('--input', metavar='MANIFEST.csv', required=True, help='CSV with traitname and vcf_files. Relative VCF paths resolve beside this manifest.')
    required.add_argument('--outdir', metavar='DIRECTORY', required=True, help='Output directory; creates gpca_inputs and optionally munge_inputs.')
    layout.add_argument('--splitby_chr', choices=['split','nosplit'], default='split', help='GPCA layout only; ignored in ldsc-only mode. Default: split (requires all chromosomes 1–22 per trait); nosplit writes one autosomal table.')
    add_prepare_options(parser, standalone=True)
    from .prepare_ldsc import add_mode_options
    add_mode_options(parser)
    return parser


def main(argv=None):
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        parser.print_help()
        return 0
    opts = parser.parse_args(argv)
    try:
        from .prepare_ldsc import validate_mode, munge_prepared
        command = validate_mode(opts, argv)
        prepare_inputs(opts.input, opts.outdir, splitby_chr=opts.splitby_chr, mode=opts.mode, **preparation_kwargs(opts))
        if command is not None:
            munge_prepared(opts, command)
    except (OSError, ValueError, pl.exceptions.PolarsError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0
