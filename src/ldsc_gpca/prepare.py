"""Prepare GPCA tables and optional HapMap-filtered munging inputs from VCFs."""
import csv
import gzip
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import polars as pl
from .postprocess import filename_component
from .helptext import HelpParser, PREPARE_INPUT_HELP

ID_CHOICES = ('chr_pos_ref_alt', 'vcf_id')
RAW_COLUMNS = ['SNP', 'CHR', 'POS', 'A1', 'A2', 'eaf_A1', 'beta', 'se', 'LP', 'N']
QUERY = r'%ID\t%CHROM\t%POS\t%ALT\t%REF[\t%AF\t%ES\t%SE\t%LP\t%NEF]\n'


def add_prepare_options(parser):
    group = parser.add_argument_group('VCF input preparation')
    group.add_argument('--p-min', type=float, default=1e-300,
                       help='Floor for P calculated from valid LP; smaller values are retained and adjusted. Default: 1e-300.')
    group.add_argument('--gpca-id-source', choices=ID_CHOICES, default='chr_pos_ref_alt',
                       help='GPCA SNPID values. Default: chr_pos_ref_alt.')
    group.add_argument('--write-munge-inputs', action='store_true',
                       help='Also write HapMap-filtered LDSC munging inputs. Default: GPCA files only; does not run munging.')
    group.add_argument('--hapmap-file', help='Tab-delimited file with a SNP header. Required only with --write-munge-inputs.')
    group.add_argument('--munge-id-source', choices=ID_CHOICES, default='vcf_id',
                       help='Munging SNP values; must match HapMap identifiers. Default: vcf_id (usually rsIDs).')
    group.add_argument('--prepare-workers', type=int, default=4,
                       help='Parallel VCF preparation workers. Default: 4.')
    group.add_argument('--bcftools', default='bcftools', help='Local bcftools executable or path. Default: bcftools on PATH.')


def manifest_inputs(path):
    path = Path(path).resolve()
    frame = pl.read_csv(path, schema_overrides={'traitname': pl.String, 'vcf_files': pl.String})
    if not {'traitname', 'vcf_files'}.issubset(frame.columns) or frame.is_empty():
        raise ValueError('Preparation requires a non-empty manifest with traitname and vcf_files columns')
    rows, seen = [], set()
    for name, vcf in frame.select('traitname', 'vcf_files').iter_rows():
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


def extract_table(vcf, executable, temporary):
    samples = subprocess.run([executable, 'query', '-l', str(vcf)], check=True,
                             capture_output=True, text=True).stdout.splitlines()
    if len(samples) != 1:
        raise ValueError(f'{vcf.name}: expected exactly one GWAS sample; found {len(samples)}')
    with open(temporary, 'w') as handle:
        subprocess.run([executable, 'query', '-f', QUERY, str(vcf)],
                       stdout=handle, stderr=subprocess.PIPE, text=True, check=True)
    if Path(temporary).stat().st_size == 0:
        raise ValueError(f'{vcf.name}: empty VCF query output')
    return pl.read_csv(temporary, separator='\t', has_header=False,
                       schema={name: pl.String for name in RAW_COLUMNS}, null_values='.')


def validate_and_transform(frame, p_min=1e-300, gpca_id_source='chr_pos_ref_alt'):
    """Filter bad rows, retaining internal row indices solely for original-record audits."""
    if not math.isfinite(p_min) or not 0 < p_min < 1:
        raise ValueError('--p-min must be finite and strictly between 0 and 1')
    original_rows = frame.height
    frame = frame.with_row_index('_row')
    frame = frame.with_columns(pl.col('CHR').str.replace(r'(?i)^chr', '').cast(pl.Int64, strict=False))
    frame = frame.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in ['POS','eaf_A1','beta','se','LP','N']])
    frame = frame.with_columns(pl.col('A1','A2').str.to_uppercase(),
                               (pl.col('beta') / pl.col('se')).alias('Z'))
    checks = [(pl.col('CHR').is_between(1,22), 'Non-autosomal or invalid CHR'),
              ((pl.col('POS') > 0) & (pl.col('POS') < 2**63) & (pl.col('POS') % 1 == 0), 'Invalid POS'),
              (pl.col('eaf_A1').is_between(0,1), 'EAF outside [0,1]'),
              (pl.col('se') > 0, 'Require SE > 0'), (pl.col('N') > 0, 'Require NEF > 0'),
              (pl.col('LP') >= 0, 'Require LP >= 0'),
              (pl.col('A1') != pl.col('A2'), 'Identical REF and ALT')]
    checks += [(pl.col(c).is_finite(), f'{c}: missing, non-numeric or non-finite')
               for c in ['POS','eaf_A1','beta','se','LP','N','Z']]
    checks += [(pl.col(c).str.contains(r'^[ACGTN]+$'), f'{c}: invalid sequence allele') for c in ['A1','A2']]
    if gpca_id_source == 'vcf_id':
        checks.append((~pl.col('SNP').str.contains(r'\s') & ~pl.col('SNP').is_in(['','.']), 'Invalid VCF identifier'))
    frame = frame.with_columns(pl.concat_str([
        pl.when(valid.fill_null(False)).then(pl.lit(None, dtype=pl.String)).otherwise(pl.lit(reason))
        for valid, reason in checks], separator='; ', ignore_nulls=True).alias('QC_reason'))
    issues = frame.filter(pl.col('QC_reason') != '').select('_row','QC_reason').with_columns(pl.lit('removed').alias('QC_action'))
    excluded = frame.filter(~pl.col('CHR').is_between(1,22).fill_null(False)).height
    frame = frame.filter(pl.col('QC_reason') == '').with_columns(pl.col('POS').cast(pl.Int64))
    frame = frame.with_columns([
        pl.concat_str(['CHR','POS','A2','A1'], separator='_').alias('coordinate_id'),
        pl.lit(10.0).pow(-pl.col('LP').clip(upper_bound=-math.log10(p_min))).clip(lower_bound=p_min).alias('p')])
    keys = ['coordinate_id'] + (['SNP'] if gpca_id_source == 'vcf_id' else [])
    duplicate_keys = [key for key in keys if frame[key].n_unique() != frame.height]
    # Only rank when necessary; use original LP and original row for exact ties.
    if duplicate_keys:
        frame = frame.sort(['LP','_row'], descending=[True,False])
    for key in duplicate_keys:
        duplicate = ~pl.col(key).is_first_distinct()
        issues = pl.concat([issues, frame.filter(duplicate).select('_row').with_columns(
            pl.lit(f'Duplicate {key}; retained largest valid LP, first original row on ties').alias('QC_reason'),
            pl.lit('duplicate_removed').alias('QC_action'))])
        frame = frame.filter(~duplicate)
    adjusted = frame.filter(pl.col('LP') > -math.log10(p_min)).select('_row').with_columns(
        pl.lit(f'P below {p_min:g}; floored to {p_min:g}').alias('QC_reason'), pl.lit('p_adjusted').alias('QC_action'))
    issues = pl.concat([issues, adjusted]).sort('_row')
    summary = {'input_rows':original_rows, 'retained_rows':frame.height,
               'removed_rows':original_rows-frame.height, 'p_adjusted_rows':adjusted.height,
               'excluded_non_autosomal_rows':excluded}
    return frame.sort(['CHR','POS','_row']), issues, summary


def write_original_issues(vcf, destination, name, issues):
    """Copy original VCF field strings, without numerical parsing or rewriting."""
    affected = {row: (action, reason) for row, reason, action in issues.iter_rows()}
    with open(vcf, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    with (gzip.open if compressed else open)(vcf, 'rt') as source, open(destination, 'w', newline='') as target:
        writer = csv.writer(target)
        index = 0
        for line in source:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                fields = line.rstrip('\r\n').split('\t')
                if set(fields) & {'traitname','QC_action','QC_reason'} or len(set(fields)) != len(fields):
                    raise ValueError('VCF column names must be unique and not collide with QC report columns')
                writer.writerow(fields + ['traitname','QC_action','QC_reason'])
                if not affected:
                    break
                continue
            if index in affected:
                writer.writerow(line.rstrip('\r\n').split('\t') + [name, *affected.pop(index)])
                if not affected:
                    break
            index += 1
        if affected:
            raise ValueError('Original VCF ended before all affected records could be reported')


def merge_issue_reports(paths, destination):
    """Stream original field strings; only the union of headers stays in memory."""
    tail = ['traitname','QC_action','QC_reason']
    columns = {}
    for path in paths:
        with open(path, newline='') as source:
            columns.update(dict.fromkeys(next(csv.reader(source))))
    header = [c for c in columns if c not in tail] + tail
    with open(destination, 'w', newline='') as target:
        writer = csv.DictWriter(target, fieldnames=header)
        writer.writeheader()
        for path in paths:
            with open(path, newline='') as source:
                writer.writerows(csv.DictReader(source))


def selected_id(source):
    return pl.col('coordinate_id' if source == 'chr_pos_ref_alt' else 'SNP')


def validate_ids(frame, column):
    values = frame[column]
    if (values.is_null().any() or values.str.strip_chars().is_in(['', '.']).any()
            or values.str.contains(r'\s').any() or values.n_unique() != frame.height):
        raise ValueError(f'{column} identifiers must be present, whitespace-free and unique')


def prepare_trait(name, vcf, stage, executable, splitby_chr, gpca_id_source, munge_id_source, hm3, p_min=1e-300):
    with tempfile.NamedTemporaryFile(dir=stage, suffix='.tsv') as temporary:
        raw = extract_table(vcf, executable, temporary.name)
        frame, issues, summary = validate_and_transform(raw, p_min, gpca_id_source)
    write_original_issues(vcf, stage/'qc'/f'{name}.csv', name, issues)
    pl.DataFrame([{'traitname':name, **summary}]).write_csv(stage/'qc'/f'{name}.summary.csv')
    if frame.is_empty():
        raise ValueError('No usable variants remain; see GPCA_Input_QC_Issues.csv for removal reasons')
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
    munge_rows = 0
    if hm3 is not None:
        munge = frame.with_columns(selected_id(munge_id_source).alias('SNP')).join(hm3, on='SNP', how='semi')
        if munge.is_empty():
            raise ValueError('No identifiers matched the HapMap SNP list; check --munge-id-source and --hapmap-file')
        validate_ids(munge, 'SNP')
        munge.select(['SNP','CHR','POS','A1','A2','eaf_A1','beta','se','N','p']).write_csv(
            stage/'munge_inputs'/f'{name}_munge_inputs.txt', separator=' ')
        munge_rows = munge.height
    return {'traitname': name, 'vcf_files': str(vcf), 'gpca_rows': gpca.height,
            'excluded_non_autosomal_rows': summary['excluded_non_autosomal_rows'], 'munge_rows': munge_rows,
            'p_underflow_rows': int((frame['p'] == 0).sum()), 'success': True, 'error': ''}


def prepare_inputs(manifest, outdir, splitby_chr='split', gpca_id_source='chr_pos_ref_alt',
                   write_munge_inputs=False, hapmap_file=None, munge_id_source='vcf_id',
                   prepare_workers=4, bcftools='bcftools', p_min=1e-300):
    if not math.isfinite(p_min) or not 0 < p_min < 1:
        raise ValueError('--p-min must be finite and strictly between 0 and 1')
    if splitby_chr not in ('split','nosplit') or gpca_id_source not in ID_CHOICES or munge_id_source not in ID_CHOICES:
        raise ValueError('Invalid layout or identifier source')
    if not isinstance(prepare_workers, int) or prepare_workers < 1:
        raise ValueError('--prepare-workers must be a positive integer')
    if write_munge_inputs != bool(hapmap_file):
        raise ValueError('--write-munge-inputs and --hapmap-file must be supplied together')
    executable = shutil.which(bcftools)
    if not executable:
        raise ValueError(f'Local bcftools executable not found: {bcftools}')
    rows = manifest_inputs(manifest)
    hm3 = None
    if write_munge_inputs:
        hm3 = pl.read_csv(hapmap_file, separator='\t', columns=['SNP'], schema_overrides={'SNP': pl.String}).drop_nulls().unique()
        if hm3.is_empty():
            raise ValueError('HapMap SNP list is empty')
    outdir = Path(outdir).resolve()
    reports = ['GPCA_Input_QC_Issues.csv', 'GPCA_Input_QC_Summary.csv']
    outputs = ['gpca_inputs', 'Preparation_Status.csv', 'Preparation_Settings.json', *reports]
    if write_munge_inputs:
        outputs.append('munge_inputs')
    for name in outputs:
        path = outdir/name
        if path.exists() or path.is_symlink():
            raise FileExistsError(f'Refusing to overwrite {path}; use a fresh preparation outdir or pass existing --gpca_input_folder to gpca')
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=outdir) as temporary:
        stage = Path(temporary)
        (stage/'gpca_inputs').mkdir()
        (stage/'qc').mkdir()
        if write_munge_inputs:
            (stage/'munge_inputs').mkdir()
        statuses = {}
        with ThreadPoolExecutor(max_workers=prepare_workers) as executor:
            futures = {executor.submit(prepare_trait, name, vcf, stage, executable, splitby_chr,
                                       gpca_id_source, munge_id_source, hm3, p_min): (name, vcf) for name,vcf in rows}
            for future in as_completed(futures):
                name, vcf = futures[future]
                try:
                    statuses[name] = future.result()
                    print(f'Prepared {name}: {statuses[name]["gpca_rows"]} GPCA variants', flush=True)
                except Exception as error:
                    detail = getattr(error, 'stderr', '') or str(error)
                    statuses[name] = {'traitname':name, 'vcf_files':str(vcf), 'gpca_rows':0,
                                      'excluded_non_autosomal_rows':0, 'munge_rows':0, 'p_underflow_rows':0,
                                      'success':False, 'error':str(detail)}
        ordered = [statuses[name] for name,_ in rows]
        summaries = []
        for status in ordered:
            path = stage/'qc'/f'{status["traitname"]}.summary.csv'
            summary = pl.read_csv(path, schema_overrides={'traitname':pl.String}).to_dicts()[0] if path.exists() else {
                'traitname':status['traitname'], 'input_rows':None, 'retained_rows':None,
                'removed_rows':None, 'p_adjusted_rows':None, 'excluded_non_autosomal_rows':None}
            summaries.append({**summary, 'success':status['success'], 'error':status['error']})
        pl.DataFrame(summaries, infer_schema_length=None).write_csv(stage/reports[1])
        merge_issue_reports([stage/'qc'/f'{name}.csv' for name,_ in rows
                             if (stage/'qc'/f'{name}.csv').exists()], stage/reports[0])
        failed = [s for s in ordered if not s['success']]
        if failed:
            pl.DataFrame(ordered).write_csv(outdir/'Preparation_Status.csv')
            for report in reports:
                os.rename(stage/report, outdir/report)
            raise ValueError('VCF preparation failed; no prepared tables published:\n' + '\n'.join(f'{s["traitname"]}: {s["error"]}' for s in failed))
        pl.DataFrame(ordered).write_csv(stage/'Preparation_Status.csv')
        settings = {'manifest':str(Path(manifest).resolve()), 'splitby_chr':splitby_chr,
                    'gpca_id_source':gpca_id_source, 'munge_id_source':munge_id_source,
                    'write_munge_inputs':write_munge_inputs, 'hapmap_file':str(Path(hapmap_file).resolve()) if hapmap_file else None,
                    'N_source':'FORMAT/NEF', 'Z_source':'FORMAT/ES / FORMAT/SE',
                    'P_source':'10 ** (-FORMAT/LP), floored at p_min for valid LP', 'p_min':p_min,
                    'duplicate_policy':'Largest valid LP, first original row on ties',
                    'QC_issues_source':'Original VCF fields, unchanged; union of headers across samples',
                    'bcftools':executable}
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
    print(f'GPCA inputs: {outdir / "gpca_inputs"}')
    return outdir/'gpca_inputs'


def preparation_kwargs(opts):
    return {name: getattr(opts,name) for name in ['gpca_id_source','write_munge_inputs','hapmap_file',
                                                'munge_id_source','prepare_workers','bcftools','p_min']}


def main(argv=None):
    parser = HelpParser(prog='ldsc-gpca prepare', description=__doc__, epilog=PREPARE_INPUT_HELP)
    parser.add_argument('--input', required=True, help='CSV with traitname and vcf_files. Relative VCF paths resolve beside this manifest.')
    parser.add_argument('--outdir', required=True, help='Output directory; creates gpca_inputs and optionally munge_inputs.')
    parser.add_argument('--splitby_chr', choices=['split','nosplit'], default='split', help='Default: split (requires all chromosomes 1–22 per trait).')
    add_prepare_options(parser)
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        parser.print_help()
        return 0
    opts = parser.parse_args(argv)
    try:
        prepare_inputs(opts.input, opts.outdir, splitby_chr=opts.splitby_chr, **preparation_kwargs(opts))
    except (OSError, ValueError, pl.exceptions.PolarsError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0
