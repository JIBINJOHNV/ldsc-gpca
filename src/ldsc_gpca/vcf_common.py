"""Shared VCF schema, sample preflight, row QC and original-record audits.

Backend filters and preparation/publication orchestration remain with callers.
"""
import csv
import gzip
import math
from pathlib import Path
import subprocess

import polars as pl

ID_CHOICES = ('chr_pos_ref_alt', 'vcf_id')
RAW_COLUMNS = ['SNP', 'CHR', 'POS', 'A1', 'A2', 'eaf_A1', 'beta', 'se', 'LP', 'N']
QUERY = r'%ID\t%CHROM\t%POS\t%ALT\t%REF[\t%AF\t%ES\t%SE\t%LP\t%NEF]\n'



def require_single_sample(vcf, executable='bcftools'):
    """Reject ambiguous GWAS VCFs before any extraction output is opened."""
    samples = subprocess.run([executable, 'query', '-l', str(vcf)], check=True,
                             capture_output=True, text=True).stdout.splitlines()
    if len(samples) != 1:
        raise ValueError(f'{Path(vcf).name}: expected exactly one GWAS sample; found {len(samples)}. '
                         'Multi-sample and sample-free VCFs are not supported.')
    return samples[0]


def extract_table(vcf, executable, temporary, *, query=QUERY, columns=RAW_COLUMNS):
    require_single_sample(vcf, executable)
    with open(temporary, 'w') as handle:
        subprocess.run([executable, 'query', '-f', query, str(vcf)],
                       stdout=handle, stderr=subprocess.PIPE, text=True, check=True)
    if Path(temporary).stat().st_size == 0:
        raise ValueError(f'{vcf.name}: empty VCF query output')
    return pl.read_csv(temporary, separator='\t', has_header=False,
                       schema={name: pl.String for name in columns}, null_values='.')


def validate_and_transform(frame, p_min=1e-300, gpca_id_source='chr_pos_ref_alt', *, require_info=False):
    """Filter bad rows, retaining internal row indices solely for original-record audits."""
    if not math.isfinite(p_min) or not 0 < p_min < 1:
        raise ValueError('--p_min must be finite and strictly between 0 and 1')
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
    if require_info:
        frame = frame.with_columns(pl.col('INFO').cast(pl.Float64, strict=False))
        checks.append((pl.col('INFO').is_finite() & pl.col('INFO').is_between(0,1),
                       'INFO: missing, non-finite or outside [0,1]'))
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
