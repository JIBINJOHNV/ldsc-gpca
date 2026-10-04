"""Convert GWAS VCFs to native munging tables without changing GenomicSEM."""
from pathlib import Path
import tempfile

import polars as pl

from .prepare import (QUERY, RAW_COLUMNS, extract_table, validate_and_transform,
                      write_original_issues, merge_issue_reports, validate_ids)
from .workers import run_parallel_jobs


def prepare_trait(row, stage, executable, p_min):
    name, vcf = row['traitname'], Path(row['source_file'])
    n = row['N']
    query = QUERY.replace(r'\t%NEF', r'\t%NEF\t%SI' if n == 'NA' else f'\t{n:.17g}' + r'\t%SI')
    with tempfile.NamedTemporaryFile(dir=stage, suffix='.tsv') as temporary:
        raw = extract_table(vcf, executable, temporary.name, query=query, columns=[*RAW_COLUMNS, 'INFO'])
        frame, issues, summary = validate_and_transform(raw, p_min, 'vcf_id', require_info=True)
    write_original_issues(vcf, stage/'qc'/f'{name}.csv', name, issues)
    summary = {'traitname': name, 'vcf_files': str(vcf), **summary,
               'N_source': 'FORMAT/NEF' if n == 'NA' else 'manifest_N',
               'N_override': None if n == 'NA' else n, 'p_min': p_min,
               'INFO_source': 'FORMAT/SI'}
    pl.DataFrame([summary]).write_csv(stage/'qc'/f'{name}.summary.csv')
    if frame.is_empty():
        raise ValueError(f'{name}: no usable VCF records remain; see GenomicSEM_VCF_QC_Issues.csv')
    validate_ids(frame, 'SNP')
    frame.select('SNP', 'CHR', 'POS', 'A1', 'A2', pl.min_horizontal('eaf_A1', 1-pl.col('eaf_A1')).alias('MAF'),
                 'beta', 'se', 'N', 'p', 'INFO').write_csv(stage/'vcf_input'/f'{name}_munge_inputs.tsv', separator='\t')
    return summary


def prepare_vcf_inputs(rows, outdir, executable, workers, p_min):
    """Publish tables only after all traits succeed; retain QC when a job fails."""
    outdir = Path(outdir)
    with tempfile.TemporaryDirectory(prefix='.native-vcf-', dir=outdir) as temporary:
        stage = Path(temporary)
        (stage/'vcf_input').mkdir()
        (stage/'qc').mkdir()
        values = run_parallel_jobs([(row['traitname'], prepare_trait, (row, stage, executable, p_min))
            for row in rows], workers, stage='GenomicSEM VCF extraction', collect_failures=True,
            status_path=outdir/'GenomicSEM_VCF_Worker_Attempts.csv')
        summaries, failures = [], []
        for row, value in zip(rows, values):
            path = stage/'qc'/f'{row["traitname"]}.summary.csv'
            summary = pl.read_csv(path, schema_overrides={'traitname': pl.String}).to_dicts()[0] if path.exists() else {
                'traitname': row['traitname'], 'vcf_files': row['source_file']}
            failed = isinstance(value, Exception)
            error = str(getattr(value, 'stderr', '') or value or type(value).__name__) if failed else ''
            summaries.append({**summary, 'success': not failed, 'error': error})
            if failed:
                failures.append(f'{row["traitname"]}: {error}')
        pl.DataFrame(summaries, infer_schema_length=None).write_csv(outdir/'GenomicSEM_VCF_QC_Summary.csv')
        merge_issue_reports([stage/'qc'/f'{row["traitname"]}.csv' for row in rows
            if (stage/'qc'/f'{row["traitname"]}.csv').exists()], outdir/'GenomicSEM_VCF_QC_Issues.csv')
        if failures:
            raise ValueError('VCF extraction failed after 2 attempts; native LDSC was not started. '
                             'See GenomicSEM_VCF_Worker_Attempts.csv.\n' + '\n'.join(failures))
        (stage/'vcf_input').rename(outdir/'vcf_input')
    return [{**row, 'source_file': str(outdir/'vcf_input'/f'{row["traitname"]}_munge_inputs.tsv')} for row in rows]
