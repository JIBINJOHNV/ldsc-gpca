"""VCF filters, extraction and prevalence preparation."""
import os
import subprocess
import shlex
import tempfile
import pandas as pd
import polars as pl
from .utils import DEFAULT_FILTERS, optional_prevalence
from .workers import run_parallel_jobs

def filter_commands(vcf_file, filters, query, bcftools='bcftools'):
    expression = '(INFO/AF=="." || INFO/EUR=="." || ABS(INFO/AF - INFO/EUR) > {d})'.format(d=filters['max_af_difference'])
    maf = filters['maf_min']
    first = [bcftools, 'view', vcf_file]
    if filters['exclude_mhc']:
        first[2:2] = ['-t', f"^{filters['mhc_chr']}:{filters['mhc_start']}-{filters['mhc_end']}"]
    commands = [
        first,
        [bcftools, 'view', '-e', expression],
        [bcftools, 'view', '-i', f'FORMAT/SI >= {filters["info_min"]} && FORMAT/AF >= {maf} && FORMAT/AF[0:0] <= {1-maf}'],
    ]
    if filters['remove_palindrome']:
        lo, hi = filters['pal_lower'], filters['pal_upper']
        pal = '((REF=="A" && ALT=="T") || (REF=="T" && ALT=="A") || (REF=="C" && ALT=="G") || (REF=="G" && ALT=="C"))'
        commands.append([bcftools, 'view', '-e', f'{pal} && FORMAT/AF[0:0] >= {lo} && FORMAT/AF[0:0] <= {hi}'])
    commands.append([bcftools, 'query', '-f', query])
    return commands

def _extract_table(sample_name, out_file, commands, header):
    """Count missing INFO frequencies in the same stream; publish QC on success."""
    audit_file = out_file.removesuffix('_mungeinput.tsv') + '_AF_Filter_QC.csv'
    if os.path.exists(audit_file):
        os.unlink(audit_file)
    with tempfile.TemporaryDirectory(prefix='.af-qc-', dir=os.path.dirname(out_file)) as tmp:
        audit = os.path.join(tmp, 'counts.csv')
        # Inspect INFO only, including absent tags and missing elements of vectors.
        # Leave filtering to bcftools; count after optional MHC selection, before QC.
        counter = r'''
            function missing(tag) {
                return $8 !~ "(^|;)" tag "=" || $8 ~ "(^|;)" tag "=([^;]*,)?[.]([,;]|$)"
            }
            !/^#/ { n++; af=missing("AF"); eur=missing("EUR"); a+=af; e+=eur; any+=(af || eur) }
            { print }
            END {
                print "records_after_mhc,missing_info_af,missing_info_eur,excluded_missing_either" > audit
                printf "%.0f,%.0f,%.0f,%.0f\n", n+0, a+0, e+0, any+0 > audit
            }
        '''
        commands = [commands[0], ['awk', '-F', '\t', '-v', f'audit={audit}', counter], *commands[1:]]
        command = ['bash', '-o', 'pipefail', '-c', ' | '.join(shlex.join(p) for p in commands)]
        with open(out_file, 'w') as handle:
            handle.write(header + '\n')
            handle.flush()
            subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=True)
        counts = pl.read_csv(audit).row(0, named=True)
        os.replace(audit, audit_file)
    print(f'{sample_name}: excluded {counts["excluded_missing_either"]} records with missing INFO/AF or INFO/EUR '
          f'from {counts["records_after_mhc"]} records after MHC selection; AF QC: {audit_file}')

def munge_input_worker(sample_name, input_path, output_folder, vcf_file, filters=None, bcftools='bcftools'):
    munge_input_folder = os.path.join(output_folder, 'munge_input')
    os.makedirs(munge_input_folder, exist_ok=True)
    out_file = os.path.join(munge_input_folder, f"{sample_name}_mungeinput.tsv")
    filters = filters or DEFAULT_FILTERS
    commands = filter_commands(vcf_file, filters, r'%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%EZ\t%LP\t%AF\t%NEF]\n', bcftools)
    commands.append(['awk', '-F', '\t', '-v', 'OFS= ', r'{ pval = ($7 == "." || $7 == "") ? "NA" : 10^(-$7); print $1, $2, $3, $4, $5, $6, pval, $8, $9 }'])
    _extract_table(sample_name, out_file, commands, 'CHROM ID POS REF ALT EZ P AF NEF')

def munge_input_worker_with_prevalence(sample_name, input_path, output_folder, vcf_file, sample_prevalence, filters=None, bcftools='bcftools'):
    """Create TSV, add total N with Polars, and return sample prevalence."""
    folder = os.path.join(output_folder, 'munge_input')
    os.makedirs(folder, exist_ok=True)
    out_file = os.path.join(folder, f'{sample_name}_mungeinput.tsv')
    query = r'%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%EZ\t%LP\t%AF\t%NEF\t%NC\t%NCO]\n'
    filters = filters or DEFAULT_FILTERS
    commands = filter_commands(vcf_file, filters, query, bcftools)
    commands.append(['awk', '-F', '\t', '-v', 'OFS=\t', r'{ $7 = ($7 == "." || $7 == "") ? "NA" : sprintf("%.17g", 10^(-$7)); print }'])
    columns = ['CHROM', 'ID', 'POS', 'REF', 'ALT', 'EZ', 'P', 'AF', 'NEF', 'N_CASES', 'N_CONTROLS']
    _extract_table(sample_name, out_file, commands, '\t'.join(columns))
    frame = pl.read_csv(out_file, separator='\t', null_values=['.', 'NA', 'nan', 'NaN', ''],
                        schema={name: pl.String if i < 5 else pl.Float64 for i, name in enumerate(columns)})
    frame = frame.with_columns((pl.col('N_CASES') + pl.col('N_CONTROLS')).alias('N_TOTAL'))
    for name in ['N_CASES', 'N_CONTROLS', 'N_TOTAL']:
        if frame.is_empty() or not ((frame[name] > 0) & frame[name].is_finite()).fill_null(False).all():
            raise ValueError(f'{sample_name}: {name} must contain positive, finite counts for every variant')
    median_cases, median_controls, calculated = frame.select(pl.col('N_CASES').median(), pl.col('N_CONTROLS').median(), (pl.col('N_CASES') / pl.col('N_TOTAL')).median().alias('prevalence')).row(0)
    supplied = optional_prevalence(sample_prevalence, sample_name)
    resolved = optional_prevalence(calculated if supplied is None else supplied, sample_name)
    print(f'{sample_name}: median cases={median_cases:g}, controls={median_controls:g}; sample prevalence={resolved:g}')
    frame.write_csv(out_file, separator='\t', null_value='NA')
    return resolved

# --- MAIN LOGIC ---
def run_vcf_to_table(n_parallel, input_df, output_folder, filters=DEFAULT_FILTERS, bcftools='bcftools'):
    os.makedirs(output_folder, exist_ok=True)
    jobs = []
    for _, row in input_df.iterrows():
        name, vcf = row['gwas_name'], row['vcf_files']
        args = (name, os.path.dirname(vcf), output_folder, vcf)
        function = munge_input_worker
        if not pd.isna(row['pop_prevalence']):
            function = munge_input_worker_with_prevalence
            args += (row['sample_prevalence'],)
        jobs.append((name, function, (*args, filters, bcftools)))
    values = run_parallel_jobs(jobs, n_parallel, stage='VCF extraction', processes=True,
        status_path=os.path.join(output_folder, 'LDSC_Extraction_Worker_Attempts.csv'))
    for index, prevalence in zip(input_df.index, values):
        if prevalence is not None:
            input_df.loc[index, 'sample_prevalence_source'] = 'median_case_fraction' if pd.isna(input_df.loc[index, 'sample_prevalence']) else 'provided'
            input_df.loc[index, 'sample_prevalence'] = prevalence
    input_df.to_csv(os.path.join(output_folder, 'LDSC_Trait_Prevalence_Metadata.csv'), index=False)
    return input_df
