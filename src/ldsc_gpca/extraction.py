"""VCF filters, extraction and prevalence preparation."""
import os
import subprocess
import shlex
import concurrent.futures
import pandas as pd
import polars as pl
from .utils import DEFAULT_FILTERS, optional_prevalence

def filter_commands(vcf_file, filters, query):
    expression = "(ABS(INFO/AF - INFO/EUR) > {d} || INFO/EUR==\".\")".format(d=filters['max_af_difference'])
    maf = filters['maf_min']
    first = ['/usr/bin/bcftools', 'view', vcf_file]
    if filters['exclude_mhc']:
        first[2:2] = ['-t', f"^{filters['mhc_chr']}:{filters['mhc_start']}-{filters['mhc_end']}"]
    commands = [
        first,
        ['/usr/bin/bcftools', 'view', '-e', expression],
        ['/usr/bin/bcftools', 'view', '-i', f'FORMAT/SI >= {filters["info_min"]} && FORMAT/AF >= {maf} && FORMAT/AF[0:0] <= {1-maf}'],
    ]
    if filters['remove_palindrome']:
        lo, hi = filters['pal_lower'], filters['pal_upper']
        pal = '((REF=="A" && ALT=="T") || (REF=="T" && ALT=="A") || (REF=="C" && ALT=="G") || (REF=="G" && ALT=="C"))'
        commands.append(['/usr/bin/bcftools', 'view', '-e', f'{pal} && FORMAT/AF[0:0] >= {lo} && FORMAT/AF[0:0] <= {hi}'])
    commands.append(['/usr/bin/bcftools', 'query', '-f', query])
    return commands

def munge_input_worker(sample_name, input_path, output_folder, vcf_file, filters=None):
    munge_input_folder = os.path.join(output_folder, 'munge_input')
    os.makedirs(munge_input_folder, exist_ok=True)
    out_file = os.path.join(munge_input_folder, f"{sample_name}_mungeinput.tsv")
    filters = filters or DEFAULT_FILTERS
    commands = filter_commands(vcf_file, filters, r'%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%EZ\t%LP\t%AF\t%NEF]\n')
    commands.append(['awk', '-F', '\t', '-v', 'OFS= ', r'{ pval = ($7 == "." || $7 == "") ? "NA" : 10^(-$7); print $1, $2, $3, $4, $5, $6, pval, $8, $9 }'])
    docker_command = ['docker', 'run', '--rm', '-v', f'{output_folder}:{output_folder}',
                      '-v', f'{input_path}:{input_path}', '--user', f'{os.getuid()}:{os.getgid()}',
                      'jibinjv/ldsc:v3', 'bash', '-o', 'pipefail', '-c',
                      ' | '.join(shlex.join(command) for command in commands)]
    try:
        with open(out_file, 'w') as handle:
            handle.write('CHROM ID POS REF ALT EZ P AF NEF\n')
            handle.flush()
            subprocess.run(docker_command, check=True, stdout=handle, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command for {sample_name}: {e.stderr.decode().strip()}")
        raise e

def munge_input_worker_with_prevalence(sample_name, input_path, output_folder, vcf_file, sample_prevalence, filters=None):
    """Create TSV, add total N with Polars, and return sample prevalence."""
    folder = os.path.join(output_folder, 'munge_input')
    os.makedirs(folder, exist_ok=True)
    out_file = os.path.join(folder, f'{sample_name}_mungeinput.tsv')
    query = r'%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%EZ\t%LP\t%AF\t%NEF\t%NC\t%NCO]\n'
    filters = filters or DEFAULT_FILTERS
    commands = filter_commands(vcf_file, filters, query)
    commands.append(['awk', '-F', '\t', '-v', 'OFS=\t', r'{ $7 = ($7 == "." || $7 == "") ? "NA" : sprintf("%.17g", 10^(-$7)); print }'])
    command = ['docker', 'run', '--rm', '-v', f'{output_folder}:{output_folder}',
               '-v', f'{input_path}:{input_path}', '--user', f'{os.getuid()}:{os.getgid()}',
               'jibinjv/ldsc:v3', 'bash', '-o', 'pipefail', '-c', ' | '.join(shlex.join(p) for p in commands)]
    columns = ['CHROM', 'ID', 'POS', 'REF', 'ALT', 'EZ', 'P', 'AF', 'NEF', 'N_CASES', 'N_CONTROLS']
    with open(out_file, 'w') as handle:
        handle.write('\t'.join(columns) + '\n')
        handle.flush()
        subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=True)
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
def run_vcf_to_table(n_parallel, input_df, output_folder, filters=DEFAULT_FILTERS):
    os.makedirs(output_folder, exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {}
        for _, row in input_df.iterrows():
            sample_name = row['gwas_name']
            vcf_file = row['vcf_files']
            input_path = "/".join(vcf_file.split("/")[:-1])
            # Submit to the un-nested top-level function
            if pd.isna(row['pop_prevalence']):
                future = executor.submit(munge_input_worker, sample_name, input_path, output_folder, vcf_file, filters)
            else:
                future = executor.submit(munge_input_worker_with_prevalence, sample_name, input_path,
                                         output_folder, vcf_file, row['sample_prevalence'], filters)
            futures[future] = row.name
        # CRITICAL: Checking results prints errors to your screen if a process crashes!
        for future in concurrent.futures.as_completed(futures):
            prevalence = future.result()
            if prevalence is not None:
                index = futures[future]
                input_df.loc[index, 'sample_prevalence_source'] = 'median_case_fraction' if pd.isna(input_df.loc[index, 'sample_prevalence']) else 'provided'
                input_df.loc[index, 'sample_prevalence'] = prevalence
    input_df.to_csv(os.path.join(output_folder, 'LDSC_Trait_Prevalence_Metadata.csv'), index=False)
    return input_df
