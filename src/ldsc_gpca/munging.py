"""LDSC summary-statistic munging and provenance."""
import os
import shlex
import json
import hashlib
import concurrent.futures
import pandas as pd
from .utils import DEFAULT_FILTERS, is_valid_gz, run_command
from .ldsc_runtime import ldsc_command

def parallel_munge_sumstats(n_parallel, input_df, output_folder, snp_include_file, filters=None, *, ldsc_input_folder=None, munge_input_folder=None, runtime=None):
    ldsc_input_folder = ldsc_input_folder or os.path.join(output_folder, "ldsc_input")
    munge_input_folder = munge_input_folder or os.path.join(output_folder, "munge_input")
    filters = {**DEFAULT_FILTERS, **(filters or {})}
    os.makedirs(ldsc_input_folder, exist_ok=True)

    def munge_sumstats(row):
        sample_name = row['gwas_name']
        n_column = row['sample_size_column']
        in_file = os.path.join(munge_input_folder, f"{sample_name}_mungeinput.tsv")
        out_prefix = os.path.join(ldsc_input_folder, sample_name)
        final_out_file = f"{out_prefix}.sumstats.gz"

        # Re-munge full runs: an existing file may have used a different N convention.
        metadata_file = out_prefix + '.prevalence.json'
        if os.path.exists(metadata_file):
            os.remove(metadata_file)
        ignore = ['--ignore', 'N_CASES,N_CONTROLS,NEF'] if n_column == 'N_TOTAL' else []
        command = shlex.join(ldsc_command('munge_sumstats.py', **(runtime or {})) + [
            '--sumstats', in_file, '--N-col', n_column, *ignore, '--snp', 'ID', '--a1', 'ALT',
            '--a2', 'REF', '--p', 'P', '--frq', 'AF', '--maf-min', str(filters['munge_maf_min']),
            '--signed-sumstats', 'EZ,0', '--merge-alleles', snp_include_file, '--out', out_prefix])

        try:
            run_command(command, f"Munge_{sample_name}", output_folder=output_folder)
            if not is_valid_gz(final_out_file):
                raise RuntimeError(f'{sample_name}: munging produced no valid gzip output')
            with open(final_out_file, 'rb') as handle:
                digest = hashlib.sha256()
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            with open(metadata_file, 'w') as handle:
                json.dump({'sample_size_column': n_column, 'sha256': digest.hexdigest(), 'filters': filters,
                           'sample_prevalence': None if pd.isna(row['sample_prevalence']) else float(row['sample_prevalence'])}, handle)
            return sample_name, True
        except RuntimeError:
            print(f"  [!] Failed: {sample_name} (See execution_errors.log)")
            return sample_name, False

    success_list = []
    with concurrent.futures.ThreadPoolExecutor(n_parallel) as executor:
        futures = {executor.submit(munge_sumstats, row): row['gwas_name'] for _, row in input_df.iterrows()}
        for future in concurrent.futures.as_completed(futures):
            name, success = future.result()
            if success: success_list.append(name)
    return success_list
