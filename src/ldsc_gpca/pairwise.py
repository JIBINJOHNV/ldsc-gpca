"""Pairwise LDSC batching and bounded execution retries."""
import os
import shlex
import math
import concurrent.futures
from .utils import run_command
from .ldsc_runtime import ldsc_regression_command
from .ldsc_export import RESULT_SUFFIX

def parallel_ldsc_analysis(n_parallel, batch_size, output_path, ld_ref_dir, input_df, ldsc_input_path, retries=1, runtime=None):
    if not isinstance(retries, int) or retries < 0:
        raise ValueError('LDSC retries must be a non-negative integer')
    os.makedirs(output_path, exist_ok=True)

    def ldsc_analysis(target_files, reference_sumstat, part, s_prevalence, p_prevalence):
        ref_prefix = os.path.basename(reference_sumstat).replace('.sumstats.gz', '')
        out_prefix = os.path.join(output_path, f"{ref_prefix}_{batch_size}_ldsc_{part}")

        prevalence_flags = []
        if any(value != 'nan' for value in p_prevalence.split(',')):
            prevalence_flags = ['--samp-prev', s_prevalence, '--pop-prev', p_prevalence]
        command = shlex.join(ldsc_regression_command(**(runtime or {})) + [
            '--rg', f'{reference_sumstat},{target_files}', '--ref-ld-chr', ld_ref_dir,
            '--w-ld-chr', ld_ref_dir, *prevalence_flags, '--out', out_prefix])
        for attempt in range(1, retries + 2):
            label = f'LDSC_{ref_prefix}_batch_{part} attempt {attempt}/{retries + 1}'
            try:
                run_command(command, label, output_folder=os.path.dirname(os.path.abspath(output_path)))
                return out_prefix + RESULT_SUFFIX
            except RuntimeError as error:
                if attempt == retries + 1:
                    raise RuntimeError(f'{label}: retries exhausted; see execution_errors.log') from error
                print(f'  [!] {label} failed; retrying this batch only.', flush=True)

    ref_names = input_df[input_df['ref'] == 'yes']['gwas_name'].unique()
    num_batches = math.ceil(len(input_df) / batch_size)
    print(f"  -> Total Batches to process: {len(ref_names) * num_batches}")

    with concurrent.futures.ThreadPoolExecutor(n_parallel) as executor:
        futures = []
        for ref_name in ref_names:
            ref_row = input_df[input_df['gwas_name'] == ref_name].iloc[0]
            ref_file = os.path.join(ldsc_input_path, f"{ref_name}.sumstats.gz")
            for i in range(0, len(input_df), batch_size):
                batch = input_df.iloc[i : i + batch_size]
                t_files = ",".join([os.path.join(ldsc_input_path, f"{n}.sumstats.gz") for n in batch['gwas_name']])
                s_prev = f"{ref_row['sample_prevalence']}," + ",".join(batch['sample_prevalence'].astype(str))
                p_prev = f"{ref_row['pop_prevalence']}," + ",".join(batch['pop_prevalence'].astype(str))
                futures.append(executor.submit(ldsc_analysis, t_files, ref_file, i//batch_size, s_prev, p_prev))
        return [future.result() for future in futures]
