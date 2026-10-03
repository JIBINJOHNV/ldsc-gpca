"""Pairwise LDSC batching, verified restarts and bounded execution retries."""
import os
import shlex
import math
import concurrent.futures
import pandas as pd
from .utils import run_command
from .ldsc_runtime import ldsc_regression_command
from .ldsc_export import RESULT_SUFFIX, write_results_csv
from .restart import (RestartContext, CHECKPOINT_SUFFIX, file_digest,
                      reusable_result, validate_batch, write_checkpoint)


class _BatchFailure(RuntimeError):
    def __init__(self, message, can_exclude, attempts):
        super().__init__(message)
        self.can_exclude = can_exclude
        self.attempts = attempts


class LDSCBatchFailures(RuntimeError):
    """Fatal execution failures, with completed paths retained for diagnostics/restart."""
    def __init__(self, message, result_files, failed_pairs, can_exclude):
        super().__init__(message)
        self.result_files = result_files
        self.failed_pairs = failed_pairs
        self.can_exclude = can_exclude


def _result_stamp(path):
    try:
        stat = os.stat(path)
        return (stat.st_ino, stat.st_mtime_ns, stat.st_size) if stat.st_size else None
    except FileNotFoundError:
        return None

def parallel_ldsc_analysis(n_parallel, batch_size, output_path, ld_ref_dir, input_df, ldsc_input_path, retries=1, runtime=None, ld_weights_dir=None, *, restart=False, checkpoint_parameters=None):
    if not isinstance(retries, int) or retries < 0:
        raise ValueError('LDSC retries must be a non-negative integer')
    os.makedirs(output_path, exist_ok=True)
    context = None
    if restart or checkpoint_parameters is not None:
        print('  -> Verifying LDSC runtime and input contents for restart checkpoints...', flush=True)
        context = RestartContext(ld_ref_dir, ld_weights_dir or ld_ref_dir,
            [os.path.join(ldsc_input_path, f'{name}.sumstats.gz') for name in input_df.gwas_name],
            runtime or {}, checkpoint_parameters or {})

    def ldsc_analysis(target_files, reference_sumstat, part, s_prevalence, p_prevalence):
        ref_prefix = os.path.basename(reference_sumstat).replace('.sumstats.gz', '')
        out_prefix = os.path.join(output_path, f"{ref_prefix}_{batch_size}_ldsc_{part}")

        prevalence_flags = []
        if any(value != 'nan' for value in p_prevalence.split(',')):
            prevalence_flags = ['--samp-prev', s_prevalence, '--pop-prev', p_prevalence]
        command_args = ldsc_regression_command(**(runtime or {})) + [
            '--rg', f'{reference_sumstat},{target_files}', '--ref-ld-chr', ld_ref_dir,
            '--w-ld-chr', ld_weights_dir if ld_weights_dir is not None else ld_ref_dir, *prevalence_flags, '--out', out_prefix]
        command = shlex.join(command_args)
        result_path = out_prefix + RESULT_SUFFIX
        targets = target_files.split(',')
        sumstats_paths = [reference_sumstat, *targets]
        request = context.request(command_args, sumstats_paths) if context else None
        reason = 'restart_not_requested'
        if restart:
            reuse, reason = reusable_result(result_path, request, reference_sumstat, targets)
            if reuse:
                context.check_unchanged(sumstats_paths)
                print(f'  -> Reusing completed batch: {ref_prefix}, batch {part}', flush=True)
                return result_path, 'reused', reason, 0
        # Invalidate completion before launching; an interruption must not leave
        # an older checkpoint endorsing this attempt's result.
        checkpoint_path = result_path + CHECKPOINT_SUFFIX
        if os.path.exists(checkpoint_path):
            os.unlink(checkpoint_path)
        for attempt in range(1, retries + 2):
            label = f'LDSC_{ref_prefix}_batch_{part} attempt {attempt}/{retries + 1}'
            stage = 'execution'
            try:
                previous = _result_stamp(result_path)
                run_command(command, label, output_folder=os.path.dirname(os.path.abspath(output_path)))
                current = _result_stamp(result_path)
                if current is None or current == previous:
                    raise RuntimeError(f'{label}: no fresh non-empty numerical CSV was produced')
                if context:
                    stage = 'validation'
                    validate_batch(result_path, reference_sumstat, targets)
                    context.check_unchanged(sumstats_paths)
                    write_checkpoint(checkpoint_path, {'status': 'completed', 'request': request,
                                                       'result_sha256': file_digest(result_path)})
                return result_path, 'completed', reason, attempt
            except (RuntimeError, OSError, ValueError) as error:
                if attempt == retries + 1:
                    raise _BatchFailure(f'{label}: retries exhausted: {error}; see execution_errors.log',
                                        can_exclude=stage == 'execution', attempts=attempt) from error
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
                future = executor.submit(ldsc_analysis, t_files, ref_file, i//batch_size, s_prev, p_prev)
                futures.append((future, ref_name, i//batch_size))
        paths, errors, failed_pairs = {}, [], {}
        can_exclude = True
        statuses = [dict(Reference=ref, Batch=part, Status='pending', Result_File='', Error='', Restart_Reason='', Attempts=0)
                    for _, ref, part in futures]
        status_path = os.path.join(output_path, 'LDSC_Batch_Status.csv')
        write_results_csv(pd.DataFrame(statuses), status_path)
        indices = {future: i for i, (future, _, _) in enumerate(futures)}
        for future in concurrent.futures.as_completed(indices):
            index = indices[future]
            try:
                path, status, reason, attempts = future.result()
                paths[index] = path
                statuses[index].update(Status=status, Result_File=path, Restart_Reason=reason, Attempts=attempts)
            except (RuntimeError, OSError, ValueError) as error:
                errors.append(str(error))
                statuses[index].update(Status='execution_failed', Error=str(error), Attempts=getattr(error, 'attempts', 0))
                can_exclude = can_exclude and isinstance(error, _BatchFailure) and error.can_exclude
                _, reference, part = futures[index]
                for target in input_df.iloc[part * batch_size:(part + 1) * batch_size].gwas_name:
                    failed_pairs[(reference, target)] = str(error)
            write_results_csv(pd.DataFrame(statuses), status_path)
        if errors:
            raise LDSCBatchFailures('Analysis stopped: LDSC execution failed; completed batch files were retained. '
                'See LDSC_Batch_Status.csv.\n' + '\n'.join(errors),
                [paths[i] for i in sorted(paths)], failed_pairs, can_exclude)
        return [paths[i] for i in sorted(paths)]
