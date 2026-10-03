"""Two-attempt execution for package-managed parallel jobs."""
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from .ldsc_export import write_results_csv


def run_parallel_jobs(jobs, workers, *, stage, status_path=None, processes=False,
                      collect_failures=False):
    """Run (label, function, args) jobs; retry only failures once in a fresh pool."""
    if not isinstance(workers, int) or isinstance(workers, bool) or workers < 1:
        raise ValueError('Parallel worker count must be a positive integer')
    jobs = list(jobs)
    values, failed, history = [None] * len(jobs), list(range(len(jobs))), []
    executor_type = ProcessPoolExecutor if processes else ThreadPoolExecutor
    if status_path is not None:
        Path(status_path).parent.mkdir(parents=True, exist_ok=True)

    for attempt in (1, 2):
        if not failed:
            break
        pending, failed, recorded = failed, [], set()

        def record(index, value, error=None):
            values[index] = value if error is None else error
            recorded.add(index)
            if error is not None:
                failed.append(index)
            detail = '' if error is None else str(getattr(error, 'stderr', '') or error)
            history.append(dict(Job=jobs[index][0], Attempt=attempt,
                                Success=error is None, Error=detail))

        count = min(workers, len(pending)) if attempt == 1 else 1
        if attempt == 2:
            print(f'[{stage}] Retrying {len(pending)} failed job(s), one at a time (attempt 2/2).', flush=True)
        try:
            with executor_type(max_workers=count) as executor:
                futures = {executor.submit(jobs[i][1], *jobs[i][2]): i for i in pending}
                for future in as_completed(futures):
                    index = futures[future]
                    try:
                        value = future.result()
                    except Exception as error:
                        record(index, None, error)
                    else:
                        record(index, value)
        except Exception as error:
            # A broken process pool may fail during submission as well as collection.
            for index in pending:
                if index not in recorded:
                    record(index, None, error)
        if status_path is not None:
            write_results_csv(pd.DataFrame(history), status_path)
    if failed and not collect_failures:
        details = '\n'.join(f'{jobs[i][0]}: {values[i]}' for i in sorted(failed))
        audit = f' See {status_path}.' if status_path is not None else ''
        raise ValueError(f'Analysis stopped: {stage} jobs failed after 2 attempts.{audit}\n{details}')
    return values
