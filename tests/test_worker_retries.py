"""Failed jobs get two attempts; successful jobs and analysis order are retained."""
import contextlib
import gzip
import io
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import munging, workers

ROOT = Path(__file__).resolve().parents[1]


def exit_once(marker):
    path = Path(marker)
    if not path.exists():
        path.write_text('first worker exited')
        os._exit(17)
    return 'recovered in a fresh process'


class WorkerRetryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_jobs(self, jobs, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return workers.run_parallel_jobs(jobs, 2, stage='test jobs',
                status_path=self.root/'attempts.csv', **kwargs)

    def test_retry_only_failed_job_preserves_order_and_history(self):
        calls = {'A': 0, 'B': 0}
        def task(name):
            calls[name] += 1
            if name == 'A' and calls[name] == 1:
                raise RuntimeError('temporary interruption')
            return name.lower()
        result = self.run_jobs([(name, task, (name,)) for name in ('A', 'B')])
        self.assertEqual(result, ['a', 'b'])
        self.assertEqual(calls, {'A': 2, 'B': 1})
        history = pd.read_csv(self.root/'attempts.csv')
        self.assertEqual(history[history.Job == 'A'].Success.tolist(), [False, True])
        self.assertEqual(history[history.Job == 'A'].Attempt.tolist(), [1, 2])

    def test_permanent_failure_stops_after_two_attempts_and_names_job(self):
        calls = []
        def task():
            calls.append(1)
            raise OSError('disk unavailable')
        with self.assertRaisesRegex(ValueError, 'Analysis stopped: test jobs jobs failed after 2 attempts'):
            self.run_jobs([('chromosome 7', task, ())])
        self.assertEqual(len(calls), 2)
        history = pd.read_csv(self.root/'attempts.csv')
        self.assertEqual(history.Job.tolist(), ['chromosome 7'] * 2)
        self.assertTrue(history.Error.str.contains('disk unavailable').all())

    def test_real_exited_process_is_retried_in_new_pool(self):
        result = self.run_jobs([('crashed child', exit_once, (str(self.root/'marker'),))], processes=True)
        self.assertEqual(result, ['recovered in a fresh process'])
        self.assertEqual(pd.read_csv(self.root/'attempts.csv').Success.tolist(), [False, True])

    def test_failed_audit_write_is_fatal_and_not_silently_swallowed(self):
        with patch.object(workers, 'write_results_csv', side_effect=OSError('audit disk full')):
            with self.assertRaisesRegex(OSError, 'audit disk full'):
                self.run_jobs([('job', str, (1,))])

    def test_empty_jobs_and_invalid_worker_counts(self):
        self.assertEqual(self.run_jobs([]), [])
        for count in (0, -1, True, 1.5):
            with self.assertRaisesRegex(ValueError, 'positive integer'):
                workers.run_parallel_jobs([], count, stage='test')

    def test_munging_retries_and_never_drops_an_execution_failure(self):
        frame = pd.DataFrame({'gwas_name': ['A', 'B'], 'sample_size_column': ['NEF']*2,
                              'sample_prevalence': [float('nan')]*2})
        calls = {'A': 0, 'B': 0}
        permanent = False
        def command(cmd, *args, **kwargs):
            prefix = shlex.split(cmd)[-1]
            name = Path(prefix).name
            calls[name] += 1
            if name == 'B' and (permanent or calls[name] == 1):
                raise RuntimeError('munge worker interrupted')
            Path(prefix+'.sumstats.gz').write_bytes(gzip.compress(('SNP A1 A2 N Z\n' +
                ''.join(f'rs{i} A G 50000 {i/10}\n' for i in range(20))).encode()))
        with patch.object(munging, 'run_command', side_effect=command), contextlib.redirect_stdout(io.StringIO()):
            names = munging.parallel_munge_sumstats(2, frame, str(self.root), 'hm3')
            self.assertEqual(names, ['A', 'B'])
            self.assertEqual(calls, {'A': 1, 'B': 2})
            permanent = True; calls = {'A': 0, 'B': 0}
            with self.assertRaisesRegex(ValueError, 'LDSC munging jobs failed after 2 attempts'):
                munging.parallel_munge_sumstats(2, frame, str(self.root), 'hm3')
            self.assertEqual(calls, {'A': 1, 'B': 2})
            self.assertFalse((self.root/'ldsc_input/B.prevalence.json').exists())
            # An unchanged old gzip cannot turn a failed rerun into success.
            with patch.object(munging, 'run_command') as no_output:
                with self.assertRaisesRegex(ValueError, 'no fresh valid gzip'):
                    munging.parallel_munge_sumstats(1, frame.iloc[:1], str(self.root), 'hm3')
                self.assertEqual(no_output.call_count, 2)


class RWorkerRetryTests(unittest.TestCase):
    def test_workers_and_both_workflows(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        result = subprocess.run([executable, str(ROOT/'tests/test_worker_retries.R'), str(ROOT)],
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
