"""Actual bundled GWAMA numerical and output-QC regressions."""
import os
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RBoundaryTests(unittest.TestCase):
    def test_actual_gwama(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        run = subprocess.run([executable, str(ROOT/'tests/test_gwama_output_qc.R'), str(ROOT)],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn('zero-weight failure and audits passed', run.stdout)
