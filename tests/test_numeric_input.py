"""Real R numeric parsing and missing-value policy regressions."""
import os
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RBoundaryTests(unittest.TestCase):
    def test_numeric_parsing(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        run = subprocess.run([executable, str(ROOT/'tests/test_numeric_input.R'), str(ROOT)],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn('Early numeric parsing and genuine missing-value policy tests passed', run.stdout)
