"""Exercise literal trait identifiers at the real R input boundaries."""
import os
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]


class LiteralTraitIDTests(unittest.TestCase):
    def run_r(self, group):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        result = subprocess.run(
            [executable, str(ROOT / 'tests/test_literal_trait_ids.R'), str(ROOT), group],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f'{group} literal trait ID checks passed', result.stdout)

    def test_manifests_preserve_names_and_reject_empty_identifiers(self):
        self.run_r('manifest')

    def test_ldsc_reading_numeric_missing_values_and_matrix_parity(self):
        self.run_r('ldsc')

    def test_genomicsem_estimator_receives_exact_identifiers(self):
        self.run_r('genomicsem')


if __name__ == '__main__':
    unittest.main()
