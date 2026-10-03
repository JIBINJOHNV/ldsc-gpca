"""Cleanup must preserve identifier, QC, ordering and missing-estimate semantics."""
import math
import os
from pathlib import Path
import shutil
import subprocess
import unittest

import pandas as pd

from ldsc_gpca.result_qc import trait_status
from ldsc_gpca.trait_selection import _self_estimates


class CleanupTests(unittest.TestCase):
    def test_indexed_summaries_preserve_duplicates_missing_values_and_order(self):
        frame = pd.DataFrame({'p1': ['B', 'A', 'A', 'B'], 'p2': ['B', 'A', 'A', 'A'],
            'h2_obs': [.4, math.nan, .2, .9], 'h2_obs_se': [.04, .02, .02, .1],
            'h2_liab': [math.nan]*4, 'h2_liab_se': [math.nan]*4,
            'h2_int': [1.1, 1., 1.2, 1.], 'h2_int_se': [.01]*4})
        pairs = pd.DataFrame({'p1': ['B', 'A', 'A', 'B', 'A'],
            'p2': ['B', 'A', 'A', 'A', 'B'],
            'Status': ['valid', 'failed_estimate', 'failed_estimate', 'failed_estimate', 'valid'],
            'Reason': ['', 'missing h2', 'missing h2', 'bad pair', '']})
        original = frame.copy(deep=True)
        actual = trait_status(frame, pairs, ['C', 'A', 'B'])
        self.assertEqual(actual.Trait.tolist(), ['C', 'A', 'B'])
        self.assertEqual(actual.Self_Status.tolist(), ['not_requested', 'failed_estimate', 'valid'])
        self.assertEqual(actual.Self_Reason.tolist(), ['', 'missing h2', ''])
        self.assertEqual(actual.Failed_Pair_Rows.tolist(), [0, 3, 1])
        self.assertTrue(math.isnan(actual.h2_obs.iloc[1]))  # First self row, not first nonmissing.
        self.assertEqual(actual.h2_int.tolist()[1:], [1., 1.1])
        estimates = _self_estimates(frame, ['C', 'A', 'B'])
        self.assertTrue(math.isnan(estimates['C'][0]) and math.isnan(estimates['A'][0]))
        self.assertEqual(estimates['B'], (.4, .04, 10.))
        pd.testing.assert_frame_equal(frame, original, check_exact=True)
        empty = trait_status(frame.iloc[:0], pairs.iloc[:0], ['C'])
        self.assertEqual(empty.Self_Status.iloc[0], 'not_requested')
        self.assertEqual(empty.Failed_Pair_Rows.iloc[0], 0)
        self.assertTrue(math.isnan(empty.h2_obs.iloc[0]))

    def test_r_identifiers_qc_and_pc1_reports(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run([executable, str(root/'tests/test_cleanup.R'), str(root)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
