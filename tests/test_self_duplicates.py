"""Conflicting self estimates must fail before normalization or trait selection."""
import contextlib
import io
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import pandas as pd

from ldsc_gpca import results
from ldsc_gpca.ldsc_export import write_results_csv
from ldsc_gpca.trait_selection import validate_duplicate_estimates
from test_results_csv import record, manifest

ROOT = Path(__file__).resolve().parents[1]
TRAITS = ('C', 'A', 'B')
POLICIES = ('error', 'report', 'drop_traits')
SELF_FIELDS = ('h2_obs', 'h2_obs_se', 'h2_liab', 'h2_liab_se', 'h2_int', 'h2_int_se')


def fixture():
    rows = pd.DataFrame([record(a, b) for a in TRAITS for b in TRAITS])
    rows['h2_liab'] = .3
    rows['h2_liab_se'] = .04
    return rows


class SelfDuplicateTests(unittest.TestCase):
    def compile(self, root, rows, policy, extra=None):
        source = root / 'first.csv'
        write_results_csv(rows, source)
        paths = [source]
        if extra is not None:
            paths.append(root / 'second.csv')
            write_results_csv(extra, paths[-1])
        with contextlib.redirect_stdout(io.StringIO()):
            return results.compile_results(root, paths, manifest(TRAITS),
                                           result_failure_action=policy)

    def test_each_self_field_conflicts_in_every_policy_before_scale_mapping(self):
        # No prevalence is supplied: compilation would otherwise erase liability
        # columns. Both supplied scales must be checked before that happens.
        for policy in POLICIES:
            for field in SELF_FIELDS:
                for replacement in (.8, math.nan, math.inf, -math.inf):
                    with self.subTest(policy=policy, field=field, value=replacement), \
                            tempfile.TemporaryDirectory() as folder:
                        root = Path(folder)
                        rows = fixture()
                        extra = rows.iloc[:1].copy()
                        extra[field] = replacement
                        with self.assertRaisesRegex(RuntimeError, f'Conflicting duplicate.*{field}') as error:
                            self.compile(root, rows, policy, extra)
                        self.assertIn('first.csv:row 2', str(error.exception))
                        self.assertIn('second.csv:row 2', str(error.exception))
                        self.assertFalse((root / 'ldsc_results.csv').exists())
                        self.assertFalse((root / 'LDSC_Retained_Traits.csv').exists())
                        self.assertEqual(pd.read_csv(root / 'LDSC_Compilation_Status.csv').Status.iloc[0],
                                         'structural_failure')

    def test_identical_self_duplicates_and_empty_alternative_scale_preserve_results(self):
        for policy in POLICIES:
            for empty_alternative in (False, True):
                with self.subTest(policy=policy, empty_alternative=empty_alternative), \
                        tempfile.TemporaryDirectory() as folder:
                    root = Path(folder)
                    rows = fixture()
                    if empty_alternative:
                        rows[['h2_liab', 'h2_liab_se']] = math.nan
                    baseline = pd.read_csv(self.compile(root, rows, policy), float_precision='round_trip')
                    actual = pd.read_csv(self.compile(root, rows, policy, rows[rows.p1.eq(rows.p2)]),
                                         float_precision='round_trip')
                    pd.testing.assert_frame_equal(actual, baseline)
                    self.assertEqual(actual.p1.drop_duplicates().tolist(), list(TRAITS))

    def test_tolerance_is_absolute_and_inclusive_for_every_self_field(self):
        rows = fixture().iloc[:1]
        order = {rows.p1.iloc[0]: 0}
        for field in SELF_FIELDS:
            for difference, accepted in ((.001, True), (.00100001, False)):
                with self.subTest(field=field, difference=difference):
                    extra = rows.copy()
                    extra[field] += difference
                    combined = pd.concat([rows, extra], ignore_index=True)
                    if accepted:
                        validate_duplicate_estimates(combined, order)
                    else:
                        with self.assertRaisesRegex(RuntimeError, field):
                            validate_duplicate_estimates(combined, order)

    def test_positive_negative_mixture_conflicts_under_every_policy(self):
        for policy in POLICIES:
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as folder:
                rows = fixture()
                rows.loc[0, 'h2_obs'] = -.2
                extra = rows.iloc[:1].copy()
                extra['h2_obs'] = .8
                with self.assertRaisesRegex(RuntimeError, 'Conflicting duplicate.*h2_obs'):
                    self.compile(Path(folder), rows, policy, extra)

    def test_invalid_source_cannot_be_hidden_by_a_positive_mean_within_tolerance(self):
        for policy in POLICIES:
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                rows = fixture()
                rows.loc[0, 'h2_obs'] = -.0001
                extra = rows.iloc[:1].copy()
                extra['h2_obs'] = .0003
                if policy == 'error':
                    with self.assertRaises(results.EstimationFailure):
                        self.compile(root, rows, policy, extra)
                else:
                    output = self.compile(root, rows, policy, extra)
                    if policy == 'drop_traits':
                        self.assertEqual(pd.read_csv(output).p1.drop_duplicates().tolist(), ['A', 'B'])
                    else:
                        self.assertEqual(Path(output).name, 'ldsc_results_diagnostic.csv')
                        self.assertFalse((root / 'ldsc_results.csv').exists())
                status = pd.read_csv(root / 'LDSC_Pair_Status.csv')
                self.assertEqual(status.Status.eq('failed_estimate').sum(), 1)

    def test_reverse_off_diagonal_h2_and_intercepts_may_differ(self):
        rows = fixture()
        for field in SELF_FIELDS:
            rows.loc[1, field] = .8
        for policy in POLICIES:
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as folder:
                self.compile(Path(folder), rows, policy)

    def test_common_field_conflicts_are_fatal_in_every_policy(self):
        rows = fixture()
        rows.loc[1, 'rg'] += .2
        for policy in POLICIES:
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as folder:
                with self.assertRaisesRegex(RuntimeError, 'Conflicting duplicate.*rg'):
                    self.compile(Path(folder), rows, policy)

    def test_r_input_validation_and_matrix_parity(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        completed = subprocess.run([executable, str(ROOT / 'tests/test_self_duplicates.R'), str(ROOT)],
                                   capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)


if __name__ == '__main__':
    unittest.main()
