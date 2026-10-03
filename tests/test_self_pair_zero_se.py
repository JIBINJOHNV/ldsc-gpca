"""A degenerate self-rg must not discard an otherwise valid trait."""
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
from ldsc_gpca.result_qc import result_status
from ldsc_gpca.trait_selection import validate_duplicate_estimates
from test_results_csv import record, manifest
from test_results_precision import HEADER, block, row


ROOT = Path(__file__).resolve().parents[1]


def zero_self():
    return dict(record(), se=0., z=math.inf, p=0.)


class SelfPairZeroSETests(unittest.TestCase):
    def test_native_compilation_preserves_values_and_traits_in_every_mode(self):
        frame = pd.DataFrame([record(a, b) for a in ('A', 'B', 'C') for b in ('A', 'B', 'C')])
        frame.loc[frame.p1.eq(frame.p2), ['se', 'z', 'p']] = [0., math.inf, 0.]
        for action in ('error', 'report', 'drop_traits'):
            with self.subTest(action=action), tempfile.TemporaryDirectory() as folder:
                path = Path(folder, 'batch.results.csv')
                write_results_csv(frame, path)
                results.compile_results(folder, [path], manifest(('C', 'A', 'B')),
                                        result_failure_action=action)
                actual = pd.read_csv(Path(folder, 'ldsc_results.csv'), float_precision='round_trip')
                self.assertEqual(len(actual), 9)
                self.assertEqual(actual.p1.drop_duplicates().tolist(), ['C', 'A', 'B'])
                selves = actual.loc[actual.p1.eq(actual.p2)]
                self.assertTrue(selves.se.eq(0).all() and selves.z.eq(math.inf).all()
                                and selves.p.eq(0).all())
                self.assertTrue(selves.rg.eq(record()['rg']).all())
                status = pd.read_csv(Path(folder, 'LDSC_Pair_Status.csv'), keep_default_na=False)
                self.assertTrue(status.Status.eq('valid').all())
                self.assertTrue(status.loc[status.p1.eq(status.p2), 'Warning'].str.contains('SE=0').all())

    def test_exception_tolerance_and_invalid_neighboring_cases(self):
        for rg in (.99, 1., 1.01):
            status = result_status(pd.DataFrame([dict(zero_self(), rg=rg)]), 'fixture')
            self.assertEqual(status.Status.iloc[0], 'valid')
        invalid = [('p2', '/inputs/B.sumstats.gz'), ('rg', .98), ('rg', 1.02),
                   ('rg', math.inf), ('rg', math.nan), ('se', -.1), ('se', math.nan),
                   ('se', math.inf), ('z', 100.), ('z', -math.inf), ('z', math.nan),
                   ('p', .05), ('p', math.nan), ('p', math.inf), ('p', -1.),
                   ('h2_obs', 0.), ('h2_obs_se', 0.), ('h2_int_se', 0.),
                   ('gcov_int_se', 0.), ('h2_int', math.nan)]
        for field, value in invalid:
            with self.subTest(field=field, value=value):
                status = result_status(pd.DataFrame([dict(zero_self(), **{field: value})]), 'fixture')
                self.assertEqual(status.Status.iloc[0], 'failed_estimate')

    def test_duplicate_infinities_agree_but_other_nonfinite_values_conflict(self):
        original = zero_self()
        order = {original['p1']: 0}
        validate_duplicate_estimates(pd.DataFrame([original, original]), order)
        for z in (100., -math.inf, math.nan):
            with self.subTest(z=z), self.assertRaisesRegex(RuntimeError, 'Conflicting duplicate'):
                validate_duplicate_estimates(pd.DataFrame([original, dict(original, z=z)]), order)

    def test_legacy_detailed_and_table_only_native_zero_are_preserved(self):
        contents = [HEADER + row('A', z='inf'),
                    'Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                    + block('A', se='0.', z='inf')
                    + 'Summary of Genetic Correlation Results\n' + HEADER + row('A')]
        for content in contents:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as folder:
                path = Path(folder, 'batch.log')
                path.write_text(content)
                results.compile_results(folder, [path], manifest())
                actual = pd.read_csv(Path(folder, 'ldsc_results.csv'))
                self.assertEqual((actual.se.iloc[0], actual.z.iloc[0], actual.p.iloc[0]),
                                 (0., math.inf, 0.))

    def test_r_validation_and_audit(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable:
            self.skipTest('Rscript unavailable')
        outcome = subprocess.run([executable, str(ROOT/'tests/test_self_pair_zero_se.R'), str(ROOT)],
                                 capture_output=True, text=True)
        self.assertEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)


if __name__ == '__main__':
    unittest.main()
