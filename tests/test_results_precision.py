"""Regression tests for precision loss in LDSC's displayed summary table."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from ldsc_gpca import results


HEADER = 'p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se\n'


def block(target, rg='1.', se='1.9647e-06', z='508990.3925', p='0.'):
    return (
        'Computing rg for phenotype 2/3\n'
        f'Reading summary statistics from /inputs/{target}.sumstats.gz ...\n'
        f'Genetic Correlation: {rg} ({se})\n'
        f'Z-score: {z}\nP: {p}\n'
    )


def row(target, rg='1.0000', se='0.0000', z='508990.3925', p='0.0000'):
    return f'/inputs/A.sumstats.gz /inputs/{target}.sumstats.gz {rg} {se} {z} {p} .2 .03 1.01 .006 1.01 .006\n'


class ResultsPrecisionTests(unittest.TestCase):
    def compile(self, content, traits=('A',)):
        with tempfile.TemporaryDirectory() as folder:
            log = Path(folder) / 'test.log'
            log.write_text(content)
            metadata = pd.DataFrame({
                'gwas_name': traits,
                'ref': ['yes'] + ['no'] * (len(traits) - 1),
                'pop_prevalence': [float('nan')] * len(traits),
            })
            results.compile_results(folder, [str(log)], metadata)
            return pd.read_csv(Path(folder) / 'ldsc_results.csv')

    def test_recovers_positive_self_se_in_scientific_notation(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A') + 'Summary of Genetic Correlation Results\n' + HEADER + row('A'))
        compiled = self.compile(content)
        self.assertEqual(compiled.se.iloc[0], 1.9647e-06)
        self.assertEqual(compiled.rg.iloc[0], 1.)
        self.assertEqual(compiled.z.iloc[0], 508990.3925)
        self.assertEqual(compiled.h2_int.iloc[0], 1.01)

    def test_matches_trait_pairs_not_summary_row_positions(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A') + block('B', '.25', '.123456', '2.025', '4.285e-02')
                   + 'Summary of Genetic Correlation Results\n' + HEADER
                   + row('B', '.2500', '.1235', '2.0250', '.0429') + row('A'))
        compiled = self.compile(content, ('A', 'B')).set_index('p2')
        self.assertEqual(compiled.loc['A', 'se'], 1.9647e-06)
        self.assertEqual(compiled.loc['B', 'se'], .123456)
        self.assertEqual(compiled.loc['B', 'p'], .04285)

    def test_preserves_tiny_p_value_and_se(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A', se='2.0071e-08', z='49823121', p='3.2e-100')
                   + 'Summary of Genetic Correlation Results\n' + HEADER + row('A'))
        compiled = self.compile(content)
        self.assertEqual(compiled.se.iloc[0], 2.0071e-08)
        self.assertEqual(compiled.p.iloc[0], 3.2e-100)

    def test_table_only_positive_se_remains_supported(self):
        compiled = self.compile(HEADER + row('A', se='.1000'))
        self.assertEqual(compiled.se.iloc[0], .1)

    def test_repeated_table_headers_are_ignored(self):
        compiled = self.compile(HEADER + HEADER + row('A', se='.1000'))
        self.assertEqual(len(compiled), 1)
        self.assertEqual(compiled.se.iloc[0], .1)

    def test_table_only_zero_fails_without_imputation(self):
        with self.assertRaisesRegex(RuntimeError, 'summary-table zero may be rounding'):
            self.compile(HEADER + row('A'))

    def test_actual_zero_negative_or_nan_se_still_fails(self):
        for se in ('0.', '-1e-06', 'nan', 'inf'):
            content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                       + block('A', se=se) + 'Summary of Genetic Correlation Results\n'
                       + HEADER + row('A', se='.1000'))
            with self.subTest(se=se), self.assertRaisesRegex(RuntimeError, 'No SE was imputed'):
                self.compile(content)

    def test_incomplete_block_fails(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A').replace('P: 0.\n', '')
                   + 'Summary of Genetic Correlation Results\n' + HEADER + row('A'))
        with self.assertRaisesRegex(RuntimeError, 'Incomplete detailed'):
            self.compile(content)

    def test_malformed_correlation_fails(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A').replace('1. (1.9647e-06)', 'bad')
                   + 'Summary of Genetic Correlation Results\n' + HEADER + row('A'))
        with self.assertRaisesRegex(RuntimeError, 'Malformed detailed'):
            self.compile(content)

    def test_nonnumeric_detail_fails(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A', se='bad') + 'Summary of Genetic Correlation Results\n'
                   + HEADER + row('A'))
        with self.assertRaisesRegex(RuntimeError, 'Non-numeric detailed'):
            self.compile(content)

    def test_wrong_pair_fails(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('B') + 'Summary of Genetic Correlation Results\n' + HEADER + row('A'))
        with self.assertRaisesRegex(RuntimeError, 'absent from summary table'):
            self.compile(content)

    def test_missing_detail_for_valid_summary_pair_fails(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A') + 'Summary of Genetic Correlation Results\n'
                   + HEADER + row('A') + row('B', se='.1000'))
        with self.assertRaisesRegex(RuntimeError, 'Missing detailed'):
            self.compile(content, ('A', 'B'))

    def test_duplicate_detail_blocks_fail(self):
        content = ('Reading summary statistics from /inputs/A.sumstats.gz ...\n'
                   + block('A') + block('A') + 'Summary of Genetic Correlation Results\n'
                   + HEADER + row('A'))
        with self.assertRaisesRegex(RuntimeError, 'Duplicate detailed'):
            self.compile(content)


if __name__ == '__main__':
    unittest.main()
