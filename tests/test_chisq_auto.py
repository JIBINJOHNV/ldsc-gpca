"""Automatic cutoffs use each trait's complete LD-matched SNPs before filtering."""
import contextlib
import gzip
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import ldsc_cli, munging


HEADER = 'SNP\tA1\tA2\tN\tZ\n'


class ChiSquareArgumentTests(unittest.TestCase):
    def test_positive_integers_auto_and_disabled_default(self):
        base = ['--input', 'traits.csv', '--outdir', 'out', '--ld_ref', 'ref']
        self.assertIsNone(ldsc_cli.parser.parse_args(base).chisq_max)
        for value, expected in [('1', 1), ('80', 80), ('500', 500), ('auto', 'auto')]:
            with self.subTest(value=value):
                self.assertEqual(ldsc_cli.parser.parse_args(base + ['--chisq_max', value]).chisq_max, expected)
        for value in ['0', '-1', '80.5', '80.0', 'nan', 'inf', '1e3', 'invalid']:
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()) as error:
                with self.assertRaises(SystemExit):
                    ldsc_cli.parser.parse_args(base + ['--chisq_max', value])
                self.assertIn('positive integer (>0) or auto', error.getvalue())
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            ldsc_cli.parser.parse_args(base + ['--chisq_max'])


class AutomaticChiSquareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ld = self.root / 'ld'; self.ld.mkdir()
        self.weights = self.root / 'weights'; self.weights.mkdir()
        self.snps = ['rs_low', 'rs_mid', 'rs_high', 'rs_missing', 'rs_incomplete', 'rs_placeholder']
        for directory, extra in [(self.ld, ['rs_ref_only']), (self.weights, ['rs_weight_only'])]:
            for chrom in range(1, 23):
                rows = self.snps + extra if chrom == 1 else []
                self.write(directory / f'{chrom}.l2.ldscore.gz', 'CHR SNP BP L2\n' +
                           ''.join(f'{chrom} {snp} {i+1} 1.5\n' for i, snp in enumerate(rows)))

    def write(self, path, text):
        # Include a gzip filename so tiny fixtures pass the existing 50-byte preflight.
        with path.open('wb') as stream, gzip.GzipFile(
                filename='automatic-chisq-test-fixture', mode='wb', fileobj=stream) as archive:
            archive.write(text.encode())

    def run_filter(self, names, cutoff='auto', **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            folder = munging.filter_munged_sumstats(2, pd.DataFrame({'gwas_name': names}),
                str(self.root), str(self.root / 'out'), cutoff,
                ld_ref_dir=str(self.ld), ld_weights_dir=str(self.weights), **kwargs)
        return Path(folder), pd.read_csv(self.root / 'out' / 'LDSC_ChiSquare_Filter_Summary.csv')

    def test_four_traits_use_independent_complete_matched_maxima_and_keep_boundary(self):
        names = ['D', 'B', 'A', 'C']
        sizes = {'A': 50000, 'B': 80000, 'C': 144000, 'D': 500250}
        original = {}
        for name, n in sizes.items():
            text = HEADER + f'rs_low\tA\tG\t{n/2}\t0\nrs_mid\tA\tG\t{n}\t-12\nrs_high\tA\tG\t{n}\t23\n'
            text += ''.join(f'{snp}\tA\tG\t9999999\t1\n' for snp in
                            ['rs_unmatched', 'rs_ref_only', 'rs_weight_only'])
            text += 'rs_missing\tA\tG\t9999999\tNA\nrs_incomplete\tNA\tG\t9999999\t1\nrs_placeholder\t\t\t\t\n'
            path = self.root / f'{name}.sumstats.gz'
            self.write(path, text); original[name] = path.read_bytes()
        with patch.object(munging, '_matched_ld_snps', wraps=munging._matched_ld_snps) as load:
            folder, summary = self.run_filter(names)
        self.assertEqual(load.call_count, 1)
        self.assertEqual(summary.gwas_name.tolist(), names)
        self.assertEqual(summary.chisq_max.tolist(), [500.25, 80, 80, 144])
        self.assertEqual(summary.maximum_matched_n.tolist(), [sizes[name] for name in names])
        self.assertEqual(summary.complete_ld_matched_rows.tolist(), [3] * 4)
        self.assertEqual(summary.variants_missing_z.tolist(), [2] * 4)
        self.assertEqual(summary.variants_removed.tolist(), [1, 2, 2, 1])
        self.assertTrue(summary.threshold_mode.eq('auto').all())
        for name in names:
            with gzip.open(folder / f'{name}.sumstats.gz', 'rt') as stream:
                kept = stream.read()
            self.assertEqual('rs_mid\t' in kept, name in ('C', 'D'))
            self.assertNotIn('rs_high\t', kept)
            self.assertIn('rs_placeholder\t\t\t\t\n', kept)
            self.assertEqual((self.root / f'{name}.sumstats.gz').read_bytes(), original[name])

    def test_final_high_n_row_sets_cutoff_even_when_that_row_is_removed(self):
        source = self.root / 'T.sumstats.gz'
        self.write(source, HEADER + 'rs_low\tA\tG\t1000\t12\nrs_high\tA\tG\t200500\t30\n')
        folder, summary = self.run_filter(['T'])
        self.assertEqual(summary.chisq_max.tolist(), [200.5])
        self.assertEqual(pd.read_csv(folder / source.name, sep='\t').SNP.tolist(), ['rs_low'])

    def test_all_missing_markers_excluded_from_maximum_and_whitespace_supported(self):
        for sep in ['\t', ' ']:
            for token in ['NA', 'nan', 'NaN', '.'] + ([''] if sep == '\t' else []):
                with self.subTest(sep=sep, token=token):
                    text = HEADER + f'rs_low\tA\tG\t50000\t1\nrs_missing\tA\tG\t999999\t{token}\n'
                    self.write(self.root / 'T.sumstats.gz', text.replace('\t', sep))
                    _, summary = self.run_filter(['T'])
                    self.assertEqual(summary.chisq_max.tolist(), [80])
                    self.assertEqual(summary.maximum_matched_n.tolist(), [50000])

    def test_invalid_n_or_missing_required_columns_fail_before_outputs(self):
        for n in ['bad', '0', '-1', 'inf', '-inf', '1e309']:
            with self.subTest(n=n):
                self.write(self.root / 'T.sumstats.gz', HEADER + f'rs_low\tA\tG\t{n}\t1\n')
                with self.assertRaisesRegex(ValueError, 'N'):
                    self.run_filter(['T'])
                self.assertFalse((self.root / 'out/ldsc_input_chisq_filtered/T.sumstats.gz').exists())
        for text in ['SNP\tA1\tA2\tZ\nrs_low\tA\tG\t1\n', 'SNP\tN\tN\tZ\n']:
            self.write(self.root / 'T.sumstats.gz', text)
            with self.assertRaisesRegex(ValueError, 'requires unique'):
                self.run_filter(['T'])

    def test_padded_missing_fields_do_not_inflate_cutoff(self):
        self.write(self.root / 'T.sumstats.gz', HEADER +
                   'rs_low\tA\tG\t50000\t1\nrs_missing\tA\tG\t999999\t NA \n' +
                   'rs_incomplete\t   \tG\t999999\t1\n')
        _, summary = self.run_filter(['T'])
        self.assertEqual(summary.chisq_max.tolist(), [80])
        self.assertEqual(summary.complete_ld_matched_rows.tolist(), [1])

    def test_empty_no_match_all_missing_and_all_excluded_fail(self):
        for body, message in [('', 'no complete rows'),
                              ('rs_other\tA\tG\t50000\t1\n', 'no complete rows'),
                              ('rs_missing\t\t\t\t\n', 'no complete rows'),
                              ('rs_low\tA\tG\t50000\t10\n', 'no variants with finite Z')]:
            with self.subTest(body=body):
                self.write(self.root / 'T.sumstats.gz', HEADER + body)
                with self.assertRaisesRegex(ValueError, message):
                    self.run_filter(['T'])

    def test_references_required_only_for_auto_and_weights_default_to_reference(self):
        expected = set(self.snps) | {'rs_ref_only'}
        self.assertEqual(set(munging._matched_ld_snps(str(self.ld), None)), expected)
        with self.assertRaisesRegex(ValueError, 'requires an LD-score reference'):
            munging._matched_ld_snps(None, None)
        for chrom in range(1, 23):
            self.write(self.weights / f'{chrom}.l2.ldscore.gz', 'CHR SNP BP L2\n1 unrelated 1 2\n')
        with self.assertRaisesRegex(ValueError, 'no shared SNPs'):
            munging._matched_ld_snps(str(self.ld), str(self.weights))
        (self.ld / '22.l2.ldscore.gz').unlink()
        with self.assertRaises(FileNotFoundError):
            munging._matched_ld_snps(str(self.ld), None)
        self.write(self.root / 'T.sumstats.gz', HEADER + 'rs_low\tA\tG\t50000\t1\n')
        with patch.object(munging, '_matched_ld_snps', side_effect=AssertionError('unused')):
            _, summary = self.run_filter(['T'], 80)
        self.assertEqual(summary.threshold_mode.tolist(), ['fixed'])
        self.assertTrue(summary.maximum_matched_n.isna().all())


if __name__ == '__main__':
    unittest.main()
