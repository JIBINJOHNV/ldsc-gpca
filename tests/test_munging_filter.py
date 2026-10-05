"""Missing munging placeholders must not prevent filtering usable variants."""
import contextlib
import gzip
import io
from pathlib import Path
import tempfile
import unittest

import pandas as pd
from ldsc_gpca.munging import _filter_one_munged_sumstats, filter_munged_sumstats


HEADER = 'SNP\tA1\tA2\tN\tZ\n'


class MungingFilterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'T.sumstats.gz'
        self.filtered = self.root/'filtered.gz'
        self.excluded = self.root/'excluded.gz'

    def write(self, text, path=None):
        (path or self.source).write_bytes(gzip.compress(text.encode()))

    def read(self, path):
        return gzip.decompress(path.read_bytes()).decode()

    def run_filter(self):
        return _filter_one_munged_sumstats('T', str(self.source), str(self.filtered), str(self.excluded), 64)

    def test_missing_tokens_and_reference_placeholders_preserve_text_and_counts(self):
        valid = 'rs1\tA\tG\t50000\t1.23456789\nrs_boundary\tA\tC\t50000\t-8\n'
        missing = 'rs_placeholder\t\t\t\t\n' + ''.join(
            f'rs_missing{i}\tC\tT\t50000\t{token}\n' for i, token in enumerate(['', 'NA', 'nan', 'NaN', '.']))
        self.write(HEADER + valid + missing + 'rs_large\tA\tG\t50000\t9\n')
        original = self.source.read_bytes()
        summary = self.run_filter()
        self.assertEqual(self.read(self.filtered), HEADER + valid + missing)
        self.assertEqual(self.read(self.excluded), 'gwas_name\tSNP\tZ\tCHISQ\nT\trs_large\t9\t81\n')
        self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual([summary[k] for k in ('variants_before','variants_removed','variants_after','variants_missing_z')],
                         [9, 1, 8, 6])
        self.assertEqual(summary['maximum_chisq'], 81)
        # Python LDSC uses this whitespace/NA parsing and drops missing rows.
        parsed = pd.read_csv(self.filtered, sep=r'\s+', na_values='.',
                             dtype={'SNP':str,'N':float,'Z':float,'A1':str,'A2':str}).dropna()
        self.assertEqual(parsed.SNP.tolist(), ['rs1','rs_boundary'])
        self.assertEqual(parsed.Z.tolist(), [1.23456789,-8.])

    def test_normal_zero_boundary_and_large_z_keep_existing_behavior(self):
        self.write(HEADER + ''.join(f'rs{i}\tA\tG\t50000\t{z}\n' for i,z in enumerate([0,8,-8,8.01])))
        summary = self.run_filter()
        self.assertEqual(summary['variants_missing_z'], 0)
        self.assertEqual(summary['variants_after'], 3)
        self.assertEqual(summary['variants_removed'], 1)
        self.assertEqual(pd.read_csv(self.filtered,sep='\t').Z.tolist(), [0,8,-8])

    def test_no_finite_survivors_fail_without_overwriting_outputs(self):
        for body, message in [('', 'contains no variants'),
                              ('rs_missing\t\t\t\t\n', 'no variants with finite Z'),
                              ('rs_large\tA\tG\t50000\t9\n', 'no variants with finite Z'),
                              ('rs_missing\t\t\t\t\nrs_large\tA\tG\t50000\t9\n', 'no variants with finite Z')]:
            with self.subTest(body=body):
                self.write(HEADER+body)
                self.filtered.write_bytes(b'previous filtered')
                self.excluded.write_bytes(b'previous excluded')
                with self.assertRaisesRegex(ValueError, message): self.run_filter()
                self.assertEqual(self.filtered.read_bytes(), b'previous filtered')
                self.assertEqual(self.excluded.read_bytes(), b'previous excluded')
                self.assertEqual(list(self.root.glob('*.tmp')), [])

    def test_invalid_values_and_malformed_rows_still_fail(self):
        for value in ['bad', 'Inf', '-inf', '1e309', '1e200']:
            with self.subTest(value=value):
                self.write(HEADER+'rs_ok\tA\tG\t50000\t1\n'+f'rs_bad\tA\tG\t50000\t{value}\n')
                with self.assertRaisesRegex(ValueError, 'row 3'): self.run_filter()
                self.assertFalse(self.filtered.exists())
        self.write(HEADER+'rs_missing\t\t\t\n')
        with self.assertRaisesRegex(ValueError, 'malformed'): self.run_filter()
        self.write('SNP\tA1\tA2\tN\nrs1\tA\tG\t50000\n')
        with self.assertRaisesRegex(ValueError, 'one Z column'): self.run_filter()

    def test_line_endings_final_newline_and_whitespace_inputs(self):
        for text, expected in [
            (HEADER.replace('\n','\r\n')+'rs1\tA\tG\t50000\t1\r\nrs2\t\t\t\t',
             HEADER.replace('\n','\r\n')+'rs1\tA\tG\t50000\t1\r\nrs2\t\t\t\t\n'),
            ('SNP A1 A2 N Z\nrs1 A G 50000 1\nrs2 A G 50000 NA\n',
             'SNP A1 A2 N Z\nrs1 A G 50000 1\nrs2 A G 50000 NA\n')]:
            with self.subTest(text=text):
                self.write(text); summary=self.run_filter()
                self.assertEqual(self.read(self.filtered), expected)
                self.assertEqual(summary['variants_missing_z'], 1)

    def test_parallel_summary_preserves_manifest_order_and_audits_missing_separately(self):
        for name in ['B','A']:
            self.write(HEADER+'rs1\tA\tG\t50000\t1\nrs2\t\t\t\t\nrs3\tA\tG\t50000\t9\n',
                       self.root/(name+'.sumstats.gz'))
        output=self.root/'out'
        with contextlib.redirect_stdout(io.StringIO()):
            folder=filter_munged_sumstats(2,pd.DataFrame({'gwas_name':['B','A']}),str(self.root),str(output),64)
        summary=pd.read_csv(output/'LDSC_ChiSquare_Filter_Summary.csv')
        self.assertEqual(summary.gwas_name.tolist(), ['B','A'])
        self.assertEqual(summary.variants_missing_z.tolist(), [1,1])
        self.assertEqual(summary.variants_before.tolist(), [3,3])
        self.assertEqual(summary.variants_after.tolist(), [2,2])
        excluded=pd.read_csv(output/'LDSC_ChiSquare_Excluded_Variants.tsv.gz',sep='\t')
        self.assertEqual(excluded.gwas_name.tolist(), ['B','A'])
        self.assertEqual(excluded.SNP.tolist(), ['rs3','rs3'])
        self.assertIn('rs2\t\t\t\t\n',self.read(Path(folder)/'B.sumstats.gz'))


if __name__ == '__main__': unittest.main()
