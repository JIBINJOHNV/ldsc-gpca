"""Munging precision, atomic publication, and optional pinned-native integration.

Set LDSC_GPCA_TEST_MUNGE_SCRIPT to native munge_sumstats.py to run real munging.
"""
import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import runpy
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from ldsc_gpca import ldsc_munge, munging


NATIVE_SCRIPT = os.environ.get('LDSC_GPCA_TEST_MUNGE_SCRIPT')


class PrecisionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='munging precision ')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.path = self.root / 'out.sumstats.gz'

    def test_exact_values_column_order_missing_and_optional_frequency(self):
        frame = pd.DataFrame({'SNP': ['rs1', 'rs2', 'rs3', 'rs4'],
            'A1': ['A', 'C', 'A', None], 'A2': ['G', 'T', 'G', None],
            'Z': [8.9444, -1.234567891234567, 1e-100, np.nan],
            'N': [50000.123456789, 50000.0, 50000.0, np.nan],
            'FRQ': [0.12345678912345678, .2, .3, np.nan], 'IGNORED': 1})
        for keep_maf in (False, True):
            with self.subTest(keep_maf=keep_maf):
                ldsc_munge.write_munged_sumstats(frame, self.path, keep_maf=keep_maf)
                observed = pd.read_csv(self.path, sep='\t', float_precision='round_trip')
                expected = frame.drop(columns=['IGNORED'] + ([] if keep_maf else ['FRQ']))
                pd.testing.assert_frame_equal(observed, expected, check_exact=True)
                with gzip.open(self.path, 'rt') as handle:
                    self.assertEqual(handle.readlines()[-1].strip('\n'), 'rs4' + '\t' * (len(expected.columns)-1))
        self.assertEqual(list(self.root.glob('.ldsc-munge-*')), [])

    def test_threshold_is_applied_to_unrounded_values(self):
        frame = pd.DataFrame({'SNP': ['above', 'below', 'negative', 'zero', 'missing'],
                              'N': 50000., 'Z': [8.9444, 8.9442, -8.9444, 0., np.nan]})
        ldsc_munge.write_munged_sumstats(frame, self.path)
        destination = self.root / 'filtered.gz'
        result = munging._filter_one_munged_sumstats('A', self.path, destination,
                                                    self.root / 'excluded.gz', 80)
        self.assertEqual(pd.read_csv(destination, sep='\t').SNP.tolist(), ['below', 'zero', 'missing'])
        self.assertEqual(result['variants_missing_z'], 1)
        self.assertLess(float('%.3f' % 8.9444)**2, 80)
        self.assertGreater(8.9444**2, 80)

    def test_failed_write_preserves_old_output_and_removes_temporary(self):
        self.path.write_bytes(b'previous output')
        def fail(frame, path, **kwargs):
            Path(path).write_bytes(b'partial gzip')
            raise OSError('disk full')
        with patch.object(pd.DataFrame, 'to_csv', fail), self.assertRaisesRegex(OSError, 'disk full'):
            ldsc_munge.write_munged_sumstats(pd.DataFrame({'SNP':['rs1'], 'Z':[1.], 'N':[100.]}), self.path)
        self.assertEqual(self.path.read_bytes(), b'previous output')
        self.assertEqual(list(self.root.glob('.ldsc-munge-*')), [])

    def test_invalid_native_schema_is_rejected_before_writing(self):
        with self.assertRaisesRegex(RuntimeError, 'expected SNP, N and Z'):
            ldsc_munge.write_munged_sumstats(pd.DataFrame({'SNP': ['rs1']}), self.path)
        self.assertFalse(self.path.exists())

    def test_missing_native_script_has_clear_error(self):
        with patch.object(ldsc_munge.shutil, 'which', return_value=None), \
             self.assertRaisesRegex(RuntimeError, 'not installed'):
            ldsc_munge.main([])

    def test_managed_route_and_format_metadata(self):
        frame = pd.DataFrame({'gwas_name': ['A'], 'sample_size_column': ['NEF'],
                              'sample_prevalence': [np.nan]})
        def run(command, *args, **kwargs):
            tokens = shlex.split(command)
            self.assertEqual(tokens[:6], ['conda', 'run', '--no-capture-output', '--prefix', '/child env', 'python'])
            self.assertTrue(tokens[6].endswith('/ldsc_munge.py'))
            data = pd.DataFrame({'SNP': [f'rs{i}' for i in range(20)], 'N':50000., 'Z':np.arange(20)/7})
            ldsc_munge.write_munged_sumstats(data, tokens[-1] + '.sumstats.gz')
        with patch.object(munging, 'run_command', side_effect=run), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(munging.parallel_munge_sumstats(1, frame, str(self.root), 'hm3',
                                                            runtime={'prefix': '/child env'}), ['A'])
        metadata = json.loads((self.root/'ldsc_input/A.prevalence.json').read_text())
        self.assertEqual(metadata['float_format'], '%.17g')


@unittest.skipUnless(NATIVE_SCRIPT, 'set LDSC_GPCA_TEST_MUNGE_SCRIPT for native integration')
class NativeMungingTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='native munging ')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        previous_path = sys.path[:]
        try:
            sys.path.insert(0, str(Path(NATIVE_SCRIPT).resolve().parent))
            self.native = runpy.run_path(NATIVE_SCRIPT)
        finally:
            sys.path[:] = previous_path
        from scipy.stats import chi2
        z = np.array([8.9444, -8.9444, 1.234567891234567, -.8, .1, -.1, 0., np.nan, .2, -.2])
        data = pd.DataFrame({'SNP': [f'rs{i}' for i in range(len(z))], 'A1':'A', 'A2':'G',
                             'N':50000.123456789, 'Z':z, 'P':chi2.sf(z*z, 1), 'FRQ':.2})
        data.loc[8, 'A2'] = 'T'  # Strand-ambiguous, removed by native QC.
        data.loc[9, 'FRQ'] = .001  # Removed by native MAF filtering.
        data.to_csv(self.root/'input.tsv', sep='\t', index=False)
        pd.DataFrame({'SNP': data.SNP.tolist()+['reference_only'], 'A1':'A', 'A2':'G'}).to_csv(
            self.root/'hm3.tsv', sep='\t', index=False)
        self.arguments = ['--sumstats', str(self.root/'input.tsv'), '--merge-alleles', str(self.root/'hm3.tsv'),
                          '--signed-sumstats', 'Z,0', '--frq', 'FRQ', '--keep-maf']

    def test_native_qc_and_computed_values_preserved_before_rounding(self):
        before = sys.path[:]
        original = self.native['parser'].parse_args(self.arguments + ['--out', str(self.root/'original')])
        with contextlib.redirect_stdout(io.StringIO()):
            expected = self.native['munge_sumstats'](original, p=True)
            with patch.object(ldsc_munge.shutil, 'which', return_value=NATIVE_SCRIPT), \
                 patch.object(pd.DataFrame, 'to_csv', autospec=True, side_effect=pd.DataFrame.to_csv) as writes:
                ldsc_munge.main(self.arguments + ['--out', str(self.root/'precise')])
                self.assertEqual(writes.call_count, 1)
        self.assertEqual(sys.path, before)
        observed = pd.read_csv(self.root/'precise.sumstats.gz', sep='\t', float_precision='round_trip')
        columns = [c for c in expected if c in ('SNP','N','Z','A1','A2')] + ['FRQ']
        pd.testing.assert_frame_equal(observed, expected[columns], check_exact=True)
        from ldscore.parse import sumstats
        parsed = sumstats(str(self.root/'precise.sumstats.gz'), alleles=True)
        pd.testing.assert_frame_equal(parsed, observed[parsed.columns].dropna(),
                                      check_exact=False, rtol=1e-15, atol=0)
        rounded = pd.read_csv(self.root/'original.sumstats.gz', sep='\t')
        self.assertTrue((observed.Z.dropna() != rounded.Z.dropna()).any())
        self.assertEqual(observed.loc[observed.Z.isna(), 'SNP'].tolist(), ['rs7','rs8','rs9','reference_only'])
        self.assertLess(rounded.Z.iloc[0]**2, 80)
        self.assertGreater(observed.Z.iloc[0]**2, 80)
        self.assertIn('%.17g', (self.root/'precise.log').read_text())

    def test_help_and_invalid_native_inputs_preserve_exit_behavior(self):
        before = sys.path[:]
        with patch.object(ldsc_munge.shutil, 'which', return_value=NATIVE_SCRIPT), \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as result:
                ldsc_munge.main(['--help'])
            self.assertEqual(result.exception.code, 0)
            with self.assertRaises(SystemExit) as result:
                ldsc_munge.main(['--not-a-native-option'])
            self.assertEqual(result.exception.code, 2)
            with self.assertRaises(ValueError):
                ldsc_munge.main(['--out', str(self.root/'invalid')])
        self.assertEqual(sys.path, before)
        self.assertFalse((self.root/'invalid.sumstats.gz').exists())


if __name__ == '__main__':
    unittest.main()
