"""Independent numerical oracle and input-preservation checks for annotations."""
import contextlib
import csv
import gzip
import io
import math
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import polars as pl
from polars.testing import assert_frame_equal
from ldsc_gpca.normalization import add_trait_wide_columns, main, DERIVED_COLUMNS
from ldsc_gpca.ldsc_export import NATIVE_COLUMNS

ROOT = Path(__file__).resolve().parents[1]


def fixture():
    # Four traits, different self/pair h2; include negative and >1 correlations.
    traits = ['D', 'NA', 'B', 'C']
    h2 = [.16, .25, .36, .49]
    rows = []
    for i, a in enumerate(traits):
        for j, b in enumerate(traits):
            rg = 1.0 if i == j else (-.3 if {i, j} == {0, 1} else 1.02 if {i, j} == {2, 3} else .4)
            rows.append(dict(p1=a, p2=b, rg=rg, h2_obs=h2[j] * (1 if i == j else .8 + i/10),
                             se=.05, z=rg/.05, p=1e-310, metadata='0001'))
    return pl.DataFrame(rows)


def native_fixture():
    rows = fixture().to_dicts()
    lookup = {(r['p1'], r['p2']): r['h2_obs'] for r in rows}
    for r in rows:
        r['h2_p1_pair_obs'] = lookup[r['p2'], r['p1']]
        r['h2_p2_pair_obs'] = r['h2_obs']
        r['gcov_native_obs'] = r['rg'] * math.sqrt(r['h2_p1_pair_obs'] * r['h2_p2_pair_obs'])
    return pl.DataFrame(rows)


class NormalizationTests(unittest.TestCase):
    def test_compiler_adds_annotations_without_changing_native_columns(self):
        import pandas as pd
        from ldsc_gpca.results import compile_results
        from ldsc_gpca.ldsc_export import write_results_csv
        from test_results_csv import record, manifest
        data = fixture()
        rows = []
        for row in data.to_dicts():
            native = record(row['p1'], row['p2'])
            native.update({key: row[key] for key in ('rg', 'se', 'z', 'p', 'h2_obs')})
            rows.append(native)
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()):
            source = Path(folder)/'batch.results.csv'
            write_results_csv(pd.DataFrame(rows), source)
            out = compile_results(folder, [source], manifest(('D', 'NA', 'B', 'C')))
            actual = pl.read_csv(out, null_values='')
            expected = add_trait_wide_columns(data)
            assert_frame_equal(actual.select('p1','p2',*DERIVED_COLUMNS),
                               expected.select('p1','p2',*DERIVED_COLUMNS))

    def test_four_trait_oracle_preserves_every_original_value_and_order(self):
        source = fixture().reverse()
        observed = add_trait_wide_columns(source)
        assert_frame_equal(observed.select(source.columns), source)
        rows = {(r['p1'], r['p2']): r for r in source.to_dicts()}
        for r in observed.to_dicts():
            a, b = r['p1'], r['p2']
            gcov = r['rg'] * math.sqrt(r['h2_obs'] * rows[b, a]['h2_obs'])
            expected = 1 if a == b else gcov / math.sqrt(rows[a, a]['h2_obs'] * rows[b, b]['h2_obs'])
            self.assertAlmostEqual(r['gcov_pair'], gcov, places=14)
            self.assertAlmostEqual(r['rg_trait_wide'], expected, places=14)
            self.assertEqual(r['normalization_status'], 'calculated')
        self.assertLess(observed.filter((pl.col('p1') == 'D') & (pl.col('p2') == 'NA'))['rg_trait_wide'][0], 0)
        self.assertGreater(observed.filter((pl.col('p1') == 'B') & (pl.col('p2') == 'C'))['rg_trait_wide'][0], 1)
        assert_frame_equal(add_trait_wide_columns(observed), observed)

    def test_native_covariance_works_without_reverse_rows_and_preserves_inputs(self):
        source = native_fixture().filter(pl.col('p1') <= pl.col('p2')).reverse()
        result = add_trait_wide_columns(source)
        assert_frame_equal(result.select(source.columns), source)
        expected = add_trait_wide_columns(fixture())
        for row in result.to_dicts():
            old = expected.filter((pl.col('p1') == row['p1']) & (pl.col('p2') == row['p2'])).row(0,named=True)
            self.assertEqual(row['normalization_status'], 'calculated')
            for column in ('gcov_pair','rg_trait_wide'):
                self.assertAlmostEqual(row[column], old[column], places=14)
        assert_frame_equal(add_trait_wide_columns(result), result)

    def test_native_normalization_uses_unconverted_h2_with_mixed_reported_scales(self):
        source = native_fixture()
        converted = source.with_columns(
            pl.when(pl.col('p2')=='D').then(pl.col('h2_obs')*2).otherwise(None).alias('h2_liab'),
            pl.when(pl.col('p2')!='D').then(pl.col('h2_obs')).otherwise(None).alias('h2_obs'))
        native_result = add_trait_wide_columns(converted)
        legacy_result = add_trait_wide_columns(converted.drop(NATIVE_COLUMNS))
        self.assertEqual(set(native_result['normalization_status']), {'calculated'})
        for column in ('gcov_pair','rg_trait_wide'):
            assert_frame_equal(native_result.select(column),legacy_result.select(column))
        assert_frame_equal(native_result.select('rg_trait_wide'),add_trait_wide_columns(source).select('rg_trait_wide'))
        bad = converted.with_columns(pl.when(pl.col('p1')!='D').then(pl.col('h2_liab')*2)
                                     .otherwise(pl.col('h2_liab')).alias('h2_liab'))
        self.assertIn('inconsistent_native_estimate',add_trait_wide_columns(bad)['normalization_status'])

    def test_native_invalid_partial_conflicting_and_inconsistent_values_are_explicit(self):
        source = native_fixture()
        for column,value in [('gcov_native_obs',float('inf')),('gcov_native_obs',None),
                             ('h2_p1_pair_obs',0),('h2_p2_pair_obs',-.1),('h2_p1_pair_obs',float('nan'))]:
            with self.subTest(column=column,value=value):
                bad = source.with_columns(pl.lit(value,dtype=pl.Float64).alias(column))
                result = add_trait_wide_columns(bad)
                self.assertEqual(set(result['normalization_status']), {'invalid_native_estimate'})
                self.assertEqual(result['rg_trait_wide'].null_count(), len(result))
        bad=source.with_columns((pl.col('gcov_native_obs')+.01).alias('gcov_native_obs'))
        self.assertEqual(set(add_trait_wide_columns(bad)['normalization_status']), {'inconsistent_native_estimate'})
        same=pl.concat([source,source.head(1)])
        self.assertEqual(set(add_trait_wide_columns(same)['normalization_status']), {'calculated'})
        conflict=pl.concat([source,source.head(1).with_columns(pl.lit(.8).alias('h2_p1_pair_obs'))])
        result=add_trait_wide_columns(conflict).filter((pl.col('p1')=='D')|(pl.col('p2')=='D'))
        self.assertEqual(set(result['normalization_status']), {'ambiguous_native_estimate'})
        missing_self=source.filter(~((pl.col('p1')=='D')&(pl.col('p2')=='D')))
        self.assertIn('missing_self_pair', add_trait_wide_columns(missing_self)['normalization_status'])

    def test_blank_native_fields_support_legacy_exports_and_zero_covariance_is_valid(self):
        source=fixture().with_columns(*[pl.lit(None,dtype=pl.Float64).alias(c) for c in NATIVE_COLUMNS])
        assert_frame_equal(add_trait_wide_columns(source).select(*DERIVED_COLUMNS),
                           add_trait_wide_columns(fixture()).select(*DERIVED_COLUMNS))
        zero=native_fixture().with_columns(*[pl.when(pl.col('p1')!=pl.col('p2')).then(0.)
                            .otherwise(pl.col(c)).alias(c) for c in ('rg','gcov_native_obs')])
        result=add_trait_wide_columns(zero)
        self.assertEqual(set(result['normalization_status']), {'calculated'})
        self.assertTrue((result.filter(pl.col('p1')!=pl.col('p2'))['rg_trait_wide']==0).all())

    def test_missing_reverse_or_self_is_explicit_and_never_imputed(self):
        source = fixture()
        for a, b, pair, status in [('D','NA',('NA','D'),'missing_reverse_pair'),
                                  ('D','D',('NA','D'),'missing_self_pair')]:
            data = source.filter(~((pl.col('p1') == a) & (pl.col('p2') == b)))
            row = add_trait_wide_columns(data).filter((pl.col('p1') == pair[0]) & (pl.col('p2') == pair[1])).row(0, named=True)
            self.assertEqual(row['normalization_status'], status)
            self.assertIsNone(row['rg_trait_wide'])

    def test_invalid_or_ambiguous_estimates_are_unavailable(self):
        for value in [0, -.1, float('inf'), float('nan'), None]:
            with self.subTest(value=value):
                data = fixture().with_columns(pl.lit(value, dtype=pl.Float64).alias('h2_obs'))
                result = add_trait_wide_columns(data)
                self.assertTrue(all(s == 'invalid_estimate' for s in result['normalization_status']))
                self.assertEqual(result['rg_trait_wide'].null_count(), len(data))
        data = fixture().with_columns(pl.col('h2_obs').alias('h2_liab'))
        self.assertEqual(set(add_trait_wide_columns(data)['normalization_status']), {'incompatible_or_ambiguous_scale'})
        for value in [float('inf'), float('nan'), None]:
            result = add_trait_wide_columns(fixture().with_columns(pl.lit(value, dtype=pl.Float64).alias('rg')))
            self.assertEqual(result['rg_trait_wide'].null_count(), len(result))

    def test_scale_conversion_cancels_and_mixed_trait_scales_are_supported(self):
        source = fixture()
        converted = source.with_columns((pl.col('h2_obs') * 2).alias('h2_liab')).drop('h2_obs')
        expected = add_trait_wide_columns(source)['rg_trait_wide']
        self.assertEqual(add_trait_wide_columns(converted)['rg_trait_wide'].round(14).to_list(), expected.round(14).to_list())
        mixed = source.with_columns(pl.when(pl.col('p2') == 'D').then(pl.col('h2_obs')*2).otherwise(None).alias('h2_liab'),
                                    pl.when(pl.col('p2') != 'D').then(pl.col('h2_obs')).otherwise(None).alias('h2_obs'))
        self.assertEqual(add_trait_wide_columns(mixed)['rg_trait_wide'].round(14).to_list(), expected.round(14).to_list())
        mismatch = mixed.with_columns(pl.when((pl.col('p1') == 'D') & (pl.col('p2') == 'D')).then(None).otherwise(pl.col('h2_liab')).alias('h2_liab'),
                                      pl.when((pl.col('p1') == 'D') & (pl.col('p2') == 'D')).then(.16).otherwise(pl.col('h2_obs')).alias('h2_obs'))
        self.assertIn('incompatible_or_ambiguous_scale', add_trait_wide_columns(mismatch)['normalization_status'])

    def test_duplicate_h2_ambiguity_does_not_choose_an_arbitrary_denominator(self):
        source = fixture()
        same = pl.concat([source, source.head(1)])
        self.assertEqual(len(add_trait_wide_columns(same)), len(same))
        self.assertEqual(set(add_trait_wide_columns(same)['normalization_status']), {'calculated'})
        bad = pl.concat([source, source.head(1).with_columns(pl.lit(.9).alias('h2_obs'))])
        observed = add_trait_wide_columns(bad).filter((pl.col('p1') == 'D') | (pl.col('p2') == 'D'))
        self.assertEqual(set(observed['normalization_status']), {'ambiguous_heritability'})

    def test_cli_in_place_preserves_numeric_text_and_rejects_bad_headers(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'file[1].csv'
            fixture().write_csv(path)
            original = list(csv.DictReader(io.StringIO(path.read_text())))
            with contextlib.redirect_stdout(io.StringIO()): main(['--input',str(path),'--out',str(path)])
            after = list(csv.DictReader(io.StringIO(path.read_text())))
            for old, new in zip(original, after):
                self.assertEqual(old, {key:new[key] for key in old})
            compressed = Path(folder)/'file[2].csv.gz'
            with gzip.open(compressed, 'wt') as handle: handle.write(path.read_text())
            with contextlib.redirect_stdout(io.StringIO()):
                main(['--input',str(compressed),'--out',str(path)])
            self.assertEqual(list(csv.DictReader(io.StringIO(path.read_text()))), after)
            for header in ['', 'p1,p1,rg\n', 'p1,,rg\n']:
                path.write_text(header)
                with self.assertRaisesRegex(ValueError, 'headers'): main(['--input',str(path),'--out',str(path)])
                self.assertEqual(path.read_text(), header)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(['--input',str(path),'--out',str(path)+'.gz'])

    def test_r_selection_and_default_compatibility(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable: self.skipTest('Rscript unavailable')
        run = subprocess.run([executable, str(ROOT/'tests/test_normalization.R'), str(ROOT)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout+run.stderr)


if __name__ == '__main__': unittest.main()
