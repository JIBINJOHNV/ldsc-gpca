"""Automatic LDSC trait exclusion, execution recovery and strict R handoff."""
import contextlib
import gzip
import io
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import ldsc_cli, pairwise, restart, results
from ldsc_gpca.ldsc_export import write_results_csv
from test_results_csv import record, manifest

ROOT = Path(__file__).resolve().parents[1]
TRAITS = ('A', 'B', 'C', 'D')


def fixture(kind='valid'):
    rows = pd.DataFrame([record(a, b) for a in TRAITS for b in TRAITS])
    names = {c: rows[c].map(results._trait_name) for c in ('p1', 'p2')}
    if kind in ('negative_self', 'all_failed', 'low_signal'):
        mask = names['p1'].eq(names['p2'])
        if kind != 'all_failed': mask &= names['p1'].eq('B')
        rows.loc[mask, 'h2_obs'] = .001 if kind == 'low_signal' else -.03
    if kind in ('negative_self', 'all_failed'):
        mask = pd.Series(True, index=rows.index) if kind == 'all_failed' else names['p1'].eq('B') | names['p2'].eq('B')
        rows.loc[mask, ['rg', 'se', 'z', 'p']] = float('nan')
    edges = {'tie': [('B', 'C')], 'lower_z': [('B', 'C')],
             'degree': [('A', 'B'), ('A', 'C')], 'missing_pair': [('B', 'C')]}.get(kind, [])
    for a, b in edges:
        mask = (names['p1'].eq(a) & names['p2'].eq(b)) | (names['p1'].eq(b) & names['p2'].eq(a))
        if kind == 'missing_pair': rows = rows.loc[~mask]
        else: rows.loc[mask, ['rg', 'se', 'z', 'p']] = float('nan')
    if kind == 'lower_z': rows.loc[names['p2'].eq('B'), 'h2_obs'] = .001
    if kind == 'missing_reverse': rows = rows.loc[~(names['p1'].eq('B') & names['p2'].eq('C'))]
    if kind == 'missing_self': rows = rows.loc[~(names['p1'].eq('B') & names['p2'].eq('B'))]
    if kind == 'out_of_range': rows.loc[names['p1'].ne(names['p2']), 'rg'] = 1.1
    return rows.copy()


class DropTraitsTests(unittest.TestCase):
    def compile(self, root, frame, names=TRAITS, metadata=None, source_manifest=None):
        source = root / 'batch.csv'; write_results_csv(frame, source)
        with contextlib.redirect_stdout(io.StringIO()):
            return results.compile_results(root, [source], metadata if metadata is not None else manifest(names),
                result_failure_action='drop_traits', input_manifest=source_manifest)

    def test_failed_self_removes_only_failed_trait_and_preserves_every_retained_value(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); original = fixture('negative_self')
            source_manifest = pd.DataFrame({'note': list('dcba'), 'traitname': ['D', 'C', 'B', 'A'],
                'ref': 'yes', 'population_prevalence': '', 'sample_prevalence': '', 'extra': 'keep'})
            path = self.compile(root, original, ('D', 'C', 'B', 'A'), source_manifest=source_manifest)
            selected = pd.read_csv(path, float_precision='round_trip')
            self.assertEqual(len(selected), 9)
            self.assertEqual(selected.p1.tolist(), ['D'] * 3 + ['C'] * 3 + ['A'] * 3)
            self.assertEqual(pd.read_csv(root / 'LDSC_Retained_Traits.csv').traitname.tolist(), ['D', 'C', 'A'])
            self.assertEqual(pd.read_csv(root / 'LDSC_Retained_Traits.csv').note.tolist(), ['d', 'c', 'a'])
            self.assertEqual(pd.read_csv(root / 'LDSC_Dropped_Traits.csv').Trait.tolist(), ['B'])
            self.assertEqual(len(pd.read_csv(root / 'ldsc_results_diagnostic.csv')), 16)
            for row in selected.to_dict('records'):
                expected = record(row['p1'], row['p2'])
                for column in expected.keys() - {'p1', 'p2'}:
                    self.assertEqual(row[column], expected[column], column)
            status = pd.read_csv(root / 'LDSC_Compilation_Status.csv').iloc[0]
            self.assertEqual(status.Status, 'compiled_with_trait_exclusions')
            self.assertEqual(status.Retained_Traits, 3)

    def test_greedy_priority_and_reverse_coverage(self):
        expected = {'valid': [], 'tie': ['C'], 'lower_z': ['B'], 'degree': ['A'],
                    'missing_pair': ['C'], 'missing_self': ['B'], 'missing_reverse': [],
                    'low_signal': [], 'out_of_range': []}
        for kind, excluded in expected.items():
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); self.compile(root, fixture(kind))
                self.assertEqual(pd.read_csv(root / 'LDSC_Dropped_Traits.csv').Trait.tolist(), excluded)
                self.assertEqual(pd.read_csv(root / 'LDSC_Retained_Traits.csv').traitname.tolist(),
                                 [t for t in TRAITS if t not in excluded])
                if kind == 'missing_reverse':
                    self.assertEqual(len(pd.read_csv(root / 'ldsc_results.csv')), 15)
                    self.assertEqual(pd.read_csv(root / 'LDSC_Pair_Status.csv').Status.eq('missing_result').sum(), 1)

    def test_same_h2_z_tie_preserves_earlier_manifest_trait(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); self.compile(root, fixture('tie'), ('A', 'C', 'B', 'D'))
            self.assertEqual(pd.read_csv(root / 'LDSC_Dropped_Traits.csv').Trait.tolist(), ['B'])

    def test_all_failed_or_only_one_remaining_saves_audits_then_fails(self):
        for kind, names in (('all_failed', TRAITS), ('negative_self', ('A', 'B'))):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); rows = fixture(kind)
                rows = rows[rows.p1.map(results._trait_name).isin(names) & rows.p2.map(results._trait_name).isin(names)]
                with self.assertRaisesRegex(results.EstimationFailure, 'at least two'):
                    self.compile(root, rows, names)
                self.assertFalse((root / 'ldsc_results.csv').exists())
                self.assertTrue((root / 'LDSC_Dropped_Traits.csv').exists())
                self.assertEqual(pd.read_csv(root / 'LDSC_Compilation_Status.csv').Status.iloc[0], 'insufficient_traits')

    def test_zero_successful_batches_still_produces_missing_pair_and_exclusion_audits(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with self.assertRaisesRegex(results.EstimationFailure, '0 trait'):
                results.compile_results(root, [], manifest(TRAITS), result_failure_action='drop_traits')
            self.assertEqual(len(pd.read_csv(root / 'LDSC_Pair_Status.csv')), 16)
            self.assertEqual(len(pd.read_csv(root / 'LDSC_Dropped_Traits.csv')), 4)
            self.assertFalse((root / 'ldsc_results.csv').exists())

    def test_conflicting_duplicates_finite_missing_and_large_z_are_fatal_before_drop(self):
        for field, change in (('rg', .2), ('z', .02), ('rg', float('nan'))):
            with self.subTest(field=field, change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); rows = fixture()
                mask = rows.p1.map(results._trait_name).eq('A') & rows.p2.map(results._trait_name).eq('B')
                rows.loc[mask, field] += change
                with self.assertRaisesRegex(RuntimeError, 'Conflicting duplicate'):
                    self.compile(root, rows)
                self.assertFalse((root / 'ldsc_results.csv').exists())
                self.assertEqual(pd.read_csv(root / 'LDSC_Compilation_Status.csv').Status.iloc[0], 'structural_failure')

    def test_consistent_duplicate_orientations_remain_accepted(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); rows = fixture()
            path = self.compile(root, pd.concat([rows, rows.iloc[:1]], ignore_index=True))
            self.assertEqual(len(pd.read_csv(path)), 16)

    def test_structural_input_errors_are_not_silently_dropped(self):
        for field, value in (('rg', 'not_numeric'), ('p1', '/inputs/UNKNOWN.sumstats.gz')):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); rows = fixture(); rows[field] = value
                with self.assertRaises(RuntimeError): self.compile(root, rows)
                self.assertFalse((root / 'ldsc_results.csv').exists())
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(RuntimeError, 'requires ref=yes'):
                self.compile(Path(folder), fixture(), metadata=manifest(TRAITS, ['yes', 'no', 'yes', 'yes']))

    def test_mixed_scales_keep_unmodified_values_and_correct_tie_break(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); rows = fixture('lower_z'); mask = rows.p2.map(results._trait_name).eq('B')
            rows.loc[mask, 'h2_liab'] = rows.loc[mask, 'h2_obs']
            rows.loc[mask, 'h2_liab_se'] = rows.loc[mask, 'h2_obs_se']
            rows.loc[mask, ['h2_obs', 'h2_obs_se']] = float('nan')
            metadata = manifest(TRAITS); metadata.loc[1, 'pop_prevalence'] = .01
            self.compile(root, rows, metadata=metadata)
            self.assertEqual(pd.read_csv(root / 'LDSC_Dropped_Traits.csv').Trait.tolist(), ['B'])

    def test_stale_selected_manifest_is_archived_on_strict_failed_recompile(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); self.compile(root, fixture('negative_self'))
            prior = (root / 'LDSC_Retained_Traits.csv').read_bytes()
            with self.assertRaises(results.EstimationFailure):
                results.compile_results(root, [root / 'batch.csv'], manifest(TRAITS))
            self.assertFalse((root / 'LDSC_Retained_Traits.csv').exists())
            self.assertEqual(next(root.glob('LDSC_Retained_Traits.csv.previous-*')).read_bytes(), prior)

    def test_selection_matches_existing_r_policy_and_retained_output_passes_strict_matrix_qc(self):
        executable = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable: self.skipTest('Rscript unavailable')
        for kind in ('negative_self', 'tie', 'lower_z', 'degree', 'missing_pair', 'missing_reverse', 'missing_self'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); self.compile(root, fixture(kind))
                write_results_csv(pd.DataFrame({'traitname': TRAITS}), root / 'original_manifest.csv')
                result = subprocess.run([executable, str(ROOT / 'tests/test_ldsc_drop_traits.R'), str(ROOT), folder],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class DropTraitsCLITests(unittest.TestCase):
    def run_case(self, native_failure=False, malformed=False):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); inputs = root / 'inputs'; inputs.mkdir()
            ld = root / 'ld'; ld.mkdir(); (ld / '1.l2.ldscore.gz').write_text('reference fixture')
            (ld / '1.l2.M_5_50').write_text('10')
            names = ('A', 'B', 'C'); calls = []
            for name in names:
                with gzip.open(inputs / f'{name}.sumstats.gz', 'wt') as handle:
                    handle.write('SNP A1 A2 N Z\n' + ''.join(f'rs{i} A G 1000 {i / 10}\n' for i in range(10)))
            source = root / 'manifest.csv'
            source.write_text('traitname,ref,population_prevalence,sample_prevalence,note\nA,yes,,,a\nB,yes,,,b\nC,yes,,,c\n')
            def native(cmd, *args, **kwargs):
                tokens = shlex.split(cmd); paths = tokens[tokens.index('--rg') + 1].split(',')
                a, b = map(results._trait_name, paths); calls.append((a, b))
                if native_failure and 'B' in (a, b): raise RuntimeError('simulated batch interruption')
                row = record(a, b); row.update(p1=paths[0], p2=paths[1])
                if 'B' in (a, b): row.update(rg=float('nan'), se=float('nan'), z=float('nan'), p=float('nan'))
                if b == 'B': row['h2_obs'] = -.03
                if malformed and a == b == 'B': row['rg'] = 'invalid'
                write_results_csv(pd.DataFrame([row]), tokens[-1] + '.results.csv')
            args = ['--input', str(source), '--outdir', str(root / 'out'), '--ld_ref', str(ld),
                    '--ldsc_only', '--munged_dir', str(inputs), '--n_cores', '9', '--ldsc_retries', '0',
                    '--chisq_max', '80', '--result_failure_action', 'drop_traits', '--restart']
            with patch.object(ldsc_cli, 'check_runtime'), \
                 patch.object(restart, 'fingerprint_runtime', return_value={'sources': {'ldsc.py': 'fixture'}}), \
                 patch.object(pairwise, 'run_command', side_effect=native), contextlib.redirect_stdout(io.StringIO()):
                if malformed:
                    with self.assertRaises(pairwise.LDSCBatchFailures) as failure: ldsc_cli.main(args)
                    self.assertFalse(failure.exception.can_exclude)
                    return
                ldsc_cli.main(args)
                retained = pd.read_csv(root / 'out/LDSC_Retained_Traits.csv')
                self.assertEqual(retained.traitname.tolist(), ['A', 'C'])
                self.assertEqual(retained.note.tolist(), ['a', 'c'])
                self.assertEqual(len(pd.read_csv(root / 'out/ldsc_results.csv')), 4)
                if native_failure:
                    audit = pd.read_csv(root / 'out/LDSC_Pair_Status.csv')
                    self.assertEqual(audit.Status.eq('missing_result').sum(), 5)
                    self.assertTrue(audit.loc[audit.Status.eq('missing_result'), 'Reason'].str.contains('interruption').all())
                count = len(calls); ldsc_cli.main(args)
                self.assertEqual(len(calls) - count, 5 if native_failure else 0)

    def test_cli_drops_numerical_failure_and_restart_reuses_all_completed_estimates(self):
        self.run_case()

    def test_cli_continues_after_exhausted_execution_retries_and_keeps_successful_batches(self):
        self.run_case(native_failure=True)

    def test_cli_keeps_structural_result_errors_fatal(self):
        self.run_case(malformed=True)


if __name__ == '__main__': unittest.main()
