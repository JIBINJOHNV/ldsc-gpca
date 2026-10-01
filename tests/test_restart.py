"""Restart must reuse only verified, complete results with matching provenance."""
import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import shlex
import tempfile
import types
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import pairwise, restart, ldsc_cli, ldsc_runtime, ldsc_provenance
from ldsc_gpca.ldsc_export import write_results_csv
from test_results_csv import record, manifest


class RestartTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root / 'out'; self.out.mkdir()
        self.inputs = self.root / 'sumstats'; self.inputs.mkdir()
        self.ld = self.root / 'ld'; self.ld.mkdir()
        self.weights = self.root / 'weights'; self.weights.mkdir()
        for folder in (self.ld, self.weights):
            (folder / '1.l2.ldscore.gz').write_bytes(b'LD score fixture')
            (folder / '1.l2.M_5_50').write_text('100\n')
        self.meta = manifest(('A', 'B')); self.meta['sample_prevalence'] = float('nan')
        self.parameters = {'chisq_max': 80, 'filters': {'maf_min': .01}}
        self.runtime = {'sources': {'ldsc.py': 'v1'}, 'packages': {'numpy': '1.26.4'}}
        self.calls = []
        for name in ('A', 'B'):
            self.write_sumstats(name)

    def write_sumstats(self, name, text='SNP A1 A2 N Z\nrs1 A G 1000 1.5\n', mtime=1):
        (self.inputs / f'{name}.sumstats.gz').write_bytes(gzip.compress(text.encode(), mtime=mtime))

    def command(self, cmd, *args, **kwargs):
        tokens = shlex.split(cmd)
        paths = tokens[tokens.index('--rg') + 1].split(',')
        self.calls.append(paths)
        rows = []
        for target in paths[1:]:
            row = record(restart._trait_name(paths[0]), restart._trait_name(target))
            row.update(p1=paths[0], p2=target)
            rows.append(row)
        write_results_csv(pd.DataFrame(rows), tokens[-1] + '.results.csv')

    def run_batches(self, *, restart_enabled=True, command=None, batch_size=1, workers=2):
        with contextlib.redirect_stdout(io.StringIO()), \
             patch.object(restart, 'fingerprint_runtime', return_value=self.runtime), \
             patch.object(pairwise, 'run_command', side_effect=command or self.command):
            return pairwise.parallel_ldsc_analysis(workers, batch_size, str(self.out), str(self.ld),
                self.meta, str(self.inputs), retries=0, ld_weights_dir=str(self.weights),
                restart=restart_enabled, checkpoint_parameters=self.parameters)

    def status(self):
        return pd.read_csv(self.out / 'LDSC_Batch_Status.csv')

    def test_identical_restart_runs_zero_commands_and_preserves_order(self):
        first = self.run_batches(restart_enabled=False)
        self.assertEqual(len(self.calls), 4)
        second = self.run_batches(workers=1)
        self.assertEqual(first, second)
        self.assertEqual(len(self.calls), 4)
        self.assertEqual(self.status().Status.tolist(), ['reused'] * 4)
        self.assertTrue(self.status().Restart_Reason.eq('verified_completed_batch').all())
        self.run_batches(restart_enabled=False)
        self.assertEqual(len(self.calls), 8)

    def test_interruption_keeps_other_checkpoints_and_reruns_only_failed_batch(self):
        def fail_one(cmd, *args, **kwargs):
            if shlex.split(cmd)[-1].endswith('B_1_ldsc_1'):
                # Even a full CSV from an unsuccessful process cannot be reused.
                self.command(cmd)
                raise RuntimeError('simulated interruption after export')
            self.command(cmd)
        with self.assertRaisesRegex(RuntimeError, 'simulated interruption'):
            self.run_batches(command=fail_one)
        self.assertEqual(len(list(self.out.glob('*.checkpoint.json'))), 3)
        self.calls.clear()
        self.run_batches()
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(self.status().Status.tolist(), ['reused'] * 3 + ['completed'])

    def test_missing_correlations_remain_diagnostic_without_repeated_regression(self):
        def missing_rg(cmd, *args, **kwargs):
            self.command(cmd)
            path = shlex.split(cmd)[-1] + '.results.csv'
            frame = pd.read_csv(path); frame['rg'] = float('nan')
            write_results_csv(frame, path)
        paths = self.run_batches(command=missing_rg)
        self.run_batches()
        self.assertEqual(len(self.calls), 4)
        self.assertTrue(pd.read_csv(paths[0]).rg.isna().all())

    def test_changed_sumstats_reruns_only_involving_batches_even_with_same_size_mtime(self):
        self.run_batches()
        path = self.inputs / 'A.sumstats.gz'; previous = path.stat()
        self.write_sumstats('A', 'SNP A1 A2 N Z\nrs1 A G 1000 2.5\n')
        os.utime(path, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        self.calls.clear(); self.run_batches()
        self.assertEqual(len(self.calls), 3)
        self.assertEqual(self.status().Status.tolist(), ['completed'] * 3 + ['reused'])

    def test_gzip_header_change_and_extra_manifest_columns_do_not_invalidate(self):
        self.run_batches()
        self.write_sumstats('A', mtime=100)
        self.meta['annotation'] = ['anything', 'else']
        self.calls.clear(); self.run_batches()
        self.assertEqual(self.calls, [])

    def test_reference_weights_parameters_prevalence_runtime_invalidate(self):
        changes = {
            'reference': lambda: (self.ld / '1.l2.ldscore.gz').write_bytes(b'changed reference'),
            'M': lambda: (self.ld / '1.l2.M_5_50').write_text('101\n'),
            'weights': lambda: (self.weights / '1.l2.ldscore.gz').write_bytes(b'changed weights'),
            'threshold': lambda: self.parameters.update(chisq_max=81),
            'filters': lambda: self.parameters['filters'].update(maf_min=.02),
            'runtime': lambda: self.runtime['sources'].update({'ldsc.py': 'v2'}),
            'prevalence': lambda: self.meta.loc.__setitem__((slice(None), 'pop_prevalence'), .01),
            'sample_prevalence': lambda: self.meta.loc.__setitem__((slice(None), 'sample_prevalence'), .2),
        }
        self.run_batches()
        for name, change in changes.items():
            with self.subTest(change=name):
                change(); self.calls.clear(); self.run_batches()
                self.assertEqual(len(self.calls), 4)
                self.assertTrue(self.status().Restart_Reason.eq('inputs_parameters_or_runtime_changed').all())

    def test_changed_batch_layout_never_reuses_incorrect_pairs(self):
        self.run_batches(batch_size=2)
        self.calls.clear(); self.run_batches(batch_size=1)
        self.assertEqual(len(self.calls), 4)
        self.meta = self.meta.iloc[::-1].reset_index(drop=True)
        self.calls.clear(); self.run_batches(batch_size=1)
        self.assertEqual(len(self.calls), 4)

    def test_missing_corrupt_legacy_or_modified_outputs_are_recomputed(self):
        for damage in ('missing_result', 'empty', 'modified', 'missing_checkpoint', 'invalid_json', 'json_list'):
            with self.subTest(damage=damage):
                paths = self.run_batches()
                path = Path(paths[0]); checkpoint = Path(str(path) + restart.CHECKPOINT_SUFFIX)
                if damage == 'missing_result': path.unlink()
                elif damage == 'empty': path.write_text('')
                elif damage == 'modified': path.write_text(path.read_text().replace('0.123', '0.456'))
                elif damage == 'missing_checkpoint': checkpoint.unlink()
                elif damage == 'invalid_json': checkpoint.write_text('{broken')
                else: checkpoint.write_text('[]')
                # Ensure modified output actually changes regardless of fixture decimals.
                if damage == 'modified': path.write_text(path.read_text() + '\n')
                self.calls.clear(); self.run_batches()
                self.assertEqual(len(self.calls), 1)

    def test_wrong_coverage_or_bad_text_cannot_receive_checkpoint(self):
        for bad in ('coverage', 'numeric', 'missing_column'):
            with self.subTest(bad=bad):
                def invalid(cmd, *args, **kwargs):
                    self.command(cmd)
                    path = shlex.split(cmd)[-1] + '.results.csv'
                    frame = pd.read_csv(path)
                    if bad == 'coverage': frame['p2'] = 'UNKNOWN'
                    elif bad == 'numeric': frame['rg'] = 'not_numeric'
                    else: frame = frame.drop(columns='se')
                    write_results_csv(frame, path)
                with self.assertRaisesRegex(RuntimeError, 'retries exhausted'):
                    self.run_batches(restart_enabled=False, command=invalid)
                self.assertEqual(list(self.out.glob('*.checkpoint.json')), [])

    def test_partial_table_rejected_even_when_checkpoint_checksum_matches(self):
        paths = self.run_batches(batch_size=2)
        path = Path(paths[0]); checkpoint = Path(str(path) + restart.CHECKPOINT_SUFFIX)
        write_results_csv(pd.read_csv(path).iloc[:1], path)
        saved = json.loads(checkpoint.read_text()); saved['result_sha256'] = restart.file_digest(path)
        restart.write_checkpoint(checkpoint, saved)
        self.calls.clear(); self.run_batches(batch_size=2)
        self.assertEqual(len(self.calls), 1)

    def test_input_mutation_during_execution_never_gets_checkpoint(self):
        def change_input(cmd, *args, **kwargs):
            self.command(cmd)
            (self.ld / '1.l2.M_5_50').write_text('changed mid run\n')
        with self.assertRaisesRegex(RuntimeError, 'inputs changed during execution'):
            self.run_batches(command=change_input, workers=1)
        self.assertEqual(list(self.out.glob('*.checkpoint.json')), [])

    def test_missing_reference_or_invalid_gzip_stops_before_launch(self):
        (self.ld / '1.l2.ldscore.gz').unlink()
        with self.assertRaisesRegex(RuntimeError, 'No LD-score files'):
            self.run_batches()
        self.assertEqual(self.calls, [])
        (self.ld / '1.l2.ldscore.gz').write_bytes(b'reference')
        (self.inputs / 'A.sumstats.gz').write_bytes(b'invalid gzip')
        with self.assertRaises(OSError): self.run_batches()
        self.assertEqual(self.calls, [])

    def test_cli_restart_with_real_chisq_filter_and_compilation(self):
        for name in ('A', 'B'):
            self.write_sumstats(name, 'SNP A1 A2 N Z\n' + ''.join(
                f'rs{i} A G 1000 {1 + i / 10}\n' for i in range(10)))
        source = self.root / 'manifest.csv'
        source.write_text('traitname,ref,population_prevalence,sample_prevalence,extra\nA,yes,,,x\nB,yes,,,y\n')
        arguments = ['--input', str(source), '--outdir', str(self.out), '--ld_ref', str(self.ld),
                     '--ld_weights', str(self.weights), '--ldsc_only', '--munged_dir', str(self.inputs),
                     '--n_cores', '4', '--chisq_max', '80', '--result_failure_action', 'report']
        with patch.object(ldsc_cli, 'check_runtime'), \
             patch.object(restart, 'fingerprint_runtime', return_value=self.runtime), \
             patch.object(pairwise, 'run_command', side_effect=self.command), \
             contextlib.redirect_stdout(io.StringIO()):
            ldsc_cli.main(arguments)
            first = (self.out / 'ldsc_results.csv').read_bytes()
            ldsc_cli.main([*arguments, '--restart'])
        self.assertEqual(len(self.calls), 4)
        self.assertTrue(pd.read_csv(self.out / 'ldsc_results' / 'LDSC_Batch_Status.csv').Status.eq('reused').all())
        self.assertEqual(first, (self.out / 'ldsc_results.csv').read_bytes())
        self.assertEqual(len(pd.read_csv(self.out / 'ldsc_results.csv')), 4)


class RuntimeAndCLITests(unittest.TestCase):
    def test_runtime_probe_parses_identity_and_rejects_failure(self):
        for code, output in ((0, '{"sources":{"ldsc.py":"abc"}}'), (1, '{}'), (0, 'garbage')):
            with self.subTest(code=code, output=output), patch.object(ldsc_runtime.subprocess, 'run',
                    return_value=types.SimpleNamespace(returncode=code, stdout=output, stderr='problem')) as run:
                if code == 0 and 'sources' in output:
                    self.assertIn('sources', ldsc_runtime.fingerprint_runtime(prefix='/child'))
                    self.assertIn('/child', run.call_args.args[0])
                else:
                    with self.assertRaisesRegex(RuntimeError, 'Cannot verify LDSC runtime'):
                        ldsc_runtime.fingerprint_runtime()

    def test_runtime_identity_changes_when_native_source_changes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); package = root / 'ldscore'; package.mkdir()
            init = package / '__init__.py'; init.write_text('')
            source = package / 'regressions.py'; source.write_text('version=1\n')
            native = root / 'ldsc.py'; native.write_text('native=1\n')
            with patch.dict('sys.modules', {'ldscore': types.SimpleNamespace(__file__=str(init))}), \
                 patch.object(ldsc_provenance.shutil, 'which', return_value=str(native)), \
                 patch.object(ldsc_provenance.importlib.metadata, 'version', return_value='1.0'):
                first = ldsc_provenance.runtime_identity()
                source.write_text('version=2\n')
                second = ldsc_provenance.runtime_identity()
            self.assertNotEqual(first['sources'], second['sources'])
            self.assertIn('ldsc_export.py', first['sources'])

    def test_restart_flag_reaches_scheduler_and_checkpoint_configuration(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); manifest_path = root / 'manifest.csv'
            manifest_path.write_text('traitname,ref,population_prevalence,sample_prevalence\nA,yes,,\n')
            with patch.object(ldsc_cli, 'check_runtime'), patch.object(ldsc_cli, 'is_valid_gz', return_value=True), \
                 patch.object(ldsc_cli, 'check_saved_filters'), patch.object(ldsc_cli, 'filter_munged_sumstats', return_value=folder), \
                 patch.object(ldsc_cli, 'parallel_ldsc_analysis', return_value=[]) as run, \
                 patch.object(ldsc_cli, 'compile_results'), contextlib.redirect_stdout(io.StringIO()):
                ldsc_cli.main(['--input', str(manifest_path), '--outdir', folder, '--ld_ref', folder,
                               '--ldsc_only', '--restart', '--chisq_max', '80'])
            self.assertTrue(run.call_args.kwargs['restart'])
            self.assertEqual(run.call_args.kwargs['checkpoint_parameters']['chisq_max'], 80)


if __name__ == '__main__': unittest.main()
