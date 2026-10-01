"""Numerical failures remain auditable; structural/execution failures remain fatal."""
import contextlib
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
from ldsc_gpca import results, pairwise, ldsc_cli
from ldsc_gpca.ldsc_export import write_results_csv
from test_results_csv import record, manifest

ROOT = Path(__file__).resolve().parents[1]


def failing_records():
    rows = [record(a, b) for a in ('A', 'B', 'C') for b in ('A', 'B', 'C')]
    for r in rows:
        if r['p1'].endswith('/A.sumstats.gz') or r['p2'].endswith('/A.sumstats.gz'):
            r.update(rg=float('nan'), se=float('nan'), z=float('nan'), p=float('nan'))
        if r['p2'].endswith('/A.sumstats.gz'):
            r['h2_obs'] = -0.03
    return rows


class FailureHandlingTests(unittest.TestCase):
    def test_strict_saves_all_estimates_and_report_preserves_valid_values(self):
        for policy in ('error', 'report'):
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                original = pd.DataFrame(failing_records())
                source = root / 'batch.results.csv'
                write_results_csv(original, source)
                with contextlib.redirect_stdout(io.StringIO()):
                    if policy == 'error':
                        with self.assertRaises(results.EstimationFailure):
                            results.compile_results(root, [source], manifest(('C', 'A', 'B')))
                    else:
                        path = results.compile_results(root, [source], manifest(('C', 'A', 'B')),
                                                       result_failure_action=policy)
                        self.assertEqual(Path(path).name, 'ldsc_results_diagnostic.csv')
                self.assertFalse((root / 'ldsc_results.csv').exists())
                observed = pd.read_csv(root / 'ldsc_results_diagnostic.csv', float_precision='round_trip')
                self.assertEqual(len(observed), 9)
                self.assertEqual(observed.p1.tolist(), ['C']*3 + ['A']*3 + ['B']*3)
                self.assertEqual(observed.rg.isna().sum(), 5)
                for a in ('B', 'C'):
                    for b in ('B', 'C'):
                        row = observed[(observed.p1 == a) & (observed.p2 == b)].iloc[0]
                        for c, value in record(a, b).items():
                            if c not in ('p1', 'p2'):
                                self.assertEqual(row[c], value, (a, b, c))
                status = pd.read_csv(root / 'LDSC_Pair_Status.csv')
                self.assertEqual(status.Status.eq('failed_estimate').sum(), 5)
                self.assertEqual(set(status.Source_Row), set(range(2, 11)))
                traits = pd.read_csv(root / 'LDSC_Trait_Status.csv', float_precision='round_trip')
                self.assertEqual(traits.Trait.tolist(), ['C', 'A', 'B'])
                self.assertEqual(traits.Self_Status.tolist(), ['valid', 'failed_estimate', 'valid'])
                self.assertEqual(traits.loc[1, 'h2_obs'], -0.03)

    def test_extra_columns_and_unused_empty_scale_are_ignored(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / 'batch.csv'
            r = record(); r.update(annotation='anything', h2_liab='NA', h2_liab_se='NA')
            write_results_csv(pd.DataFrame([r])[list(reversed(r))], source)
            path = results.compile_results(folder, [source], manifest())
            data = pd.read_csv(path, float_precision='round_trip')
            self.assertEqual(data.rg.iloc[0], r['rg'])
            self.assertNotIn('annotation', data)

    def test_low_signal_and_out_of_range_rg_warn_without_exclusion(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / 'batch.csv'
            rows = [record(a,b) for a in ('A','B') for b in ('A','B')]
            for r in rows:
                r['h2_obs'] = 0.001
                if r['p1'] != r['p2']: r['rg'] = 1.1
            write_results_csv(pd.DataFrame(rows), source)
            results.compile_results(folder, [source], manifest(('A','B')))
            status = pd.read_csv(Path(folder) / 'LDSC_Pair_Status.csv')
            self.assertTrue(status.Status.eq('valid').all())
            self.assertTrue(status.Warning.notna().all())

    def test_mixed_scale_rows_ignore_only_the_unused_alternative(self):
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'batch.csv'
            a=record('A','A'); b=record('A','B')
            b['h2_liab']=b.pop('h2_obs'); b['h2_liab_se']=b.pop('h2_obs_se')
            write_results_csv(pd.DataFrame([a,b]),source)
            meta=manifest(('A','B'),['yes','no']); meta.loc[1,'pop_prevalence']=.01
            results.compile_results(folder,[source],meta)
            status=pd.read_csv(Path(folder)/'LDSC_Pair_Status.csv')
            self.assertTrue(status.Status.eq('valid').all())

    def test_report_never_hides_structural_failures(self):
        for modify, message in (
            (lambda d: d.assign(rg='not_a_number'), 'Non-numeric'),
            (lambda d: d.drop(columns='se'), 'Missing LDSC result columns'),
            (lambda d: d.assign(p1='/input/UNKNOWN.sumstats.gz'), 'unexpected='),
        ):
            with self.subTest(message=message), tempfile.TemporaryDirectory() as folder:
                source = Path(folder) / 'batch.csv'
                write_results_csv(modify(pd.DataFrame([record()])), source)
                with self.assertRaisesRegex(RuntimeError, message):
                    results.compile_results(folder, [source], manifest(), result_failure_action='report')
                self.assertEqual(pd.read_csv(Path(folder)/'LDSC_Compilation_Status.csv').Status.iloc[0], 'structural_failure')

    def test_stale_success_is_archived_on_failed_rerun(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); source=root/'batch.csv'
            write_results_csv(pd.DataFrame([record()]), source)
            results.compile_results(root, [source], manifest())
            previous=(root/'ldsc_results.csv').read_bytes()
            write_results_csv(pd.DataFrame([dict(record(), rg=float('nan'))]), source)
            with self.assertRaises(results.EstimationFailure):
                results.compile_results(root, [source], manifest())
            self.assertFalse((root/'ldsc_results.csv').exists())
            self.assertEqual(next(root.glob('ldsc_results.csv.previous-*')).read_bytes(), previous)

    def test_all_missing_h2_and_all_failed_rows_stay_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'batch.csv'
            r=record(); r.update(h2_obs=float('nan'), h2_obs_se=float('nan'), rg=float('nan'))
            write_results_csv(pd.DataFrame([r]), source)
            results.compile_results(folder,[source],manifest(),result_failure_action='report')
            self.assertFalse((Path(folder)/'ldsc_results.csv').exists())
            observed=pd.read_csv(Path(folder)/'ldsc_results_diagnostic.csv')
            self.assertTrue(observed[['h2_obs','h2_obs_se','rg']].isna().all().all())

    def test_cli_option_and_extra_alias_columns_do_not_override_canonical(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); source=root/'manifest.csv'
            source.write_text('traitname,ref,population_prevalence,sample_prevalence,gwas_name,pop_prevalence\nA,yes,,,WRONG,0.9\n')
            with patch.object(ldsc_cli,'check_runtime'), patch.object(ldsc_cli,'is_valid_gz',return_value=True), \
                 patch.object(ldsc_cli,'check_saved_filters'), patch.object(ldsc_cli,'parallel_ldsc_analysis',return_value=[]) as run, \
                 patch.object(ldsc_cli,'compile_results') as compile:
                ldsc_cli.main(['--input',str(source),'--outdir',str(root/'out'),'--ld_ref',str(root),
                               '--ldsc_only','--result_failure_action','report'])
            self.assertEqual(run.call_args.args[4].gwas_name.tolist(), ['A'])
            self.assertTrue(run.call_args.args[4].pop_prevalence.isna().all())
            self.assertEqual(compile.call_args.kwargs['result_failure_action'], 'report')


class BatchFailureTests(unittest.TestCase):
    def test_failed_job_does_not_discard_independent_success(self):
        with tempfile.TemporaryDirectory() as folder:
            frame=manifest(('A','B')); frame['sample_prevalence']=float('nan')
            def command(cmd,*a,**k):
                output=Path(shlex.split(cmd)[-1]+'.results.csv')
                if output.name.startswith('B_'): raise RuntimeError('interrupted')
                output.write_text('successful result')
            with patch.object(pairwise,'run_command',side_effect=command), \
                 self.assertRaisesRegex(RuntimeError,'completed batch files were retained'):
                pairwise.parallel_ldsc_analysis(2,100,folder,folder,frame,folder,retries=0)
            self.assertEqual(Path(folder,'A_100_ldsc_0.results.csv').read_text(),'successful result')
            self.assertEqual(pd.read_csv(Path(folder)/'LDSC_Batch_Status.csv').Status.tolist(),
                             ['completed','execution_failed'])

    def run_batch(self, folder, command):
        frame=manifest(); frame['sample_prevalence']=float('nan')
        with patch.object(pairwise, 'run_command', side_effect=command) as run:
            value=pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder,retries=1)
        return value, run.call_count

    def test_missing_output_is_retried_then_succeeds(self):
        attempts=[]
        def command(cmd,*a,**k):
            attempts.append(cmd)
            if len(attempts)==2: Path(shlex.split(cmd)[-1]+'.results.csv').write_text('fresh')
        with tempfile.TemporaryDirectory() as folder:
            _,count=self.run_batch(folder,command)
            self.assertEqual(count,2)
            self.assertEqual(pd.read_csv(Path(folder)/'LDSC_Batch_Status.csv').Status.iloc[0],'completed')

    def test_execution_error_and_stale_output_fail_after_bounded_retries(self):
        for stale in (False,True):
            with tempfile.TemporaryDirectory() as folder:
                if stale: Path(folder,'A_1_ldsc_0.results.csv').write_text('old')
                command=None if stale else RuntimeError('interrupted')
                with self.assertRaisesRegex(RuntimeError,'retries exhausted'):
                    self.run_batch(folder,command)
                status=pd.read_csv(Path(folder)/'LDSC_Batch_Status.csv')
                self.assertEqual(status.Status.iloc[0],'execution_failed')


class RFailureTests(unittest.TestCase):
    def test_compiled_diagnostic_table_reaches_strict_and_drop_policies(self):
        executable=os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable: self.skipTest('Rscript unavailable')
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); source=root/'batch.csv'
            write_results_csv(pd.DataFrame(failing_records()),source)
            diagnostic=results.compile_results(root,[source],manifest(('C','A','B')),result_failure_action='report')
            (root/'manifest.csv').write_text('traitname,extra,gwas_name\nC,note,ignored\nA,note,ignored\nB,note,ignored\n')
            for policy, matrix in (('error','correlation'), ('drop_traits','correlation'), ('drop_traits','covariance')):
                out=root/(policy+'_'+matrix)
                outcome=subprocess.run([executable,str(ROOT/'src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r'),
                    '--input',str(root/'manifest.csv'),'--ldsc_results',str(diagnostic),'--outdir',str(out),
                    '--validate_only','--failed_ldsc_action',policy,'--pca_matrix',matrix],capture_output=True,text=True)
                status=pd.read_csv(out/'GenomicPCA_Run_Status.csv')
                if policy=='error':
                    self.assertNotEqual(outcome.returncode,0)
                    self.assertEqual(status.Status.iloc[0],'failed')
                    self.assertEqual(pd.read_csv(out/'GWAMA_Run_Status.csv').Chromosome.iloc[0],'not_run')
                else:
                    self.assertEqual(outcome.returncode,0,outcome.stdout+outcome.stderr)
                    self.assertEqual(status.Status.iloc[0],'validated')
                    self.assertEqual(pd.read_csv(out/'Python_LDSC_Retained_Traits.csv').traitname.tolist(),['C','B'])
                    self.assertEqual(pd.read_csv(out/'Python_LDSC_Dropped_Failed_Traits.csv').Trait.tolist(),['A'])
                self.assertTrue((out/'Python_LDSC_Self_Pair_QC.csv').exists())

    def test_matrix_and_input_acceptance(self):
        executable=os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not executable: self.skipTest('Rscript unavailable')
        outcome=subprocess.run([executable,str(ROOT/'tests/test_failure_handling.R'),str(ROOT)],capture_output=True,text=True)
        self.assertEqual(outcome.returncode,0,outcome.stdout+outcome.stderr)


if __name__ == '__main__': unittest.main()
