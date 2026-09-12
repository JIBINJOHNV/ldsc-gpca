import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import concurrent.futures
import pandas as pd
import shlex
import subprocess
real_run = subprocess.run

from ldsc_gpca import extraction, munging, pairwise, results, utils

class PrevalenceTests(unittest.TestCase):
    def worker(self, cases=('100', '200'), controls=('900', '800'), supplied=None):
        def extract(command, stdout, **kwargs):
            self.assertIn('pipefail', command)
            raw = ''.join(f'1\trs{i}\t100\tA\tG\t2\t8\t0.2\t360\t{a}\t{b}\n' for i, (a,b) in enumerate(zip(cases,controls)))
            tokens = shlex.split(command[-1])
            return real_run(tokens[tokens.index('awk'):], input=raw, stdout=stdout, stderr=subprocess.PIPE, text=True, check=True)
        with tempfile.TemporaryDirectory() as folder, patch.object(extraction.subprocess, 'run', side_effect=extract):
            result = extraction.munge_input_worker_with_prevalence('trait', folder, folder, folder+'/input.vcf', supplied)
            table = pd.read_csv(Path(folder)/'munge_input/trait_mungeinput.tsv', sep='\t')
            self.assertEqual(table.N_TOTAL.tolist(), [float(a)+float(b) for a,b in zip(cases,controls)])
            self.assertTrue((table.P == 1e-8).all())
            return result

    def test_calculated(self):
        self.assertAlmostEqual(self.worker(), .15)

    def test_supplied(self):
        self.assertEqual(self.worker(supplied=.3), .3)

    def test_constant(self):
        self.assertEqual(self.worker(cases=('100','100'), controls=('900','900')), .1)

    def test_invalid_counts(self):
        for counts in [('.', '.'), ('100', '.'), ('0','100'), ('-1','100'), ('inf','100')]:
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                self.worker(cases=counts)

    def test_prevalence_validation(self):
        for value in [0,1,-.1,float('inf'),'bad']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                utils.optional_prevalence(value,'test')
        for value in [None, float('nan'), 'NA', '.']:
            self.assertIsNone(utils.optional_prevalence(value,'test'))

    def test_routing(self):
        frame = pd.DataFrame({'gwas_name':['q','b'], 'vcf_files':['/x/q','/x/b'],
                              'pop_prevalence':[float('nan'),.01], 'sample_prevalence':[float('nan'),float('nan')]})
        with tempfile.TemporaryDirectory() as folder, patch.object(extraction.concurrent.futures, 'ProcessPoolExecutor', concurrent.futures.ThreadPoolExecutor), patch.object(extraction,'munge_input_worker', return_value=None) as old, patch.object(extraction,'munge_input_worker_with_prevalence',return_value=.2) as new:
            result = extraction.run_vcf_to_table(1,frame,folder)
            old.assert_called_once()
            new.assert_called_once()
            self.assertEqual(result.loc[1,'sample_prevalence'],.2)

    def test_pairwise_flags_and_failures(self):
        frame = pd.DataFrame({'gwas_name':['q','b'], 'ref':['yes','no'],
                              'pop_prevalence':[float('nan'),.01], 'sample_prevalence':[float('nan'),.2]})
        with tempfile.TemporaryDirectory() as folder, patch.object(pairwise,'run_command') as run:
            pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder)
            commands = [call.args[0] for call in run.call_args_list]
            self.assertNotIn('--samp-prev', commands[0])
            self.assertIn('--samp-prev nan,0.2 --pop-prev nan,0.01',commands[1])
        with tempfile.TemporaryDirectory() as folder, patch.object(pairwise,'run_command',side_effect=RuntimeError('failed')), self.assertRaises(RuntimeError):
            pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder)

    def test_compile_mixed_headers(self):
        with tempfile.TemporaryDirectory() as folder:
            paths=[]
            for i, (target, scale) in enumerate([('q','obs'),('q','liab'),('b','liab')]):
                path=Path(folder)/f'{i}.log'
                path.write_text(f'p1 p2 rg se z p h2_{scale} h2_{scale}_se gcov_int_se\nq {target} 1 .1 10 .01 .2 .03 .01\n\n')
                paths.append(str(path))
            results.compile_results(folder,paths,pd.DataFrame({'gwas_name':['q','b'],'ref':['yes','no'],'pop_prevalence':[float('nan'),.01]}))
            result=pd.read_csv(Path(folder)/'ldsc_results.csv')
            self.assertTrue(result.loc[result.p2=='q','h2_liab'].isna().all())
            self.assertTrue((result.loc[result.p2=='q','h2_obs']==.2).all())
            self.assertTrue((result.loc[result.p2=='b','h2_liab']==.2).all())

    def test_munge_selects_n_and_records_provenance(self):
        import gzip
        with tempfile.TemporaryDirectory() as folder:
            frame=pd.DataFrame({'gwas_name':['q','b'], 'sample_size_column':['NEF','N_TOTAL'],
                                'sample_prevalence':[float('nan'),.2]})
            def run(command, label, **kwargs):
                name=label.removeprefix('Munge_')
                if name=='q':
                    self.assertIn('--N-col NEF',command)
                    self.assertNotIn('--ignore',command)
                else:
                    self.assertIn('--N-col N_TOTAL --ignore N_CASES,N_CONTROLS,NEF',command)
                with gzip.open(Path(folder)/(name+'.sumstats.gz'),'wt') as handle:
                    handle.write('SNP A1 A2 Z N\n' + '\n'.join(f'rs{i} A G 2 1000' for i in range(50)))
            with patch.object(munging,'run_command',side_effect=run):
                self.assertEqual(set(munging.parallel_munge_sumstats(1,frame,folder,folder+'/hm3', ldsc_input_folder=folder, munge_input_folder=folder)),{'q','b'})
            import json
            self.assertEqual(json.loads((Path(folder)/'b.prevalence.json').read_text())['sample_size_column'],'N_TOTAL')
            self.assertEqual(json.loads((Path(folder)/'b.prevalence.json').read_text())['filters'], munging.DEFAULT_FILTERS)

    def test_original_worker_shell_and_failure(self):
        def extract(command, stdout, **kwargs):
            self.assertIn('pipefail', command)
            tokens = shlex.split(command[-1])
            self.assertIn('(ABS(INFO/AF - INFO/EUR) > 0.2 || INFO/EUR==".")', tokens)
            raw = '1\trs1\t100\tA\tG\t2\t8\t0.2\t360\n'
            return real_run(tokens[tokens.index('awk'):], input=raw, stdout=stdout, stderr=subprocess.PIPE, text=True, check=True)
        with tempfile.TemporaryDirectory(prefix='ldsc space ') as folder:
            with patch.object(extraction.subprocess, 'run', side_effect=extract):
                extraction.munge_input_worker('trait', folder, folder, folder+'/input file.vcf')
            table = pd.read_csv(Path(folder)/'munge_input/trait_mungeinput.tsv', sep=r'\s+')
            self.assertEqual(table.P.tolist(), [1e-8])
            def fail(command, **kwargs):
                return real_run(['bash', '-o', 'pipefail', '-c', 'false | awk "{print}"'], **kwargs)
            with patch.object(extraction.subprocess, 'run', side_effect=fail), self.assertRaises(subprocess.CalledProcessError):
                extraction.munge_input_worker('trait', folder, folder, folder+'/input.vcf')

    def test_incomplete_and_invalid_results(self):
        metadata = pd.DataFrame({'gwas_name':['q','b'], 'ref':['yes','no'], 'pop_prevalence':[float('nan'),.01]})
        header = 'p1 p2 rg se z p h2_obs h2_obs_se gcov_int_se\n'
        row = 'q q 1 .1 10 .01 .2 .03 .01\n'
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'test.log'
            for content, error in [('', 'No correlation'), (header+row, 'missing='),
                                   (header+'q q 1\n', 'Malformed'),
                                   (header+row+'q b NA .1 10 .01 .2 .03 .01\n', 'Non-finite')]:
                path.write_text(content)
                with self.subTest(error=error), self.assertRaisesRegex(RuntimeError, error):
                    results.compile_results(folder, [str(path)], metadata)
                self.assertFalse((Path(folder)/'ldsc_results.csv').exists())

    def test_saved_filters(self):
        results.check_saved_filters({'filters':dict(utils.DEFAULT_FILTERS)}, utils.DEFAULT_FILTERS, 'q')
        with self.assertRaisesRegex(ValueError, 'filter settings differ'):
            results.check_saved_filters({'filters':{**utils.DEFAULT_FILTERS,'maf_min':.02}},utils.DEFAULT_FILTERS,'q')
        with patch('builtins.print') as log:
            results.check_saved_filters({}, utils.DEFAULT_FILTERS, 'q')
            self.assertIn('unknown', log.call_args.args[0])

    def test_custom_munge_maf(self):
        import gzip
        with tempfile.TemporaryDirectory(prefix='ldsc space ') as folder:
            frame=pd.DataFrame({'gwas_name':['q'], 'sample_size_column':['NEF'], 'sample_prevalence':[float('nan')]})
            def run(command, label, **kwargs):
                tokens = shlex.split(command)
                self.assertEqual(tokens[tokens.index('--maf-min')+1], '0.02')
                self.assertEqual(tokens[tokens.index('--out')+1], folder+'/q')
                with gzip.open(Path(folder)/'q.sumstats.gz','wt') as handle:
                    handle.write('SNP A1 A2 Z N\n'+'\n'.join(f'rs{i} A G 2 1000' for i in range(50)))
            with patch.object(munging, 'run_command', side_effect=run):
                self.assertEqual(munging.parallel_munge_sumstats(1,frame,folder,folder+'/hm3',{'munge_maf_min':.02}, ldsc_input_folder=folder, munge_input_folder=folder), ['q'])

    def test_retries_only_failed_batch(self):
        frame = pd.DataFrame({'gwas_name':['q','b'], 'ref':['yes','no'],
                              'pop_prevalence':[float('nan'),float('nan')],
                              'sample_prevalence':[float('nan'),float('nan')]})
        attempts = {}
        def run(command, label, **kwargs):
            attempts[command] = attempts.get(command, 0) + 1
            if 'batch_1 ' in label and attempts[command] == 1:
                raise RuntimeError('temporary failure')
        with tempfile.TemporaryDirectory() as folder, patch.object(pairwise, 'run_command', side_effect=run):
            logs = pairwise.parallel_ldsc_analysis(2,1,folder,folder,frame,folder)
            self.assertEqual(len(logs),2)
            self.assertEqual(sorted(attempts.values()),[1,2])

    def test_retry_limit_and_validation(self):
        frame = pd.DataFrame({'gwas_name':['q'], 'ref':['yes'],
                              'pop_prevalence':[float('nan')], 'sample_prevalence':[float('nan')]})
        with tempfile.TemporaryDirectory() as folder:
            for retries in (0,1,2):
                with self.subTest(retries=retries), patch.object(pairwise,'run_command',side_effect=RuntimeError('failed')) as run:
                    with self.assertRaisesRegex(RuntimeError,'retries exhausted'):
                        pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder,retries=retries)
                    self.assertEqual(run.call_count,retries+1)
            for retries in (-1,1.5,None):
                with self.subTest(retries=retries), self.assertRaises(ValueError):
                    pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder,retries=retries)

if __name__ == '__main__':
    unittest.main()
