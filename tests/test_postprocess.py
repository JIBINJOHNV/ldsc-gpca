"""GWAMA export values, sorting, compression and transactional failure checks."""
import contextlib
import csv
import gzip
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

import pandas as pd
import polars as pl
from ldsc_gpca import gpca, postprocess as pp


COLUMNS = pp.GPCA_COLUMNS + ['Direction', 'Extra']


def row(snp='rs01', chrom='1', bp='123', direction='+?-'):
    return dict(zip(COLUMNS, [snp, chrom, bp, 'A', 'G', '0.1234567890123456789',
        '32817.123456789', '-0.0000000000000000001', '0.0012345678901234567',
        '-12.345678901234567', '1.234567890123456e-310', 'NA', direction, 'NA']))


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        redirect = contextlib.redirect_stdout(io.StringIO())
        redirect.__enter__()
        self.addCleanup(redirect.__exit__, None, None, None)

    def source(self, name, rows, columns=COLUMNS):
        path = self.root / (name + pp.RESULT_SUFFIX)
        with gzip.open(path, 'wt', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
        return path

    def status(self, names, success='TRUE'):
        with (self.root/'GWAMA_Run_Status.csv').open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['Output','Success'])
            writer.writerows((name, success) for name in names)

    def read(self, path):
        return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)

    def test_values_headers_stable_order_counts_and_override_scope(self):
        first = [row('0001','chrX','1'), row('rsTie1','chr2','020'), row('rsMT','MT','1')]
        second = [row('rsTie2','2.0','2e1','++??--'*30), row('rsEarly','1','5')]
        first[0]['Extra'] = 'a\tb "quoted"\nnew line'; second[1]['Extra'] = ''
        paths = [self.source('a', first), self.source('b', second)]
        combined, summary = pp.combine_results(paths, 33000.0, .9)
        self.assertEqual(combined['SNPID'].to_list(), ['rsEarly','rsTie1','rsTie2','0001','rsMT'])
        self.assertEqual(combined.columns, COLUMNS + ['count_question','count_plus','count_minus'])
        self.assertEqual(summary.columns, pp.GPCA_COLUMNS)
        originals = {r['SNPID']: r for r in first + second}
        for result in combined.to_dicts():
            original = originals[result['SNPID']]
            for key, value in original.items(): self.assertEqual(result[key], value)
            for label, char in [('question','?'),('plus','+'),('minus','-')]:
                self.assertEqual(result['count_'+label], original['Direction'].count(char))
        self.assertTrue(summary['N_eff'].eq('33000.0').all()); self.assertTrue(summary['INFO'].eq('0.9').all())
        self.assertEqual(summary['PVAL'].to_list(), combined['PVAL'].to_list())
        self.assertTrue(combined['N_eff'].eq('32817.123456789').all()); self.assertTrue(combined['INFO'].eq('NA').all())
        for backend in ('fallback', 'pigz'):
            if backend == 'pigz' and not shutil.which('pigz'): continue
            kwargs = {'return_value': None} if backend == 'fallback' else {'wraps': shutil.which}
            with patch.object(pp.shutil, 'which', **kwargs):
                output = self.root/backend; output.mkdir()
                audit = {}
                pp.save_outputs(combined, summary, output/'c.gz', output/'s.gz', output/'a.json', audit, n_cores=2)
                # Compare the entire decompressed TSV, including precise literals,
                # empty values, quoting and the original column order.
                for path, data in [(output/'c.gz',combined),(output/'s.gz',summary)]:
                    with gzip.open(path,'rt') as stream: self.assertEqual(stream.read(), pd.DataFrame(data.to_dict(as_series=False)).to_csv(sep='\t',index=False))
                saved = json.loads((output/'a.json').read_text())
                self.assertEqual(saved['compression']['workers'], 1 if backend == 'fallback' else 2)
                self.assertEqual(saved['compression']['level'], 1)
                self.assertEqual(set(saved['timings_seconds']), {'write_combined','write_summary'})

    def test_preserve_metadata_and_union_extra_columns(self):
        a = row('a'); b = row('b','2'); a['INFO']='0.876543210987654321'; b['Other']='extra'
        files = [self.source('a',[a]),self.source('b',[b],COLUMNS+['Other'])]
        combined, summary = pp.combine_results(files,None,None)
        self.assertEqual(summary['INFO'].to_list(),[a['INFO'],b['INFO']])
        self.assertEqual(summary['N_eff'].to_list(),[a['N_eff'],b['N_eff']])
        self.assertEqual(combined['Other'][1],'extra');self.assertIsNone(combined['Other'][0])

    def test_missing_metadata_requires_explicit_corresponding_override(self):
        for absent in ('INFO','N_eff'):
            path=self.source('a',[row()], [c for c in COLUMNS if c != absent])
            with self.assertRaisesRegex(ValueError,'missing columns'): pp.combine_results([path],None,None)
            combined,summary=pp.combine_results([path],33000 if absent=='N_eff' else None,.9 if absent=='INFO' else None)
            self.assertNotIn(absent,combined.columns);self.assertIn(absent,summary.columns)
        path=self.source('b',[row()],[c for c in COLUMNS if c not in ('N_eff','INFO')])
        combined,summary=pp.combine_results([path],33000,.9)
        self.assertEqual(summary.columns,pp.GPCA_COLUMNS)

    def test_invalid_rows_and_empty_inputs_fail(self):
        for column, value in [('SNPID',' '),('Direction',''),('Direction','+-0'),('Direction','+-\n'),
                              ('CHR','0'),('CHR','NA'),('BP','nan'),('BP','inf'),('BP','-1'),('BP','1.5'),('BP','')]:
            with self.subTest(column=column,value=value):
                value_row=row();value_row[column]=value
                with self.assertRaises(ValueError): pp.combine_results([self.source('bad',[value_row])],None,None)
        with self.assertRaisesRegex(ValueError,'empty'): pp.combine_results([self.source('empty',[])],None,None)
        with self.assertRaisesRegex(ValueError,'unique'): pp.combine_results([self.source('dup',[row(),row()])],None,None)
        with self.assertRaisesRegex(ValueError,'unique'): pp.combine_results([self.source('a',[row()]),self.source('b',[row()])],None,None)

    def test_numerical_failures_are_audited_before_overrides_and_export(self):
        for field in ('BETA', 'Z', 'SE', 'N_eff', 'PVAL'):
            values = ['NA', 'NaN', 'inf', '-inf', 'broken', '']
            if field in ('SE', 'N_eff'):
                values += ['0', '-1']
            if field == 'PVAL':
                values += ['-0.1', '1.1']
            for value in values:
                with self.subTest(field=field, value=value):
                    item = row('invalid'); item[field] = value
                    path = self.source('numerical', [item])
                    before = path.read_bytes()
                    # An explicit summary N override cannot conceal invalid raw N.
                    with self.assertRaisesRegex(ValueError, 'GWAMA output QC failed'):
                        pp.combine_results([path], 33000, .9)
                    issues = pd.read_csv(self.root/'numerical.GWAMA_Export_QC_Issues.csv')
                    self.assertEqual(issues.SNPID.tolist(), ['invalid'])
                    self.assertTrue(issues.Reason.str.contains(field).all())
                    self.assertEqual(path.read_bytes(), before)
        self.status(['numerical'])
        with self.assertRaisesRegex(ValueError, 'GWAMA output QC failed'):
            pp.process_gwama_results(self.root, self.root/'summary', name='failed', archive=True)
        self.assertFalse((self.root/'chromosome_wise').exists())
        self.assertFalse(list(self.root.glob('*combined*')))
        self.assertFalse((self.root/'summary').exists())
        attempts = pd.read_csv(self.root/'GWAMA_Export_Worker_Attempts.csv')
        self.assertTrue((~attempts.Success).all())

    def test_valid_probability_underflow_and_zero_effect_are_preserved(self):
        items = [row('underflow'), row('zero_effect')]
        items[0].update(PVAL='0', Z='1000', BETA='1')
        items[1].update(PVAL='1', Z='0', BETA='0')
        combined, _ = pp.combine_results([self.source('valid_numeric', items)], None, None)
        self.assertEqual(combined['PVAL'].to_list(), ['0', '1'])
        self.assertEqual(combined['Z'].to_list(), ['1000', '0'])
        summary = pd.read_csv(self.root/'valid_numeric.GWAMA_Export_QC_Summary.csv')
        self.assertEqual(summary.Status.tolist(), ['passed'])
        self.assertEqual(summary.Invalid_Rows.tolist(), [0])

    def test_chromosome_aliases_and_large_integer_sort_are_preserved(self):
        inputs=[row('a','X'),row('b','Y'),row('c','XY'),row('d','M'),row('e','chrMT'),
                row('large2','1','9007199254740993'),row('large1','1','9007199254740992')]
        combined,_=pp.combine_results([self.source('aliases',inputs)],None,None)
        self.assertEqual(combined['SNPID'].to_list(),['large1','large2','a','b','c','d','e'])

    def test_already_sorted_fast_path_preserves_ties(self):
        paths=[self.source('a',[row('a','1','10'),row('b','1','10')]),self.source('b',[row('c','2','1')])]
        with patch.object(pl.DataFrame,'sort',side_effect=AssertionError('redundant sort')):
            combined,_=pp.combine_results(paths,None,None)
        self.assertEqual(combined['SNPID'].to_list(),['a','b','c'])

    def test_current_run_selection_freshness_archive_and_audit(self):
        good=self.source('good',[row()]);stale=self.source('unlisted',[row()]);self.status(['good'])
        previous=pp.snapshot_outputs(self.root)
        with self.assertRaisesRegex(ValueError,'status was not updated'): pp.current_run_files(self.root,previous)
        previous.pop('GWAMA_Run_Status.csv')
        with self.assertRaisesRegex(ValueError,'result was not updated'): pp.current_run_files(self.root,previous)
        previous.pop(good.name)
        with patch.object(pp.shutil,'which',return_value=None):
            combined,summary=pp.process_gwama_results(self.root,self.root/'harmonisation_input',name='run',
                n_eff=33000,info_value=.9,archive=True,previous_files=previous,n_cores=22)
        self.assertTrue(stale.exists());self.assertFalse(good.exists());self.assertTrue((self.root/'chromosome_wise'/good.name).exists())
        self.assertEqual(self.read(combined).N_eff.iloc[0],'32817.123456789')
        self.assertEqual(self.read(summary).N_eff.iloc[0],'33000')
        audit=json.loads((self.root/'run_postprocess.json').read_text())
        self.assertEqual(audit['rows'],1);self.assertEqual(audit['sources'],[str(good.resolve())])
        self.assertEqual(audit['table_engine'],'polars');self.assertEqual(audit['read_workers'],1)
        self.assertEqual(audit['csv_threads_per_file'],1)
        self.assertIn('combine_and_prepare',audit['timings_seconds'])
        with self.assertRaises(FileExistsError): pp.process_gwama_results(self.root,self.root/'harmonisation_input',name='run')

    def test_status_failures_missing_files_and_archive_collision(self):
        self.source('a',[row()])
        with self.assertRaisesRegex(ValueError,'Missing run status'):pp.current_run_files(self.root,None)
        for names,success in [(['a'],'FALSE'),(['a','a'],'TRUE'),(['absent'],'TRUE'),(['../a'],'TRUE')]:
            self.status(names,success)
            with self.assertRaises(ValueError):pp.current_run_files(self.root,None)
        self.status(['a']);archive=self.root/'chromosome_wise';archive.mkdir();(archive/('a'+pp.RESULT_SUFFIX)).touch()
        with self.assertRaises(FileExistsError):pp.process_gwama_results(self.root,self.root/'summary',archive=True)
        self.assertFalse(list(self.root.glob('*combined*')))

    def test_write_failure_leaves_no_published_or_temporary_outputs(self):
        combined,summary=pp.combine_results([self.source('a',[row()])],None,None)
        targets=[self.root/'c.gz',self.root/'s.gz',self.root/'a.json']
        original=pp.write_compressed_table
        calls=[]
        def fail_second(data,destination,settings):
            calls.append(destination)
            if len(calls)==2:raise OSError('simulated full disk')
            return original(data,destination,settings)
        with patch.object(pp,'write_compressed_table',side_effect=fail_second),patch.object(pp.shutil,'which',return_value=None):
            with self.assertRaisesRegex(OSError,'full disk'):pp.save_outputs(combined,summary,*targets,{})
        self.assertTrue(all(not p.exists() for p in targets));self.assertFalse(list(self.root.glob('.gpca-*')))

    def test_publish_race_preserves_existing_destination_and_rolls_back(self):
        combined,summary=pp.combine_results([self.source('a',[row()])],None,None)
        c,s,a=[self.root/name for name in ('c.gz','s.gz','a.json')];s.write_text('other process')
        with patch.object(pp.shutil,'which',return_value=None),self.assertRaises(FileExistsError):pp.save_outputs(combined,summary,c,s,a,{})
        self.assertEqual(s.read_text(),'other process');self.assertFalse(c.exists());self.assertFalse(a.exists())
        self.assertFalse(list(self.root.glob('.gpca-*')))

    def test_compressor_failure_is_fatal_and_environment_cannot_change_format(self):
        combined,summary=pp.combine_results([self.source('a',[row()])],None,None)
        false=shutil.which('false')
        if false:
            with patch.object(pp.shutil,'which',return_value=false),self.assertRaises(OSError):
                pp.save_outputs(combined,summary,self.root/'c.gz',self.root/'s.gz',self.root/'a.json',{})
            self.assertFalse((self.root/'c.gz').exists());self.assertFalse(list(self.root.glob('.gpca-*')))
        if shutil.which('pigz'):
            with patch.dict(os.environ,{'PIGZ':'--zip','GZIP':'--best'}):
                pp.save_outputs(combined,summary,self.root/'c.gz',self.root/'s.gz',self.root/'a.json',{},n_cores=2)
            self.assertEqual((self.root/'c.gz').read_bytes()[:2],b'\x1f\x8b')
            self.assertEqual(self.read(self.root/'c.gz').SNPID.tolist(),['rs01'])

    def test_large_output_crosses_csv_chunk_boundary(self):
        data=pl.DataFrame({'SNPID':[f'rs{i}' for i in range(100003)],'PVAL':['1.23456789e-310']*100003})
        with patch.object(pp.shutil,'which',return_value=None):
            settings=pp.compression_settings(4,1);pp.write_compressed_table(data,self.root/'chunk.gz',settings)
        observed=self.read(self.root/'chunk.gz')
        pd.testing.assert_frame_equal(pd.DataFrame(data.to_dict(as_series=False)),observed)

    def test_invalid_override_and_compression_settings(self):
        for n,info in [(0,None),(-1,None),(float('inf'),None),(None,-.1),(None,1.1),(None,float('nan'))]:
            with self.assertRaises(ValueError):pp.validate_overrides(n,info)
        pp.validate_overrides(1,0);pp.validate_overrides(1,1)
        for cores,level in [(-1,1),(1.5,1),(True,1),(1,0),(1,10),(1,1.5),(1,True)]:
            with self.assertRaises(ValueError):pp.compression_settings(cores,level)
        with patch.object(pp.shutil,'which',return_value='/pigz'):
            self.assertEqual(pp.compression_settings(22,9)['workers'],22)
            self.assertLessEqual(pp.compression_settings(0,1)['workers'],22)

    def test_both_gpca_backends_forward_shared_export_settings(self):
        for script in ('gpsca_gwama_python_ldsc.r','gpsca_gwama_v2.r'):
            with patch.object(gpca.shutil,'which',return_value='/Rscript'),patch.object(gpca.subprocess,'run') as run,patch.object(pp,'process_gwama_results') as export:
                run.return_value.returncode=0
                args=['--input','traits','--outdir',str(self.root),'--gpca_input_folder','prepared','--n_cores','22','--gzip_level','3']
                self.assertEqual(gpca.main(args,r_script=script),0)
                self.assertEqual(export.call_args.kwargs['n_cores'],22);self.assertEqual(export.call_args.kwargs['gzip_level'],3)
                self.assertNotIn('--gzip_level',run.call_args.args[0])
        for value in ('0','10','fast'):
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):gpca.postprocess_parser().parse_args(['--gzip_level',value])

    def test_concurrent_read_budget_and_completion_order_do_not_change_ties(self):
        files=[self.source(name,[row(name,'1','123')]) for name in ('a','b','c','d')]
        original=pp.read_chromosome
        lock=threading.Lock();barrier=threading.Barrier(2);second_finished=threading.Event()
        state={'active':0,'peak':0,'finished':[]}
        def read(path,required):
            with lock:
                state['active']+=1;state['peak']=max(state['peak'],state['active'])
            try:
                if path in files[:2]:barrier.wait(timeout=10)
                if path==files[0]:self.assertTrue(second_finished.wait(timeout=10))
                result=original(path,required)
                with lock:state['finished'].append(path)
                if path==files[1]:second_finished.set()
                return result
            finally:
                with lock:state['active']-=1
        with patch.object(pp,'read_chromosome',side_effect=read):
            combined,_=pp.combine_results(files,None,None,n_cores=2)
        self.assertEqual(state['peak'],2);self.assertEqual(state['active'],0)
        self.assertLess(state['finished'].index(files[1]),state['finished'].index(files[0]))
        self.assertEqual(combined['SNPID'].to_list(),['a','b','c','d'])

    def test_empty_file_list_and_malformed_headers_fail(self):
        with self.assertRaisesRegex(ValueError,'No current-run'):pp.combine_results([],None,None)
        path=self.root/('invalid'+pp.RESULT_SUFFIX)
        for header in ('','SNPID\tSNPID\n','SNPID\t\n'):
            with gzip.open(path,'wt') as stream:stream.write(header)
            with self.assertRaisesRegex(ValueError,'headers'):pp.combine_results([path],None,None,n_cores=1)

    def test_corrupted_gzip_and_failed_worker_do_not_publish_outputs(self):
        path=self.source('bad',[row(f'rs{i}') for i in range(10000)])
        data=bytearray(path.read_bytes());data[-8]^=0xff;path.write_bytes(data)
        self.status(['bad'])
        with self.assertRaises((OSError,ValueError,EOFError)):
            pp.process_gwama_results(self.root,self.root/'summary',name='bad_crc',n_cores=2)
        self.assertFalse((self.root/'bad_crc_GWAMA_combined_results.txt.gz').exists())
        self.assertFalse(list(self.root.glob('.gpca-*')))

    def test_csv_parser_is_single_threaded_and_preserves_literal_input(self):
        path=self.source('file[1]',[row('0001')])
        original=pl.read_csv
        with patch.object(pl,'read_csv',wraps=original) as reader:
            combined,_=pp.combine_results([path],None,None,n_cores=4)
        kwargs=reader.call_args.kwargs
        self.assertEqual(kwargs['n_threads'],1);self.assertEqual(kwargs['infer_schema_length'],0)
        self.assertIs(kwargs['glob'],False);self.assertTrue(kwargs['missing_utf8_is_empty_string'])
        self.assertEqual(combined['SNPID'].to_list(),['0001'])
        self.assertEqual(combined['PVAL'][0],'1.234567890123456e-310')

    def test_polars_thread_pool_is_set_before_import_and_environment_override_is_kept(self):
        source_root=str(Path(pp.__file__).resolve().parents[1])
        code='''import os,sys
expected=int(sys.argv[1])
def check(event,args):
    if event=='import' and args[0]=='polars':
        assert os.environ.get('POLARS_MAX_THREADS')==str(expected)
sys.addaudithook(check)
from ldsc_gpca.gpca import main
assert main(['--postprocess_help',*sys.argv[2:]])==0
import polars as pl
assert pl.thread_pool_size()==expected,(pl.thread_pool_size(),expected)
'''
        for args,override,expected in [(['--n_cores','3'],None,3),(['--n_cores=4'],None,4),(['--n_cores','3'],'2',2),(['--n_cores','1'],None,1)]:
            env={k:v for k,v in os.environ.items() if k!='POLARS_MAX_THREADS'}
            env.update({'PYTHONPATH':source_root,'PYTHONDONTWRITEBYTECODE':'1'})
            if override is not None:env['POLARS_MAX_THREADS']=override
            result=subprocess.run([sys.executable,'-c',code,str(expected),*args],capture_output=True,text=True,env=env)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main()
