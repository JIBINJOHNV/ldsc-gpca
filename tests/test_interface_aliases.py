"""Single-name interface acceptance and rejection of removed spellings."""
import contextlib
import csv
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca.interfaces import read_manifest, REMOVED_MANIFEST_COLUMNS, REMOVED_OPTIONS
from ldsc_gpca import genomicsem_ldsc, gpca, ldsc_cli, pairwise, prepare

ROOT = Path(__file__).resolve().parents[1]


def write_csv(path, columns, rows):
    with path.open('w', newline='') as stream:
        writer=csv.writer(stream);writer.writerow(columns);writer.writerows(rows)


class CliNamesTests(unittest.TestCase):
    def test_python_ldsc_single_names(self):
        x=ldsc_cli.parser.parse_args(['--input=a','--outdir=b','--ld_ref=ld','--ld_weights=w',
            '--hm3=hm3','--munged_dir=m','--n_cores=22'])
        self.assertEqual((x.input,x.outdir,x.ld_ref,x.ld_weights,x.hm3,x.munged_dir,x.n_cores),('a','b','ld','w','hm3','m',22))
        for action in ldsc_cli.parser._actions:
            if action.dest!='help':self.assertLessEqual(len(action.option_strings),1)

    def test_removed_names_are_rejected_even_beside_valid_names(self):
        base=['--input','in','--outdir','out','--ld_ref','ld']
        for flag in sorted(REMOVED_OPTIONS | {'--info-min','-input_file','--ldsc-env','--n_core'}):
            with self.subTest(flag=flag),contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                ldsc_cli.parser.parse_args(base+[flag,'x'])

    def test_native_modes_single_names(self):
        parser=genomicsem_ldsc.build_parser();base=['--input','in','--outdir','out','--ld_ref','ld']
        for mode in (['--hm3','h'],['--munged_dir','m'],['--munged_input']):
            x=parser.parse_args(base+mode+['--ld_weights','w','--n_cores','4'])
            self.assertEqual((x.ld_ref,x.ld_weights,x.n_cores),('ld','w',4))
        for extra in (['--ld','ld'],['--munge-output','m'],['--munged_dir','m','--munged_input']):
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):parser.parse_args(base+extra)

    def test_wrapper_and_preparation_have_one_name_per_setting(self):
        parser=gpca.postprocess_parser(include_prepare=True)
        x,remaining=parser.parse_known_args(['--dataset_id','set','--gwama_output_n_eff','33000',
            '--gwama_output_info','.9','--hm3','hm3','--write_munge_inputs','--prepare_workers','2','--n_cores','22'])
        self.assertEqual((x.dataset_id,x.n_eff,x.info_value,x.hapmap_file,x.prepare_workers),('set',33000,.9,'hm3',2))
        self.assertEqual(remaining,['--n_cores','22'])
        for action in parser._actions:self.assertLessEqual(len(action.option_strings),1)
        with patch.object(prepare,'prepare_inputs') as run:
            self.assertEqual(prepare.main(['--input','a','--outdir','b','--n_cores','3']),0)
            self.assertEqual(run.call_args.kwargs['prepare_workers'],3)
        for flag in ('--prepare_workers','--prepare-workers','--cores','--hapmap_file'):
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                prepare.main(['--input','a','--outdir','b',flag,'3'])

    def test_conflicting_repeats_invalid_missing_and_abbreviated_values(self):
        base=['--input','a','--outdir','b','--ld_ref','ld']
        for more in (['--input','different'],['--n_cores','2','--n_cores','3'],['--hm3'],['--n_cores','bad'],['--ld_wei','w']):
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):ldsc_cli.parser.parse_args(base+more)
        self.assertEqual(ldsc_cli.parser.parse_args(base+['--n_cores','022','--n_cores','22']).n_cores,22)

    def test_old_gpca_options_fail_before_preparation_or_r(self):
        for script in ('gpsca_gwama_python_ldsc.r','gpsca_gwama_v2.r'):
            for flag in ('--cores','--python_ldsc','--ldsc_path','--input_file','--dataset-id','--tolerance'):
                with patch.object(prepare,'prepare_inputs') as prep,patch.object(gpca.subprocess,'run') as run, \
                     contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                    gpca.main(['--input','a','--outdir','b',flag,'x'],r_script=script)
                prep.assert_not_called();run.assert_not_called()

    def test_gpca_r_arguments_and_exports(self):
        from ldsc_gpca import postprocess
        with tempfile.TemporaryDirectory() as folder:
            for script in ('gpsca_gwama_python_ldsc.r','gpsca_gwama_v2.r'):
                args=['--input','in.csv','--outdir',folder,'--ldsc_results','results','--gpca_input_folder','prepared',
                      '--n_cores','22','--dataset_id','test','--gwama_output_n_eff','33000','--gwama_output_info','.9']
                with patch.object(gpca.shutil,'which',return_value='/fake/Rscript'),patch.object(gpca.subprocess,'run') as run, \
                     patch.object(postprocess,'process_gwama_results') as export:
                    run.return_value.returncode=0
                    self.assertEqual(gpca.main(args,r_script=script),0)
                    self.assertIn('--ldsc_results',run.call_args.args[0]);self.assertNotIn('--dataset_id',run.call_args.args[0])
                    self.assertEqual(export.call_args.kwargs['name'],'test');self.assertEqual(export.call_args.kwargs['n_eff'],33000)


class ManifestNamesTests(unittest.TestCase):
    def test_master_manifest_preserves_strings_rows_and_columns(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'in.csv';cols=['traitname','vcf_files','ref','sample_prevalence','population_prevalence','sumstats_file','munged_file']
            rows=[['002','vcf','yes','NA','NA','raw','munged'],['001','vcf2','yes','.2','.05','raw2','munged2']]
            write_csv(path,cols,rows);observed,parsed=read_manifest(path)
            self.assertEqual(observed,cols);self.assertEqual([list(row.values()) for row in parsed],rows)

    def test_removed_headers_rejected_including_when_canonical_is_present(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'in.csv'
            for old in REMOVED_MANIFEST_COLUMNS:
                write_csv(path,['traitname',old],[['A','A'],['B','B']])
                with self.subTest(old=old),self.assertRaisesRegex(ValueError,'Unsupported manifest headers'):read_manifest(path)

    def test_duplicate_empty_and_malformed_csv_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'in.csv'
            for text in ('','traitname,traitname\nA,A\n','traitname,vcf_files\nA\n','traitname\nA,extra\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):read_manifest(path)

    def test_prepare_paths_order_missing_names_and_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve();(root/'in.vcf').touch();path=root/'in.csv'
            write_csv(path,['traitname','vcf_files'],[['002','in.vcf'],['001','in.vcf']])
            self.assertEqual(prepare.manifest_inputs(path),[('002',root/'in.vcf'),('001',root/'in.vcf')])
            for row in (['A','missing.vcf'],['','in.vcf']):
                write_csv(path,['traitname','vcf_files'],[row])
                with self.assertRaises(ValueError):prepare.manifest_inputs(path)

    def test_native_all_modes_and_prevalence_rules(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve()
            for name in ('002','001'):(root/(name+'.sumstats')).write_text('SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100\t1\n')
            path=root/'in.csv'
            for mode in ('unmunged','paths','directory'):
                columns=['traitname','sample_prevalence','population_prevalence'];rows=[['002','NA','NA'],['001','.2','.05']]
                if mode!='directory':
                    columns+=['munged_file' if mode=='paths' else 'sumstats_file'];rows=[row+[row[0]+'.sumstats'] for row in rows]
                extra=['--munged_dir',str(root)] if mode=='directory' else ['--munged_input'] if mode=='paths' else []
                write_csv(path,columns,rows)
                opts=genomicsem_ldsc.build_parser().parse_args(['--input',str(path),'--outdir','out','--ld_ref','ld',*extra])
                result,backend=genomicsem_ldsc.resolve_manifest(opts)
                self.assertEqual([row['traitname'] for row in result],['002','001'])
                self.assertEqual(result[1]['sampleprevalence'],.2);self.assertEqual(result[1]['populationprevalence'],.05)
                self.assertEqual(backend,'munge' if mode=='unmunged' else 'existing')
            write_csv(path,['traitname','sample_prevalence','population_prevalence','munged_file'],[['A','.2','','002.sumstats'],['B','','','001.sumstats']])
            opts=genomicsem_ldsc.build_parser().parse_args(['--input',str(path),'--outdir','out','--ld_ref','ld','--munged_input'])
            with self.assertRaisesRegex(ValueError,'both prevalences'):genomicsem_ldsc.resolve_manifest(opts)

    def test_python_ldsc_preserves_backend_inputs_and_weights(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve();path=root/'in.csv';out=root/'out'
            write_csv(path,['traitname','population_prevalence','sample_prevalence','ref'],[['002','','','yes'],['001','','','yes']])
            with patch.object(ldsc_cli,'check_runtime'),patch.object(ldsc_cli,'is_valid_gz',return_value=True),patch.object(ldsc_cli,'check_saved_filters'), \
                 patch.object(ldsc_cli,'parallel_ldsc_analysis',return_value=[]) as run,patch.object(ldsc_cli,'compile_results'),contextlib.redirect_stdout(io.StringIO()):
                ldsc_cli.main(['--input',str(path),'--outdir',str(out),'--ld_ref',str(root/'ld'),'--ld_weights',str(root/'w'),
                               '--munged_dir',str(root/'munged'),'--ldsc_only'])
            frame=run.call_args.args[4];self.assertEqual(frame.gwas_name.tolist(),['002','001']);self.assertEqual(frame.sample_size_column.tolist(),['NEF','NEF'])
            self.assertTrue(frame.pop_prevalence.isna().all());self.assertEqual(run.call_args.kwargs['ld_weights_dir'],str(root/'w')+os.sep)
            self.assertEqual(json.loads((out/'LDSC_Runtime.json').read_text())['ld_weights'],str(root/'w')+os.sep)

    def test_native_runner_receives_canonical_options(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve();ld=root/'ld';weights=root/'weights';munged=root/'munged'
            for path in (ld,weights,munged):path.mkdir()
            for path in (ld/'1.l2.ldscore.gz',ld/'1.l2.M_5_50',weights/'1.l2.ldscore.gz'):path.write_text('fixture')
            for name in ('002','001'):(munged/(name+'.sumstats')).write_text('SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100\t1\n')
            manifest=root/'in.csv';out=root/'out'
            write_csv(manifest,['traitname','sample_prevalence','population_prevalence'],[['002','NA','NA'],['001','NA','NA']])
            with patch.object(genomicsem_ldsc.shutil,'which',return_value='/fake/Rscript'),patch.object(genomicsem_ldsc.subprocess,'run') as run:
                run.return_value.returncode=0
                self.assertEqual(genomicsem_ldsc.main(['--input',str(manifest),'--outdir',str(out),'--ld_ref',str(ld),
                    '--ld_weights',str(weights),'--munged_dir',str(munged),'--n_cores','4','--chromosomes','1']),0)
                args=run.call_args.args[0]
                self.assertEqual(args[4:9],[str(ld),str(weights),'','existing','4'])
            with (out/'Resolved_Manifest.csv').open() as stream:rows=list(csv.DictReader(stream))
            self.assertEqual([row['traitname'] for row in rows],['002','001'])
            self.assertTrue(all(row['sampleprevalence']=='NA' for row in rows))

    def test_upstream_weights_flags_unchanged(self):
        frame=pd.DataFrame({'gwas_name':['A'],'ref':['yes'],'sample_prevalence':[float('nan')],'pop_prevalence':[float('nan')]})
        with tempfile.TemporaryDirectory() as folder:
            for weights in (None,'/separate weights/'):
                with patch.object(pairwise,'run_command') as run:
                    pairwise.parallel_ldsc_analysis(1,1,folder,'/ref/',frame,folder,ld_weights_dir=weights)
                    tokens=shlex.split(run.call_args.args[0]);self.assertEqual(tokens[tokens.index('--ref-ld-chr')+1],'/ref/')
                    self.assertEqual(tokens[tokens.index('--w-ld-chr')+1],weights or '/ref/')


class RInterfaceTests(unittest.TestCase):
    def test_real_r_parsers_and_manifest_names(self):
        rscript=os.environ.get('LDSC_GPCA_TEST_RSCRIPT')
        if not rscript:self.skipTest('Set LDSC_GPCA_TEST_RSCRIPT to R with argparse/data.table')
        result=subprocess.run([rscript,str(ROOT/'tests/test_interface_aliases.R'),str(ROOT)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('R single-name checks passed',result.stdout)


if __name__=='__main__':unittest.main()
