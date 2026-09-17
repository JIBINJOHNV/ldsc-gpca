import os
import io
import shlex
import shutil
import sys
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pandas as pd
from ldsc_gpca import ldsc_runtime, extraction, pairwise

class RuntimeTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('bcftools'), 'local bcftools not installed')
    def test_real_bcftools_extraction(self):
        header = '##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=6>\n'
        for name in ['AF','EUR']:
            header += f'##INFO=<ID={name},Number=1,Type=Float,Description="test">\n'
        for name in ['SI','AF','EZ','LP','NEF','NC','NCO']:
            header += f'##FORMAT=<ID={name},Number=1,Type=Float,Description="test">\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait\n'
        rows = ''.join(f'{chrom}\t{pos}\t{name}\tA\tG\t.\tPASS\tAF=0.2;EUR=0.2\tSI:AF:EZ:LP:NEF:NC:NCO\t{si}:0.2:2:8:360:100:900\n'
                       for chrom,pos,name,si in [('1',100,'rs1',.9),('1',200,'rs2',.1),('6',30000000,'rs3',.9)])
        with tempfile.TemporaryDirectory(prefix='real bcftools ') as root:
            vcf=Path(root)/'input file.vcf'; vcf.write_text(header+rows)
            extraction.munge_input_worker('q',root,root,str(vcf))
            result=pd.read_csv(Path(root)/'munge_input/q_mungeinput.tsv',sep=r'\s+')
            self.assertEqual(result.ID.tolist(),['rs1','rs3'])
            self.assertTrue((result.P==1e-8).all())
            filters={**extraction.DEFAULT_FILTERS,'exclude_mhc':True}
            prev=extraction.munge_input_worker_with_prevalence('b',root,root,str(vcf),None,filters)
            binary=pd.read_csv(Path(root)/'munge_input/b_mungeinput.tsv',sep='\t')
            self.assertEqual(binary.ID.tolist(),['rs1'])
            self.assertEqual(prev,.1)
            self.assertEqual(binary.N_TOTAL.tolist(),[1000])
            with self.assertRaises(subprocess.CalledProcessError):
                extraction.munge_input_worker('bad',root,root,str(Path(root)/'missing.vcf'))

    def test_default_command(self):
        self.assertEqual(ldsc_runtime.ldsc_command('ldsc.py'),
            ['conda', 'run', '--no-capture-output', '--name', 'ldsc-cbiit', 'ldsc.py'])

    def test_prefix_and_spaces(self):
        cmd = ldsc_runtime.ldsc_command('munge_sumstats.py', conda='/tools with spaces/conda', prefix='/child env')
        self.assertIn('--prefix', cmd)
        self.assertNotIn('--name', cmd)
        self.assertEqual(shlex.split(shlex.join(cmd)), cmd)

    def test_raw_script_passthrough_and_exit_status(self):
        arguments = ['--rg', '/a path/a.sumstats.gz,/b.sumstats.gz', '--chisq-max', '80']
        with patch.dict(os.environ, {
                'CONDA_EXE': '/tools with spaces/conda',
                'LDSC_GPCA_LDSC_PREFIX': '/child env'}, clear=False), \
                patch.object(ldsc_runtime.subprocess, 'run',
                             return_value=SimpleNamespace(returncode=7)) as run:
            self.assertEqual(ldsc_runtime.run_ldsc_script('ldsc.py', arguments), 7)
        self.assertEqual(run.call_args.args[0], [
            '/tools with spaces/conda', 'run', '--no-capture-output',
            '--prefix', '/child env', 'ldsc.py', *arguments])

        with self.assertRaisesRegex(ValueError, 'Unsupported raw LDSC script'):
            ldsc_runtime.run_ldsc_script('other.py', [])
        with patch.object(ldsc_runtime.subprocess, 'run', side_effect=FileNotFoundError('missing')), \
                patch('sys.stderr', new_callable=io.StringIO) as error:
            self.assertEqual(ldsc_runtime.run_ldsc_script('munge_sumstats.py', []), 127)
            self.assertIn('Unable to launch', error.getvalue())

    def test_missing_tools(self):
        with patch.object(ldsc_runtime.shutil, 'which', return_value=None), self.assertRaisesRegex(ValueError, 'Conda executable'):
            ldsc_runtime.check_runtime()
        with patch.object(ldsc_runtime.shutil, 'which', side_effect=lambda x: None if x=='absent' else x), self.assertRaisesRegex(ValueError, 'extraction executable'):
            ldsc_runtime.check_runtime(bcftools='absent')

    def test_preflight_success_and_failure(self):
        with patch.object(ldsc_runtime.shutil, 'which', return_value='/tool'), patch.object(ldsc_runtime.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            ldsc_runtime.check_runtime(prefix='/child env')
            self.assertEqual(run.call_count, 2)
            self.assertEqual(run.call_args.args[0][-2:], ['munge_sumstats.py', '--help'])
        with patch.object(ldsc_runtime.shutil, 'which', return_value='/tool'), patch.object(ldsc_runtime.subprocess, 'run', return_value=SimpleNamespace(returncode=1, stderr='broken import')), self.assertRaisesRegex(ValueError, 'broken import'):
            ldsc_runtime.check_runtime()

    def test_extraction_executable(self):
        commands = extraction.filter_commands('/input file.vcf', extraction.DEFAULT_FILTERS, '%ID', '/tools/bcftools')
        self.assertTrue(all(c[0]=='/tools/bcftools' for c in commands))
        self.assertNotIn('docker', shlex.join(sum(commands, [])))

    def test_custom_runtime_reaches_pairwise(self):
        frame=pd.DataFrame({'gwas_name':['q'], 'ref':['yes'], 'pop_prevalence':[float('nan')], 'sample_prevalence':[float('nan')]})
        with tempfile.TemporaryDirectory(prefix='ldsc runtime ') as folder, patch.object(pairwise, 'run_command') as run:
            pairwise.parallel_ldsc_analysis(1,1,folder,folder,frame,folder,runtime={'prefix':'/child env'})
            args=shlex.split(run.call_args.args[0])
            self.assertEqual(args[:6], ['conda','run','--no-capture-output','--prefix','/child env','python'])
            self.assertTrue(args[6].endswith('/ldsc_export.py'))
            self.assertNotIn('--pop-prev',args)

    def test_setup_help_invalid_and_existing(self):
        script=Path(__file__).resolve().parents[1]/'scripts/setup_environments.sh'
        self.assertEqual(subprocess.run(['bash','-n',str(script)]).returncode,0)
        self.assertEqual(subprocess.run(['bash',str(script),'--help'],capture_output=True).returncode,0)
        self.assertEqual(subprocess.run(['bash',str(script),'a','b'],capture_output=True).returncode,2)
        with tempfile.TemporaryDirectory() as root:
            (Path(root)/'main').mkdir()
            result=subprocess.run(['bash',str(script),root],env={**os.environ,'CONDA_EXE':'true'},capture_output=True,text=True)
            self.assertEqual(result.returncode,2)
            self.assertIn('already exists',result.stderr)

    def test_setup_routes_and_stops_on_failure(self):
        script=Path(__file__).resolve().parents[1]/'scripts/setup_environments.sh'
        with tempfile.TemporaryDirectory(prefix='setup routing ') as root:
            fake=Path(root)/'conda'; log=Path(root)/'commands.jsonl'
            hook=Path(root)/'etc/profile.d/conda.sh'; hook.parent.mkdir(parents=True)
            hook.write_text('conda() { export CONDA_PREFIX="$2"; }\n')
            fake.write_text('#!'+sys.executable+'\nimport sys,os,json\n'
                            'with open(os.environ["COMMAND_LOG"],"a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\n'
                            'if sys.argv[1:]==["info","--base"]: print(os.path.dirname(os.path.abspath(sys.argv[0])))\n'
                            'sys.exit(int(os.environ.get("COMMAND_FAIL","0")))\n')
            fake.chmod(0o755)
            env={**os.environ,'PATH':'/usr/bin:/bin','CONDA_EXE':str(fake),'COMMAND_LOG':str(log)}
            result=subprocess.run(['bash',str(script),str(Path(root)/'environments')],env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            calls=[json.loads(line) for line in log.read_text().splitlines()]
            self.assertEqual(sum(c[:2]==['env','create'] for c in calls),2)
            config=[c for c in calls if c[:4]==['env','config','vars','set']][0]
            self.assertTrue(config[-1].endswith('/environments/ldsc'))
            self.assertTrue((Path(root)/'environments/activate.sh').exists())
            log.unlink()
            result=subprocess.run(['bash',str(script),str(Path(root)/'failure')],env={**env,'COMMAND_FAIL':'7'},capture_output=True)
            self.assertEqual(result.returncode,7)
            self.assertEqual(len(log.read_text().splitlines()),1)

if __name__=='__main__':
    unittest.main()
