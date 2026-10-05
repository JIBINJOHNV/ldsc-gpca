"""Preparation modes with real bcftools and optional real Python/native munging."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

import polars as pl
from ldsc_gpca import prepare, prepare_ldsc, ldsc_munge, ldsc_cli, genomicsem_ldsc, results
import test_genomicsem_vcf as fixtures
from ldsc_gpca.restart import file_digest


@unittest.skipUnless(shutil.which('bcftools'), 'bcftools not installed')
class PreparationModesTests(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.NativeVCFTests()
        self.case.setUp()
        self.addCleanup(self.case.doCleanups)
        self.root, self.manifest, self.out, self.hm3 = self.case.root, self.case.manifest, self.case.out, self.case.hm3
        self.hm3.write_text('SNP\tA1\tA2\nrs1\tG\tA\nrs2\tG\tA\n')
        self.capture = io.StringIO()

    def args(self, *extra, mode='ldsc', raw=True):
        return ['--input',str(self.manifest),'--outdir',str(self.out),'--mode',mode,
                '--hm3',str(self.hm3),*(['--raw_only'] if raw else []),*extra]

    def run_cli(self, args):
        with contextlib.redirect_stdout(self.capture), contextlib.redirect_stderr(self.capture):
            return prepare.main(args)

    def raw(self, **kwargs):
        self.assertEqual(self.run_cli(self.args(**kwargs)),0,self.capture.getvalue())

    def test_shared_raw_columns_order_qc_and_no_gpca(self):
        before=hashlib.sha256(self.case.vcf.read_bytes()).hexdigest()
        self.raw()
        frame=pl.read_csv(self.out/'munge_inputs/002_munge_inputs.tsv',separator='\t')
        self.assertEqual(frame.columns,['SNP','CHR','BP','A1','A2','EAF','BETA','SE','P','N','INFO'])
        self.assertEqual(frame['SNP'].to_list(),['rs1','rs2'])
        self.assertEqual(frame['N'].to_list(),[20000,22000])
        self.assertAlmostEqual(frame['EAF'][0],.8,places=6)
        self.assertAlmostEqual(frame['INFO'][1],.8,places=6)
        self.assertEqual(frame['P'][1],1e-300)
        self.assertFalse((self.out/'gpca_inputs').exists())
        self.assertFalse((self.out/'munged').exists())
        self.assertFalse((self.out/'GPCA_Input_QC_Summary.csv').exists())
        manifest=pl.read_csv(self.out/'Prepared_LDSC_Manifest.csv',schema_overrides={'traitname':pl.String})
        self.assertEqual(manifest['traitname'].to_list(),['002','001'])
        self.assertEqual(manifest['ref'].to_list(),['yes','yes'])
        self.assertEqual(before,hashlib.sha256(self.case.vcf.read_bytes()).hexdigest())

    def test_both_preserves_gpca_rows_and_nef_despite_ldsc_qc_and_n_override(self):
        baseline=self.root/'baseline'
        with contextlib.redirect_stdout(self.capture):
            prepare.prepare_inputs(self.manifest,baseline,splitby_chr='nosplit')
        self.case.write_manifest(n=35000)
        self.assertEqual(self.run_cli(self.args('--splitby_chr','nosplit',mode='both')),0,self.capture.getvalue())
        for name in ('002','001'):
            filename=f'gpca_inputs/{name}_GenomicPCA_inputs.tsv'
            self.assertEqual((baseline/filename).read_bytes(),(self.out/filename).read_bytes())
            raw=pl.read_csv(self.out/f'munge_inputs/{name}_munge_inputs.tsv',separator='\t')
            self.assertEqual(raw['N'].to_list(),[35000]*2)
        self.assertTrue((self.out/'GPCA_Input_QC_Summary.csv').exists())
        self.assertTrue((self.out/'LDSC_Input_QC_Summary.csv').exists())

    def test_gpca_default_and_legacy_export_work_without_si(self):
        self.case.write_vcf(missing='SI')
        base=['--input',str(self.manifest),'--outdir',str(self.out),'--splitby_chr','nosplit']
        self.assertEqual(self.run_cli(base),0,self.capture.getvalue())
        self.out=self.root/'legacy'
        base[3]=str(self.out)
        self.assertEqual(self.run_cli(base+['--write_munge_inputs','--hm3',str(self.hm3)]),1) # duplicate VCF IDs are diagnosed in legacy mode
        # The existing legacy ID convention also supports coordinate IDs.
        self.out=self.root/'legacy2'; base[3]=str(self.out)
        ids=['1_1_A_G','1_2_A_G']
        self.hm3.write_text('SNP\n'+'\n'.join(ids)+'\n')
        self.assertEqual(self.run_cli(base+['--write_munge_inputs','--hm3',str(self.hm3),'--munge_id_source','chr_pos_ref_alt']),0,self.capture.getvalue())
        self.assertEqual((self.out/'munge_inputs/002_munge_inputs.txt').read_text().splitlines()[0],
                         'SNP CHR POS A1 A2 eaf_A1 beta se N p')

    def test_n_override_without_nef_and_one_trait(self):
        self.case.write_vcf(missing='NEF')
        self.case.write_manifest(binary=True,n=40000)
        self.manifest.write_text('\n'.join(self.manifest.read_text().splitlines()[:2])+'\n')
        self.raw()
        self.assertEqual(pl.read_csv(self.out/'munge_inputs/002_munge_inputs.tsv',separator='\t')['N'].to_list(),[40000]*2)

    def test_invalid_modes_and_options_fail_before_outputs(self):
        cases=[self.args('--munge_backend','python'), self.args(mode='gpca'),
               self.args('--write_munge_inputs'), self.args('--info_filter','.8'),
               self.args('--maf_filter','nan',raw=False), self.args('--info_filter','2',raw=False),
               self.args('--rscript','Rscript',raw=False)]
        for args in cases:
            with self.subTest(args=args):
                self.assertEqual(self.run_cli(args),1,self.capture.getvalue())
                self.assertFalse(self.out.exists())

    def test_bad_metadata_rejected_before_outputs(self):
        for n in (None,0,-1,'inf'):
            self.case.write_manifest(binary=True,n=n)
            self.assertEqual(self.run_cli(self.args()),1,self.capture.getvalue())
            self.assertFalse(self.out.exists())
        self.case.write_manifest(binary=True,n=50000)
        self.manifest.write_text(self.manifest.read_text().replace(',.1,',',,'))
        self.assertEqual(self.run_cli(self.args()),1)

    def test_missing_si_retries_and_publishes_no_tables(self):
        self.case.write_vcf(missing='SI')
        self.assertEqual(self.run_cli(self.args()),1)
        attempts=pl.read_csv(self.out/'Preparation_Worker_Attempts.csv')
        self.assertEqual(attempts.height,4)
        self.assertFalse((self.out/'munge_inputs').exists())

    def test_both_flushes_independent_original_reports_when_ldsc_fails(self):
        self.hm3.write_text('SNP\nno_matching_snp\n')
        self.assertEqual(self.run_cli(self.args('--splitby_chr', 'nosplit', mode='both')), 1)
        self.assertIn('no usable LDSC records', self.capture.getvalue())
        gpca = (self.out/'GPCA_Input_QC_Issues.csv').read_text()
        ldsc = (self.out/'LDSC_Input_QC_Issues.csv').read_text()
        self.assertIn('P below', gpca)
        self.assertIn('SNP not in HapMap reference', ldsc)
        self.assertIn('#CHROM', gpca)
        self.assertIn('#CHROM', ldsc)
        self.assertFalse((self.out/'gpca_inputs').exists())
        self.assertFalse((self.out/'munge_inputs').exists())
        self.assertEqual(pl.read_csv(self.out/'Preparation_Worker_Attempts.csv').height, 4)

    def test_both_flushes_gpca_report_when_gpca_write_fails(self):
        with patch.object(prepare, 'write_gpca', side_effect=OSError('disk failure')):
            self.assertEqual(self.run_cli(self.args('--splitby_chr', 'nosplit', mode='both')), 1)
        self.assertIn('disk failure', self.capture.getvalue())
        self.assertIn('P below', (self.out/'GPCA_Input_QC_Issues.csv').read_text())
        self.assertFalse((self.out/'gpca_inputs').exists())

    def test_existing_outputs_are_preserved(self):
        self.raw()
        before=(self.out/'munge_inputs/002_munge_inputs.tsv').read_bytes()
        self.assertEqual(self.run_cli(self.args()),1)
        self.assertEqual(before,(self.out/'munge_inputs/002_munge_inputs.tsv').read_bytes())

    def test_hapmap_removals_are_audited(self):
        self.hm3.write_text('SNP\trs_note\nrs1\tx\n')
        self.raw()
        qc=pl.read_csv(self.out/'LDSC_Input_QC_Summary.csv')
        self.assertEqual(qc['retained_rows'].to_list(),[1,1])
        self.assertEqual(qc['hapmap_removed_rows'].to_list(),[1,1])
        self.assertEqual(qc['p_adjusted_rows'].to_list(),[0,0])
        self.assertIn('SNP not in HapMap reference',(self.out/'LDSC_Input_QC_Issues.csv').read_text())

    def test_munging_failure_retains_raw_and_two_attempt_audit(self):
        self.raw()
        opts=prepare.build_parser().parse_args(self.args(raw=False))
        with contextlib.redirect_stdout(self.capture), self.assertRaisesRegex(ValueError,'Munging failed after 2 attempts'):
            prepare_ldsc.munge_prepared(opts,[sys.executable,'-c','raise SystemExit(7)'])
        self.assertTrue((self.out/'munge_inputs').is_dir())
        self.assertFalse((self.out/'munged').exists())
        self.assertEqual(pl.read_csv(self.out/'Preparation_Munging_Worker_Attempts.csv').height,4)
        self.assertEqual(json.loads((self.out/'Preparation_Settings.json').read_text())['munging_status'],'failed')

    def test_existing_filter_provenance_checks_remain_active(self):
        with self.assertRaisesRegex(ValueError,'filter settings differ'):
            results.check_saved_filters({'filters': {'info_min': .7}}, {'info_min': .9}, 'A')
        results.check_saved_filters({'filters': {'info_min': .7}}, {'info_min': .7}, 'A')

    @unittest.skipUnless(shutil.which('munge_sumstats.py') or (Path(sys.executable).parent/'munge_sumstats.py').exists(), 'Python LDSC not installed')
    def test_real_python_cli_matches_direct_munging_and_sidecars(self):
        # A one-SNP fixture must still pass the original BETA median sanity check.
        self.case.vcf.write_text(self.case.vcf.read_text().replace('-0.123456789', '-0.08'))
        self.case.write_manifest(binary=True,n=40000)
        command=[sys.executable,str(Path(ldsc_munge.__file__).resolve())]
        env={'PATH':str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH']}
        with patch.dict(os.environ,env), patch('ldsc_gpca.ldsc_runtime.ldsc_munging_command',return_value=command):
            self.assertEqual(self.run_cli(self.args(raw=False)),0,self.capture.getvalue())
            log=json.loads((self.out/'munge_logs/002.command.json').read_text())
            log[log.index('--out')+1]=str(self.root/'direct')
            subprocess.run(log,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
        got=pl.read_csv(self.out/'munged/002.sumstats.gz',separator='\t')
        expected=pl.read_csv(self.root/'direct.sumstats.gz',separator='\t')
        self.assertTrue(got.equals(expected))
        self.assertEqual(got.drop_nulls()['SNP'].to_list(),['rs1'])
        sidecar=json.loads((self.out/'munged/002.prevalence.json').read_text())
        self.assertEqual(sidecar['sample_size_column'],'N_TOTAL')
        self.assertEqual(sidecar['sample_prevalence'],.5)
        self.assertEqual(sidecar['sha256'],file_digest(self.out/'munged/002.sumstats.gz'))
        with contextlib.redirect_stdout(self.capture): results.check_saved_filters(sidecar,{},'002')
        self.assertIn('VCF extraction filters are not applied',self.capture.getvalue())
        manifest=pl.read_csv(self.out/'Prepared_LDSC_Manifest.csv')
        self.assertTrue(all(Path(path).exists() for path in manifest['munged_file']))
        self.check_reuse()
        settings_path = self.out/'Preparation_Settings.json'
        original_settings = settings_path.read_text()
        settings = json.loads(original_settings)
        settings['munging_status'] = 'failed'
        settings_path.write_text(json.dumps(settings))
        with self.assertRaisesRegex(ValueError, 'not a completed'):
            self.check_reuse()
        settings_path.write_text(original_settings)
        sidecar['sha256']='changed'
        (self.out/'munged/002.prevalence.json').write_text(json.dumps(sidecar))
        with self.assertRaisesRegex(ValueError,'munged file or N convention changed'):
            self.check_reuse()

    @unittest.skipUnless(os.environ.get('LDSC_GPCA_TEST_RSCRIPT'), 'set LDSC_GPCA_TEST_RSCRIPT for installed native munging')
    def test_native_cli_matches_direct_munging_and_allele_alignment(self):
        rscript = os.environ['LDSC_GPCA_TEST_RSCRIPT']
        self.hm3.write_text('SNP\tA1\tA2\nrs1\tA\tG\nrs2\tG\tA\n')
        self.assertEqual(self.run_cli(self.args('--munge_backend','genomicsem','--rscript',rscript,
                                               '--splitby_chr','nosplit',mode='both',raw=False)),0,self.capture.getvalue())
        script = self.root/'direct.R'
        script.write_text('a <- commandArgs(TRUE)\nGenomicSEM::munge(files=a[1], hm3=a[2], trait.names=a[3], N=NA, info.filter=.9, maf.filter=.01, parallel=FALSE, cores=1L, overwrite=TRUE)\n')
        subprocess.run([rscript,str(script),str(self.out/'munge_inputs/002_munge_inputs.tsv'),
                        str(self.hm3),'direct'],cwd=self.root,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
        got=pl.read_csv(self.out/'munged/002.sumstats.gz',separator='\t')
        self.assertTrue(got.equals(pl.read_csv(self.root/'direct.sumstats.gz',separator='\t')))
        self.assertEqual(got['SNP'].to_list(),['rs1'])
        self.assertEqual(got['A1'].to_list(),['A'])
        self.assertGreater(got['Z'][0],0)  # Native output aligns the negative ALT effect to reference A1.
        self.assertTrue((self.out/'gpca_inputs/002_GenomicPCA_inputs.tsv').exists())
        self.assertEqual(json.loads((self.out/'munged/002.prevalence.json').read_text())['munge_backend'],'genomicsem')
        self.check_reuse()

    def check_reuse(self):
        manifest=str(self.out/'Prepared_LDSC_Manifest.csv')
        folder=str(self.out/'munged')
        args=['--input',manifest,'--outdir',str(self.root/'reuse'),'--ld_ref',str(self.case.ld),'--munged_dir',folder]
        # Exercise real file/provenance validation, stopping at the regression boundary.
        with patch.object(ldsc_cli,'check_runtime'), patch.object(ldsc_cli,'parallel_ldsc_analysis',return_value=[]) as run, \
                patch.object(ldsc_cli,'compile_results'), contextlib.redirect_stdout(self.capture):
            ldsc_cli.main([*args,'--ldsc_only'])
            self.assertEqual(run.call_args.args[4].gwas_name.tolist(),['002','001'])
        opts=genomicsem_ldsc.build_parser().parse_args(args)
        if json.loads((self.out/'munged/002.prevalence.json').read_text())['munge_backend'] == 'python':
            with self.assertRaisesRegex(ValueError, 'backend must be genomicsem'):
                genomicsem_ldsc.resolve_manifest(opts)
            return
        rows, mode=genomicsem_ldsc.resolve_manifest(opts)
        self.assertEqual(mode,'existing')
        self.assertEqual([r['traitname'] for r in rows],['002','001'])

if __name__=='__main__': unittest.main()
