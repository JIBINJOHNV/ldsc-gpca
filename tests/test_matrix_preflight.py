"""Exercise the installed CLI paths and R implementation, not mocked matrices alone."""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')


@unittest.skipUnless(RSCRIPT, 'Rscript unavailable')
class MatrixPreflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.manifest = self.folder / 'traits.csv'
        self.manifest.write_text('traitname\nC\nA\nB\n')
        self.source = self.folder / 'ldsc.csv'
        self.traits = ['C', 'A', 'B']
        fields = ['p1','p2','rg','se','z','p','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se','extra']
        with self.source.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for i, a in enumerate(self.traits):
                for j, b in enumerate(self.traits):
                    rg = 1 if i == j else .2
                    intercept = 1 if i == j else (-.9 if {i,j} == {1,2} else .9)
                    writer.writerow(dict(p1=a,p2=b,rg=rg,se=.1,z=rg/.1,p=math.erfc(abs(rg/.1)/math.sqrt(2)),
                        h2_obs=.2,h2_obs_se=.02,h2_int=1,h2_int_se=.01,gcov_int=intercept,gcov_int_se=.01,extra='keep, exactly'))

    def command(self, out, *extra, backend='python'):
        route = ['gpca'] if backend == 'python' else ['genomicsem','gpca']
        result = subprocess.run([sys.executable,'-m','ldsc_gpca',*route,
            '--input',str(self.manifest),'--ldsc_results',str(self.source),
            '--outdir',str(out),*extra], capture_output=True,text=True,
            env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        return result

    def test_r_numerical_cases(self):
        result = subprocess.run([RSCRIPT,str(ROOT/'tests/test_matrix_preflight.R'),str(ROOT)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_pc1_only_retains_every_trait_and_skips_preparation_export(self):
        out = self.folder/'pc1'
        result = self.command(out,'--pc1_only')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        with (out/'GenomicPCA_Matrix_Validation.csv').open() as f: row = next(csv.DictReader(f))
        self.assertEqual(row['PC1_Status'],'PASS');self.assertEqual(row['CTI_Status'],'INDEFINITE')
        self.assertEqual(row['Traits'],'3')
        self.assertIn('not_run_pc1_only',(out/'GWAMA_Run_Status.csv').read_text())
        self.assertFalse((out/'gpca_inputs').exists());self.assertFalse((out/'harmonisation_input').exists())
        self.assertEqual((out/'Python_LDSC_Retained_Results.csv').read_text(),self.source.read_text())

    def test_strict_combined_validation_reports_both_and_fails(self):
        out = self.folder/'strict'
        result = self.command(out,'--validate_only')
        self.assertNotEqual(result.returncode,0)
        with (out/'GenomicPCA_Matrix_Validation.csv').open() as f: row = next(csv.DictReader(f))
        self.assertEqual(row['PC1_Status'],'PASS');self.assertEqual(row['CTI_Status'],'INDEFINITE')
        self.assertFalse((out/'GenomicPCA_Retained_Traits.csv').exists())

    def test_explicit_selection_preserves_rows_and_order_and_rechecks(self):
        out = self.folder/'selected'
        result = self.command(out,'--validate_only','--cti_action','explore_drop','--max_cti_drop_fraction','.34')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        with (out/'GenomicPCA_Retained_Traits.csv').open() as f: traits = [r['traitname'] for r in csv.DictReader(f)]
        self.assertEqual(len(traits),2)
        self.assertEqual(traits,[t for t in self.traits if t in traits])
        with (out/'GenomicPCA_Matrix_Validation_History.csv').open() as f: history=list(csv.DictReader(f))
        self.assertEqual(len(history),2);self.assertEqual(history[-1]['CTI_Status'],'PASS')
        original=self.source.read_text().splitlines(keepends=True)
        expected=[original[0]]+[line for line in original[1:] if all(t in traits for t in next(csv.reader([line]))[:2])]
        self.assertEqual((out/'Python_LDSC_Retained_Results.csv').read_text(),''.join(expected))
        # Exported retained source must independently pass strict validation.
        self.source=out/'Python_LDSC_Retained_Results.csv';self.manifest=out/'GenomicPCA_Retained_Traits.csv'
        rerun=self.command(self.folder/'rerun','--validate_only')
        self.assertEqual(rerun.returncode,0,rerun.stdout+rerun.stderr)

    def test_existing_assessment_cannot_be_overwritten(self):
        out=self.folder/'existing';out.mkdir()
        report=out/'GenomicPCA_Matrix_Validation.csv';report.write_text('previous result\n')
        result=self.command(out,'--pc1_only')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('fresh --outdir',result.stderr)
        self.assertEqual(report.read_text(),'previous result\n')

    def test_covariance_selection_retains_scale_and_source_provenance(self):
        out=self.folder/'covariance'
        result=self.command(out,'--validate_only','--pca_matrix','covariance',
            '--cti_action','explore_drop','--max_cti_drop_fraction','.34')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        with (out/'Python_LDSC_Input_Validation_Summary.csv').open() as f: rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),2)
        self.assertTrue(all(r['PCA_Matrix_Type']=='covariance' and r['Covariance_Matrix_Source'] for r in rows))

    def test_limit_exhaustion_publishes_no_usable_subset(self):
        out=self.folder/'limited'
        result=self.command(out,'--validate_only','--cti_action','explore_drop','--max_cti_drop_fraction','.1')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('exclusion limit',result.stderr)
        self.assertFalse((out/'GenomicPCA_Retained_Traits.csv').exists())
        self.assertFalse((out/'Python_LDSC_Retained_Results.csv').exists())

    def test_invalid_options_fail(self):
        for i,extra in enumerate([['--pc1_only','--validate_only'],['--pc1_only','--cti_action','explore_drop','--max_cti_drop_fraction','.5'],['--validate_only','--cti_action','explore_drop'],['--validate_only','--max_cti_drop_fraction','NaN']]):
            result=self.command(self.folder/str(i),*extra)
            self.assertNotEqual(result.returncode,0,result.stdout)

    def test_auto_preparation_is_blocked_by_failed_matrix_preflight(self):
        out=self.folder/'auto'
        result=self.command(out)
        self.assertNotEqual(result.returncode,0)
        self.assertTrue((out/'matrix_preflight/GenomicPCA_Matrix_Validation.csv').exists(),result.stdout+result.stderr)
        self.assertFalse((out/'gpca_inputs').exists())

    def test_genomicsem_uses_same_policy_and_subsets_sampling_covariance(self):
        # This is an RData fixture, not a native GenomicSEM LDSC fit.
        script=self.folder/'fixture.R';self.source=self.folder/'input.RData'
        script.write_text('''traits<-c("C","A","B")
S<-matrix(.04,3,3,dimnames=list(traits,traits));diag(S)<-.2
I<-matrix(c(1,.9,.9,.9,1,-.9,.9,-.9,1),3,3,dimnames=list(traits,traits))
LDSCoutput<-list(S=S,S_Stand=S/.2,I=I,V=diag(.0004,6),V_Stand=diag(.01,6),N=matrix(1000,1,6))
save(LDSCoutput,file=commandArgs(TRUE)[1])
''')
        subprocess.run([RSCRIPT,str(script),str(self.source)],check=True,capture_output=True)
        for name,extra in [('pc1',['--pc1_only']),('selected',['--validate_only','--cti_action','explore_drop','--max_cti_drop_fraction','.34'])]:
            out=self.folder/('gsem_'+name)
            result=self.command(out,*extra,backend='genomicsem')
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            retained=3 if name=='pc1' else 2
            check='load(commandArgs(TRUE)[1]);stopifnot(nrow(LDSCoutput$S)=='+str(retained)+', nrow(LDSCoutput$V)=='+str(retained*(retained+1)//2)+')'
            test=subprocess.run([RSCRIPT,'-e',check,str(out/'GenomicSEM_LDSC_Used.RData')],capture_output=True,text=True)
            self.assertEqual(test.returncode,0,test.stderr)

if __name__=='__main__':unittest.main()
