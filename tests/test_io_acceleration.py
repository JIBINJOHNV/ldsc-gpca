"""Fast I/O preserves parsing, scientific literals and transactional failures."""
import contextlib
import gzip
import io
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import compression, munging, results
from test_results_csv import record
from ldsc_gpca.ldsc_export import write_results_csv


class AcceleratedIOTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_ld_tsv_identifiers_and_missing_tokens_match_pandas_without_pandas_read(self):
        path = self.root/'LD[1].gz'
        tokens = ['001', 'rs001', 'NA', 'N/A', 'NaN', 'nan', 'NULL', 'null', 'None',
                  '<NA>', '#NA', '#N/A', '-NaN', '-nan', '1.#IND', '-1.#IND',
                  '1.#QNAN', '-1.#QNAN', 'rs001']
        text = 'CHR\tSNP\tBP\tL2\n' + ''.join(f'1\t{snp}\t{i}\t1.5\n' for i,snp in enumerate(tokens))
        path.write_bytes(gzip.compress(text.encode()))
        expected = pd.read_csv(path,sep=r'\s+',usecols=['SNP'],dtype=str).SNP.dropna().tolist()
        with patch.object(munging.pd,'read_csv',side_effect=AssertionError('slow TSV reader')):
            self.assertEqual(munging._read_ld_snps(path),expected)

    def test_ld_legacy_whitespace_quotes_and_empty_fields_match_pandas(self):
        path = self.root/'ld.gz'
        for text in [
            'CHR SNP BP L2\n1 001 2 3\n1 NA 4 5\n',
            'CHR\tSNP\tBP\tL2\n1 rsMixed\t2\t3\n1\trsPlain\t3\t4\n',
            'CHR\tSNP\tBP\tL2\n1\t\trsShifted\t2\t3\n',
            'CHR\tSNP\tBP\tL2\r\n1\trsCRLF\t2\t3\r\n',
            'CHR\tSNP\tBP\tL2\n1\t"rs quoted"\t2\t3\n',
            'CHR\tSNP\tBP\tL2\n1\t#N/A N/A\t2\t3\n',
        ]:
            with self.subTest(text=text):
                path.write_bytes(gzip.compress(text.encode()))
                expected=pd.read_csv(path,sep=r'\s+',usecols=['SNP'],dtype=str).SNP.dropna().tolist()
                self.assertEqual(list(munging._read_ld_snps(path)),expected)

    def test_csv_special_names_bom_and_precise_values_survive(self):
        path=self.root/'results[1].csv'
        names=['001','NA','None','trait,with comma','trait "quote"','trait\nnewline']
        rows=[dict(record(),p1=name,p2=name) for name in names]
        write_results_csv(pd.DataFrame(rows),path)
        path.write_bytes(b'\xef\xbb\xbf'+path.read_bytes())
        with patch.object(results.pd,'read_csv',side_effect=AssertionError('pandas result reader')):
            actual=results._read_numerical_csv(path)
        self.assertEqual(actual.p1.tolist(),names)
        self.assertEqual(actual.p2.tolist(),names)
        for col,value in record().items():
            if col not in ('p1','p2'): self.assertTrue(actual[col].eq(value).all(),col)

    def test_csv_bad_quotes_duplicate_header_and_ragged_extra_values_fail(self):
        path=self.root/'results.csv'
        frame=pd.DataFrame([record()]); text=frame.to_csv(index=False)
        for bad in [text+'"unfinished', text.replace('p1,p2,','p1,p1,'),
                    text.rstrip('\n')+',extra\n']:
            with self.subTest(bad=bad):
                path.write_text(bad)
                with self.assertRaisesRegex(RuntimeError,'Cannot read numerical LDSC CSV'):
                    results._read_numerical_csv(path)

    def test_filter_backends_preserve_missing_rows_boundaries_and_source(self):
        source=self.root/'input.gz'
        text='SNP\tA1\tA2\tN\tZ\r\nrsBoundary\tA\tG\t12345\t-8\r\nrsHigh\tA\tG\t12345\t9\r\nrsMissing\t\t\t\t'
        source.write_bytes(gzip.compress(text.encode()))
        original=source.read_bytes(); summaries=[]; bodies=[]
        for executable in [None,shutil.which('pigz')]:
            if executable is None and summaries: continue
            with patch.object(compression.shutil,'which',return_value=executable), \
                 patch.dict(os.environ,{'PIGZ':'--zip','GZIP':'--best'}):
                target=self.root/'filtered.gz';excluded=self.root/'excluded.gz'
                summaries.append(munging._filter_one_munged_sumstats('T',str(source),str(target),str(excluded),64))
                bodies.append((gzip.decompress(target.read_bytes()),gzip.decompress(excluded.read_bytes())))
        self.assertEqual(source.read_bytes(),original)
        self.assertEqual(bodies[0][0],text.replace('rsHigh\tA\tG\t12345\t9\r\n','').encode()+b'\n')
        self.assertTrue(all(x==bodies[0] for x in bodies))
        self.assertTrue(all(x==summaries[0] for x in summaries))
        self.assertEqual(summaries[0]['variants_removed'],1)
        self.assertEqual(summaries[0]['variants_missing_z'],1)

    def test_compressor_failure_retries_twice_and_preserves_existing_outputs(self):
        fake=self.root/'pigz'
        fake.write_text('#!'+sys.executable+'\nimport sys\nsys.stderr.write("compressor test failure\\n")\nsys.exit(27)\n')
        fake.chmod(0o755)
        source=self.root/'T.sumstats.gz'
        # Large enough to exercise a broken pipe, not just a nonzero final exit.
        source.write_bytes(gzip.compress(('SNP\tA1\tA2\tN\tZ\n'+'rs1\tA\tG\t50000\t1\n'*30000).encode()))
        output=self.root/'out'; folder=output/'ldsc_input_chisq_filtered';folder.mkdir(parents=True)
        targets=[folder/'T.sumstats.gz',folder/'T.chisq_excluded.tsv.gz',output/'LDSC_ChiSquare_Excluded_Variants.tsv.gz']
        for target in targets: target.write_bytes(b'previous valid output')
        with patch.object(compression.shutil,'which',return_value=str(fake)),contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,'failed after 2 attempts') as error:
                munging.filter_munged_sumstats(2,pd.DataFrame({'gwas_name':['T']}),str(self.root),str(output),80)
        self.assertIn('pigz failed (exit 27): compressor test failure',str(error.exception))
        history=pd.read_csv(output/'LDSC_ChiSquare_Worker_Attempts.csv')
        self.assertEqual(history.Attempt.tolist(),[1,2]); self.assertFalse(history.Success.any())
        self.assertTrue(all(p.read_bytes()==b'previous valid output' for p in targets))
        self.assertFalse(list(output.rglob('*.tmp')))

    def test_filter_divides_compression_workers_across_traits(self):
        for name in ('B','A'):
            with open(self.root/(name+'.sumstats.gz'),'wb') as f, gzip.GzipFile(filename='fixture-with-long-name',fileobj=f,mode='wb') as g:
                g.write(b'SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t10000\t1\n')
        with patch.object(munging,'compression_settings',wraps=compression.compression_settings) as settings,contextlib.redirect_stdout(io.StringIO()):
            munging.filter_munged_sumstats(8,pd.DataFrame({'gwas_name':['B','A']}),str(self.root),str(self.root/'out'),80)
        self.assertEqual(settings.call_args_list[0].args,(2,9))
        self.assertEqual(settings.call_args_list[-1].args,(8,9))
        self.assertEqual(pd.read_csv(self.root/'out/LDSC_ChiSquare_Filter_Summary.csv').gwas_name.tolist(),['B','A'])


if __name__=='__main__': unittest.main()
