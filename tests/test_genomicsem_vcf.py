"""GWAS VCF conversion with real bcftools and native CLI contract tests."""
import contextlib
import csv
import gzip
import hashlib
import io
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import polars as pl
from ldsc_gpca import genomicsem_ldsc, genomicsem_vcf

REAL_WHICH = shutil.which


@unittest.skipUnless(shutil.which('bcftools'), 'bcftools not installed')
class NativeVCFTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='native VCF test ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root/'output'
        self.manifest = self.root/'traits.csv'
        self.vcf = self.root/'trait.vcf'
        self.ld = self.root/'ld'
        self.ld.mkdir()
        for name in ('1.l2.ldscore.gz', '1.l2.M_5_50'):
            (self.ld/name).write_text('reference fixture\n')
        self.hm3 = self.root/'hm3.tsv'
        self.hm3.write_text('SNP A1 A2\nrs1 G A\nrs2 T C\n')
        self.write_vcf()
        self.write_manifest()

    def write_vcf(self, *, missing=None, samples=1, empty=False):
        columns = ['AF', 'ES', 'SE', 'LP', 'NEF', 'SI']
        if missing:
            columns.remove(missing)
        lines = ['##fileformat=VCFv4.2', '##contig=<ID=1>']
        lines += [f'##FORMAT=<ID={c},Number=1,Type=Float,Description="fixture">' for c in columns]
        lines += ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                  '\t'.join(f'sample{i}' for i in range(samples))]
        # Low but valid INFO is retained here: native munging applies its threshold.
        records = [('rs1', 1, '.8', '-0.123456789', '.04', '4', '20000', '.95'),
                   ('rs2', 2, '.2', '.25', '.05', '400', '22000', '.8'),
                   ('rs2', 3, '.2', '.1', '.05', '1', '22000', '.9'),
                   ('rs_missing_info', 4, '.2', '.1', '.05', '1', '22000', '.'),
                   ('rs_bad_af', 5, '1.2', '.1', '.05', '1', '22000', '.9'),
                   ('rs_bad_se', 6, '.2', '.1', '0', '1', '22000', '.9'),
                   ('rs_bad_info', 7, '.2', '.1', '.05', '1', '22000', '1.2')]
        for snp, pos, *values in ([] if empty else records):
            values = dict(zip(['AF','ES','SE','LP','NEF','SI'], values))
            sample = ':'.join(values[c] for c in columns)
            lines.append(f'1\t{pos}\t{snp}\tA\tG\t.\tPASS\t.\t' + ':'.join(columns) + '\t' +
                         '\t'.join([sample]*samples))
        self.vcf.write_text('\n'.join(lines)+'\n')

    def write_manifest(self, *, binary=False, n=None, vcf=None):
        columns = ['traitname','vcf_files','sample_prevalence','population_prevalence','N']
        with self.manifest.open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            for name in ('002', '001'):
                writer.writerow([name, (vcf or self.vcf).name, '.5' if binary else '',
                                 '.1' if binary else '', '' if n is None else n])

    def args(self, *extra):
        return ['--input', str(self.manifest), '--outdir', str(self.out), '--ld_ref', str(self.ld),
                '--hm3', str(self.hm3), '--chromosomes', '1', '--vcf_input', *extra]

    def resolve(self, *extra):
        return genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(self.args(*extra)))

    def prepare(self, *extra):
        rows, mode = self.resolve(*extra)
        self.assertEqual(mode, 'vcf')
        self.out.mkdir()
        return genomicsem_vcf.prepare_vcf_inputs(rows, self.out, shutil.which('bcftools'), 2, 1e-300)

    def test_quantitative_values_order_precision_qc_and_source_preserved(self):
        digest = hashlib.sha256(self.vcf.read_bytes()).hexdigest()
        rows = self.prepare()
        self.assertEqual([r['traitname'] for r in rows], ['002','001'])
        self.assertTrue(all(r['N'] == 'NA' for r in rows))
        frame = pl.read_csv(rows[0]['source_file'], separator='\t')
        self.assertEqual(frame.columns, ['SNP','CHR','POS','A1','A2','MAF','beta','se','N','p','INFO'])
        self.assertEqual(frame['SNP'].to_list(), ['rs1','rs2'])
        self.assertEqual(frame['A1'].to_list(), ['G','G'])
        self.assertEqual(frame['A2'].to_list(), ['A','A'])
        self.assertEqual(frame['beta'][0], -.123457)  # bcftools' decoded float; no further rounding
        self.assertAlmostEqual(frame['MAF'][0], .2, places=7)
        self.assertEqual(frame['N'].to_list(), [20000,22000])
        self.assertEqual(frame['p'].to_list(), [1e-4,1e-300])
        self.assertAlmostEqual(frame['INFO'][1], .8, places=7)
        qc = pl.read_csv(self.out/'GenomicSEM_VCF_QC_Summary.csv', schema_overrides={'traitname':pl.String})
        self.assertEqual(qc['traitname'].to_list(), ['002','001'])
        self.assertEqual(qc['input_rows'].to_list(), [7,7])
        self.assertEqual(qc['retained_rows'].to_list(), [2,2])
        self.assertEqual(qc['p_adjusted_rows'].to_list(), [1,1])
        self.assertEqual(qc['N_source'].to_list(), ['FORMAT/NEF']*2)
        issues = pl.read_csv(self.out/'GenomicSEM_VCF_QC_Issues.csv')
        self.assertTrue(any('INFO:' in value for value in issues['QC_reason']))
        self.assertTrue(any('Duplicate SNP' in value for value in issues['QC_reason']))
        self.assertEqual(hashlib.sha256(self.vcf.read_bytes()).hexdigest(), digest)
        self.assertFalse((self.out/'gpca_inputs').exists())

    def test_binary_explicit_n_does_not_require_or_use_nef(self):
        self.write_vcf(missing='NEF')
        self.write_manifest(binary=True, n=35000)
        rows = self.prepare()
        self.assertEqual(rows[0]['sampleprevalence'], .5)
        self.assertEqual(rows[0]['populationprevalence'], .1)
        self.assertEqual(pl.read_csv(rows[0]['source_file'], separator='\t')['N'].to_list(), [35000]*2)
        self.assertEqual(pl.read_csv(self.out/'GenomicSEM_VCF_QC_Summary.csv')['N_source'].to_list(), ['manifest_N']*2)

    def test_binary_missing_n_or_partial_prevalence_rejected(self):
        self.write_manifest(binary=True)
        with self.assertRaisesRegex(ValueError, 'binary VCF input requires.*N'):
            self.resolve()
        text = self.manifest.read_text().replace(',0.1,', ',,').replace(',.1,', ',,')
        self.manifest.write_text(text)
        with self.assertRaisesRegex(ValueError, 'both prevalences'):
            self.resolve()

    def test_mutually_exclusive_modes_and_explicit_selection(self):
        parser = genomicsem_ldsc.build_parser()
        for flags in (['--munged_input'], ['--munged_dir','data']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                parser.parse_args(self.args(*flags))
        args = self.args(); args.remove('--vcf_input')
        with self.assertRaisesRegex(ValueError, 'sumstats_file'):
            genomicsem_ldsc.resolve_manifest(parser.parse_args(args))

    def test_invalid_sample_sizes_rejected_before_output(self):
        for n in (0,-1,'NaN','inf'):
            self.write_manifest(binary=True,n=n)
            with self.assertRaises(ValueError): self.resolve()
            self.assertFalse(self.out.exists())

    def test_invalid_p_floors_rejected_before_conversion(self):
        for value in ('0','1','-1','nan','inf'):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                genomicsem_ldsc.main(self.args('--p_min',value))
            self.assertFalse(self.out.exists())

    def test_gzip_vcf_and_quantitative_n_override(self):
        compressed = self.root/'trait.vcf.gz'
        with gzip.open(compressed,'wt') as stream: stream.write(self.vcf.read_text())
        self.write_manifest(n=12345,vcf=compressed)
        rows = self.prepare()
        self.assertEqual(pl.read_csv(rows[0]['source_file'], separator='\t')['N'].to_list(), [12345]*2)

    def test_failed_fields_empty_and_multiple_samples_retry_without_publishing(self):
        for i, options in enumerate(({'missing':'SI'}, {'missing':'NEF'}, {'samples':2}, {'empty':True})):
            with self.subTest(options=options):
                self.out = self.root/f'failed{i}'
                self.write_vcf(**options)
                with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError,'failed after 2 attempts'):
                    self.prepare()
                attempts = pl.read_csv(self.out/'GenomicSEM_VCF_Worker_Attempts.csv')
                self.assertEqual(attempts.height,4)
                self.assertEqual(attempts['Success'].to_list(),[False]*4)
                self.assertFalse((self.out/'vcf_input').exists())

    def test_cli_handoff_uses_converted_table_and_unchanged_native_mode(self):
        run = subprocess.run
        calls = []
        def dispatch(command, *args, **kwargs):
            if len(command)>1 and str(command[1]).endswith('ldsc_entry.R'):
                calls.append(command)
                return subprocess.CompletedProcess(command,0)
            return run(command,*args,**kwargs)
        # Only native R execution is stubbed; bcftools and table conversion are real.
        with patch.object(genomicsem_ldsc.shutil,'which',side_effect=lambda name: '/verified/Rscript' if name=='Rscript' else REAL_WHICH(name)), \
             patch.object(genomicsem_ldsc.subprocess,'run',side_effect=dispatch):
            self.assertEqual(genomicsem_ldsc.main(self.args('--n_cores','2')),0)
        self.assertEqual(len(calls),1)
        command = calls[0]
        self.assertEqual(command[7],'munge')
        self.assertEqual(command[6],str(self.hm3.resolve()))
        with (self.out/'Resolved_Manifest.csv').open() as stream: rows=list(csv.DictReader(stream))
        self.assertEqual([r['traitname'] for r in rows],['002','001'])
        self.assertTrue(all('/vcf_input/' in r['source_file'] for r in rows))

    def test_failed_vcf_and_existing_output_never_start_r(self):
        self.write_vcf(missing='SI')
        run = subprocess.run
        def dispatch(command,*args,**kwargs):
            self.assertFalse(len(command)>1 and str(command[1]).endswith('ldsc_entry.R'))
            return run(command,*args,**kwargs)
        with patch.object(genomicsem_ldsc.shutil,'which',side_effect=lambda name: '/verified/Rscript' if name=='Rscript' else REAL_WHICH(name)), \
             patch.object(genomicsem_ldsc.subprocess,'run',side_effect=dispatch), \
             contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit): genomicsem_ldsc.main(self.args())
            with self.assertRaises(SystemExit): genomicsem_ldsc.main(self.args())

if __name__ == '__main__':
    unittest.main()
