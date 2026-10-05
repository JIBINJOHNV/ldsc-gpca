"""Opt-in full installed pinned-package preparation, alignment and regression test.

Set LDSC_GPCA_TEST_RSCRIPT to an Rscript (or wrapper) with the complete pinned
GenomicSEM package. The estimator and munger are never mocked or modified.
"""
import contextlib
import csv
import gzip
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import polars as pl
from scipy.stats import chi2

from genomicsem_fixtures import seal_bundle
from ldsc_gpca import genomicsem_ldsc, prepare
from ldsc_gpca.genomicsem_inputs import GENOMICSEM_COMMIT, InputCompatibilityError
from ldsc_gpca.restart import file_digest


@unittest.skipUnless(os.environ.get('LDSC_GPCA_TEST_RSCRIPT') and shutil.which('bcftools'),
                     'requires bcftools and LDSC_GPCA_TEST_RSCRIPT with full pinned GenomicSEM')
class FullGenomicSEMReuseTests(unittest.TestCase):
    def test_preparation_retained_snps_signed_products_estimates_and_ses(self):
        rscript = os.environ['LDSC_GPCA_TEST_RSCRIPT']
        with tempfile.TemporaryDirectory(prefix='full GenomicSEM ') as temporary:
            root = Path(temporary).resolve()
            count, sample_size = 2400, 20000
            rng = np.random.default_rng(731)
            ld_scores = rng.uniform(5, 50, count)
            z_a = rng.normal(size=count)*np.sqrt(1+.04*ld_scores)
            z_b = .4*z_a + rng.normal(size=count)*np.sqrt(.84+.05*ld_scores)
            snps = [f'rs{i+1}' for i in range(count)]
            hm3 = root/'reference.tsv'
            hm3.write_text('SNP\tA1\tA2\n'+''.join(f'{snp}\tA\tG\n' for snp in snps))
            for name, values in [('A', z_a), ('B', z_b)]:
                lines = ['##fileformat=VCFv4.2', '##contig=<ID=1>']
                columns = ['AF', 'ES', 'SE', 'LP', 'NEF', 'SI']
                lines += [f'##FORMAT=<ID={c},Number=1,Type=Float,Description="fixture">' for c in columns]
                lines += ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS']
                for i, (snp, z) in enumerate(zip(snps, values)):
                    kind = i % 4 if name == 'B' else 0
                    a1, a2 = [('A','G'), ('G','A'), ('T','C'), ('A','C')][kind]
                    effect = -.01*z if kind == 1 else .01*z
                    lp = -np.log10(chi2.sf(z*z, 1))
                    lines.append(f'1\t{i+1}\t{snp}\t{a2}\t{a1}\t.\tPASS\t.\tAF:ES:SE:LP:NEF:SI\t.2:{effect:.12g}:.01:{lp:.12g}:{sample_size}:1')
                (root/f'{name}.vcf').write_text('\n'.join(lines)+'\n')
            manifest = root/'input.csv'
            manifest.write_text('traitname,vcf_files,sample_prevalence,population_prevalence\nA,A.vcf,,\nB,B.vcf,,\n')
            bundle = root/'prepared'
            with contextlib.redirect_stdout(io.StringIO()):
                status = prepare.main(['--input', str(manifest), '--outdir', str(bundle), '--mode', 'ldsc',
                    '--munge_backend', 'genomicsem', '--rscript', rscript, '--hm3', str(hm3), '--n_cores', '1'])
            self.assertEqual(status, 0)
            frame_a = pl.read_csv(bundle/'munged/A.sumstats.gz', separator='\t')
            frame_b = pl.read_csv(bundle/'munged/B.sumstats.gz', separator='\t')
            retained = [snp for i, snp in enumerate(snps) if i % 4 < 2]
            self.assertEqual(set(frame_a['SNP']), set(snps))
            self.assertEqual(set(frame_b['SNP']), set(retained))
            # Use actual prepared P text, since bcftools query rounds FORMAT floats.
            # This independently checks qchisq(P) and the reference-order sign rule.
            raw_a = pl.read_csv(bundle/'munge_inputs/A_munge_inputs.tsv', separator='\t')
            raw_b = pl.read_csv(bundle/'munge_inputs/B_munge_inputs.tsv', separator='\t')
            expected_a = np.sign(raw_a['BETA'])*np.sqrt(chi2.isf(raw_a['P'], 1))
            expected_b = np.sign(raw_b['BETA'])*np.sqrt(chi2.isf(raw_b['P'], 1))*np.where(raw_b['A1']=='G', -1, 1)
            expected = raw_a.select('SNP').with_columns(pl.Series('za', expected_a)).join(
                raw_b.select('SNP').with_columns(pl.Series('zb', expected_b)), on='SNP')
            products = frame_a.join(frame_b, on='SNP', suffix='_b').join(expected, on='SNP')
            np.testing.assert_allclose(products['Z']*products['Z_b'], products['za']*products['zb'], rtol=1e-12, atol=1e-12)
            self.assertLess(float((frame_a['Z']**2).max()), 80)
            self.assertLess(float((frame_b['Z']**2).max()), 80)
            self.assertTrue((frame_b['A1'] == 'A').all() and (frame_b['A2'] == 'G').all())

            ld = root/'ld'; ld.mkdir()
            with gzip.open(ld/'1.l2.ldscore.gz', 'wt') as stream:
                stream.write('CHR\tSNP\tBP\tL2\n')
                stream.writelines(f'1\t{snp}\t{i+1}\t{score:.17g}\n' for i, (snp, score) in enumerate(zip(snps, ld_scores)))
            (ld/'1.l2.M_5_50').write_text('100000\n')
            output = root/'regression'
            args = ['--input', str(bundle/'Prepared_LDSC_Manifest.csv'), '--outdir', str(output),
                    '--ld_ref', str(ld), '--munged_dir', str(bundle/'munged'), '--rscript', rscript,
                    '--chromosomes', '1', '--n_blocks', '20', '--chisq_max', '80', '--invalid_h2_action', 'error']
            before = {p.name: file_digest(p) for p in (bundle/'munged').iterdir()}
            self.assertEqual(genomicsem_ldsc.main(args), 0)
            comparison = root/'compare.R'
            comparison.write_text('''a <- commandArgs(TRUE)
stopifnot(packageDescription("GenomicSEM")$RemoteSha == a[4])
setwd(a[3])
paths <- file.path(a[1], "munged", paste0(c("A", "B"), ".sumstats.gz"))
direct <- GenomicSEM::ldsc(traits=paths, sample.prev=c(NA,NA), population.prev=c(NA,NA),
  trait.names=c("A","B"), ld=a[2], wld=a[2], sep_weights=FALSE,
  chr=1, n.blocks=20, chisq.max=80, stand=TRUE, ldsc.log="direct")
load(file.path(a[3], "genomicsem_LDSC.RData"))
for (field in c("S","I","V","S_Stand","V_Stand")) {
  # The wrapper adds missing matrix dimnames; numerical estimates are untouched.
  stopifnot(all(is.finite(LDSCoutput[[field]])),
            identical(dim(LDSCoutput[[field]]), dim(direct[[field]])),
            identical(as.numeric(LDSCoutput[[field]]), as.numeric(direct[[field]])))
}
stopifnot(all(diag(direct$V)>0), all(diag(direct$V_Stand)>0),
          identical(sqrt(diag(LDSCoutput$V)),sqrt(diag(direct$V))),
          identical(sqrt(diag(LDSCoutput$V_Stand)),sqrt(diag(direct$V_Stand))))
cat("Full pinned GenomicSEM estimates, intercepts, covariance matrices and SEs match exactly.\\n")
''')
            result = subprocess.run([rscript, str(comparison), str(bundle), str(ld), str(output), GENOMICSEM_COMMIT],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertEqual(before, {p.name: file_digest(p) for p in (bundle/'munged').iterdir()})
            with (output/'GenomicSEM_Input_QC.csv').open() as stream:
                audit = list(csv.DictReader(stream))
            self.assertEqual([int(row['n_snps']) for row in audit], [2400, 1200])

            # Full-package outputs still pass through the independent strict gate.
            path = bundle/'munged/B.sumstats.gz'
            original = path.read_bytes()
            for pair in [('G','A'), ('T','C'), ('A','C'), None]:
                with gzip.open(path, 'rt') as stream: lines = stream.readlines()
                fields = lines[1].rstrip().split('\t')
                if pair is None:
                    lines.append(lines[1])
                else:
                    header = lines[0].rstrip().split('\t')
                    fields[header.index('A1')], fields[header.index('A2')] = pair
                    lines[1] = '\t'.join(fields)+'\n'
                with gzip.open(path, 'wt') as stream: stream.writelines(lines)
                seal_bundle(bundle)
                with self.assertRaises(InputCompatibilityError):
                    genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(args))
                path.write_bytes(original)
                seal_bundle(bundle)
            sidecar = bundle/'munged/B.prevalence.json'
            original_sidecar = sidecar.read_text()
            for field, value in [('sha256', 'corrupt'), ('sample_prevalence', .2),
                                 ('population_prevalence', .05), ('N_convention', 'total')]:
                metadata = json.loads(original_sidecar)
                metadata[field] = value
                sidecar.write_text(json.dumps(metadata))
                with self.assertRaises(InputCompatibilityError):
                    genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(args))
            sidecar.write_text(original_sidecar)


if __name__ == '__main__': unittest.main()
