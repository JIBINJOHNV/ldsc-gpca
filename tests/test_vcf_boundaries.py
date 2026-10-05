"""Sample ambiguity must fail before any extracted files are replaced."""
import contextlib
import io
from pathlib import Path
import shutil
import unittest

import pandas as pd
from ldsc_gpca import extraction, prepare, genomicsem_vcf
import test_extraction_af as extraction_fixture
import test_genomicsem_vcf as preparation_fixture


def set_samples(path, count):
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith('##'):
            lines.append(line)
            continue
        cells = line.split('\t')
        if count == 0:
            cells = cells[:8]  # Valid sites-only VCF.
        else:
            cells = cells[:10] + (['different_GWAS'] if line.startswith('#') else
                                  ['.9:.2:-6:10:50000:200:800'])
        lines.append('\t'.join(cells))
    path.write_text('\n'.join(lines) + '\n')


@unittest.skipUnless(shutil.which('bcftools'), 'bcftools not installed')
class VCFBoundaryTests(unittest.TestCase):
    def test_original_routes_reject_zero_and_multiple_samples_before_writes(self):
        case = extraction_fixture.ExtractionAFTests(); case.setUp(); self.addCleanup(case.doCleanups)
        for count in (0, 2):
            for binary in (False, True):
                with self.subTest(samples=count, binary=binary):
                    path = Path(case.vcf([('rs1', 'AF=.125;EUR=.125')]))
                    # An existing valid result and sidecar must stay byte-identical.
                    with contextlib.redirect_stdout(io.StringIO()):
                        case.extract('saved', str(path), binary)
                    files = list((case.root/'munge_input').iterdir())
                    before = {p.name: p.read_bytes() for p in files}
                    set_samples(path, count)
                    with self.assertRaisesRegex(ValueError, f'exactly one GWAS sample; found {count}'):
                        case.extract('saved', str(path), binary)
                    self.assertEqual(before, {p.name: p.read_bytes() for p in files})
                    fresh = case.root/f'fresh-{count}-{binary}'
                    args = ('fresh', str(case.root), str(fresh), str(path))
                    with self.assertRaises(ValueError):
                        if binary:
                            extraction.munge_input_worker_with_prevalence(*args, None)
                        else:
                            extraction.munge_input_worker(*args)
                    self.assertFalse(fresh.exists())

    def test_original_controller_records_failure_and_never_publishes_prevalences(self):
        case = extraction_fixture.ExtractionAFTests(); case.setUp(); self.addCleanup(case.doCleanups)
        path = Path(case.vcf([('rs1', 'AF=.125;EUR=.125')]))
        set_samples(path, 2)
        manifest = pd.DataFrame({'gwas_name': ['A'], 'vcf_files': [str(path)],
            'pop_prevalence': [float('nan')], 'sample_prevalence': [float('nan')]})
        with self.assertRaisesRegex(ValueError, 'failed after 2 attempts'):
            extraction.run_vcf_to_table(1, manifest, str(case.root))
        attempts = pd.read_csv(case.root/'LDSC_Extraction_Worker_Attempts.csv')
        self.assertEqual(attempts.Attempt.tolist(), [1, 2])
        self.assertTrue((~attempts.Success).all())
        self.assertTrue(attempts.Error.str.contains('found 2').all())
        self.assertFalse((case.root/'munge_input').exists())
        self.assertFalse((case.root/'LDSC_Trait_Prevalence_Metadata.csv').exists())

    def test_preparation_and_genomicsem_routes_share_single_sample_contract(self):
        case = preparation_fixture.NativeVCFTests(); case.setUp(); self.addCleanup(case.doCleanups)
        case.hm3.write_text('SNP\tA1\tA2\nrs1\tG\tA\n')
        for count in (0, 2):
            case.write_vcf(); set_samples(case.vcf, count)
            for mode in ('gpca', 'ldsc', 'both'):
                out = case.root/f'{mode}-{count}'
                with self.subTest(samples=count, mode=mode), contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(ValueError, f'found {count}'):
                        prepare.prepare_inputs(case.manifest, out, splitby_chr='nosplit', mode=mode,
                            hapmap_file=case.hm3 if mode != 'gpca' else None, prepare_workers=1)
                    self.assertFalse((out/'gpca_inputs').exists())
                    self.assertFalse((out/'munge_inputs').exists())
            out = case.root/f'genomicsem-{count}'; out.mkdir()
            with self.assertRaisesRegex(ValueError, f'found {count}'):
                genomicsem_vcf.prepare_vcf_inputs([{'traitname': 'A', 'source_file': str(case.vcf), 'N': 'NA'}],
                                                 out, shutil.which('bcftools'), 1, 1e-300)
            self.assertFalse((out/'vcf_input').exists())
