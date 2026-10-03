"""Exercise INFO-frequency missingness through real bcftools, Bash and awk."""
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import extraction


@unittest.skipUnless(shutil.which('bcftools'), 'local bcftools not installed')
class ExtractionAFTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='AF extraction test ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.filters = {**extraction.DEFAULT_FILTERS, 'max_af_difference': .25}

    def vcf(self, rows, missing_header=None):
        header = '##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=6>\n'
        for name in ('AF', 'EUR', 'EAF'):
            if name != missing_header:
                header += f'##INFO=<ID={name},Number=A,Type=Float,Description="test">\n'
        for name in ('SI', 'AF', 'EZ', 'LP', 'NEF', 'NC', 'NCO'):
            header += f'##FORMAT=<ID={name},Number=1,Type=Float,Description="test">\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait\n'
        body = []
        for i, row in enumerate(rows, 1):
            name, info, *options = row
            opt = options[0] if options else {}
            body.append(f'{opt.get("chrom", 1)}\t{opt.get("pos", i)}\t{name}\tA\t{opt.get("alt", "G")}\t.\tPASS\t{info}\t'
                        f'SI:AF:EZ:LP:NEF:NC:NCO\t{opt.get("si", .9)}:{opt.get("af", .125)}:2:8:360:100:900\n')
        path = self.root/'input file.vcf'
        path.write_text(header + ''.join(body))
        return str(path)

    def extract(self, name, vcf, binary=False, filters=None):
        args = (name, str(self.root), str(self.root), vcf)
        filters = filters or self.filters
        if binary:
            self.assertEqual(extraction.munge_input_worker_with_prevalence(*args, None, filters), .1)
        else:
            self.assertIsNone(extraction.munge_input_worker(*args, filters))
        return pd.read_csv(self.root/f'munge_input/{name}_mungeinput.tsv', sep=r'\s+')

    def audit(self, name):
        return pd.read_csv(self.root/f'munge_input/{name}_AF_Filter_QC.csv').iloc[0].to_dict()

    def test_missing_records_excluded_and_counted_in_both_routes(self):
        rows = [('valid', 'AF=.125;EUR=.125', {'af': .875}),
                ('boundary', 'AF=.375;EUR=.125'), ('zero', 'AF=0;EUR=0'),
                ('difference', 'AF=.875;EUR=.125'),
                ('absent_af', 'EUR=.125'), ('dot_af', 'AF=.;EUR=.125'),
                ('absent_eur', 'AF=.125'), ('dot_eur', 'AF=.125;EUR=.'),
                ('both', '.'), ('other_tag', 'EAF=.125;EUR=.125'),
                ('vector_first', 'AF=.,.125;EUR=.125,.125', {'alt': 'G,T'}),
                ('vector_last', 'AF=.125,.;EUR=.125,.125', {'alt': 'G,T'}),
                ('low_si', 'AF=.125;EUR=.125', {'si': .1}),
                ('low_maf', 'AF=.125;EUR=.125', {'af': .001}),
                ('mhc_missing', 'EUR=.125', {'chrom': 6, 'pos': 30000000})]
        vcf = self.vcf(rows)
        for binary in (False, True):
            name = f'trait_{binary}'
            result = self.extract(name, vcf, binary, {**self.filters, 'exclude_mhc': True})
            self.assertEqual(result.ID.tolist(), ['valid', 'boundary', 'zero'])
            self.assertTrue((result.P == 1e-8).all())
            self.assertEqual(self.audit(name), dict(records_after_mhc=14, missing_info_af=6,
                                                  missing_info_eur=3, excluded_missing_either=8))
            if binary:
                self.assertEqual(result.N_TOTAL.tolist(), [1000]*3)

    def test_valid_inputs_match_previous_filter_exactly(self):
        vcf = self.vcf([('valid', 'AF=.125;EUR=.125', {'af': .875}),
                        ('boundary', 'AF=.375;EUR=.125'), ('zero', 'AF=0;EUR=0'),
                        ('different', 'AF=.875;EUR=.125'),
                        ('palindrome', 'AF=.5;EUR=.5', {'alt': 'T', 'af': .5}),
                        ('low_si', 'AF=.125;EUR=.125', {'si': .1})])
        original = extraction.filter_commands
        def legacy(*args, **kwargs):
            commands = original(*args, **kwargs)
            commands[1][-1] = '(ABS(INFO/AF - INFO/EUR) > 0.25 || INFO/EUR==".")'
            return commands
        for binary in (False, True):
            for remove_palindrome in (False, True):
                filters = {**self.filters, 'remove_palindrome': remove_palindrome}
                new = self.extract('new', vcf, binary, filters)
                with patch.object(extraction, 'filter_commands', side_effect=legacy):
                    old = self.extract('old', vcf, binary, filters)
                pd.testing.assert_frame_equal(new, old)
                self.assertEqual((self.root/'munge_input/new_mungeinput.tsv').read_bytes(),
                                 (self.root/'munge_input/old_mungeinput.tsv').read_bytes())
                self.assertEqual(self.audit('new')['excluded_missing_either'], 0)

    def test_absent_headers_fail_and_remove_stale_audit(self):
        for field, info in [('AF', 'EUR=.125'), ('EUR', 'AF=.125')]:
            for binary in (False, True):
                self.extract('trait', self.vcf([('valid', 'AF=.125;EUR=.125')]), binary)
                vcf = self.vcf([('invalid', info)], missing_header=field)
                with self.assertRaises(subprocess.CalledProcessError) as error:
                    self.extract('trait', vcf, binary)
                self.assertRegex(error.exception.stderr,
                                 f'No such INFO field: {field}|the tag "{field}" is not defined in the VCF header')
                self.assertFalse((self.root/'munge_input/trait_AF_Filter_QC.csv').exists())
                self.assertFalse(list((self.root/'munge_input').glob('.af-qc-*')))

    def test_empty_stream_and_rerun_replace_counts(self):
        self.extract('trait', self.vcf([('missing', 'EUR=.125')]))
        self.assertEqual(self.audit('trait')['excluded_missing_either'], 1)
        result = self.extract('trait', self.vcf([]))
        self.assertTrue(result.empty)
        self.assertTrue(all(value == 0 for value in self.audit('trait').values()))
        self.extract('trait', self.vcf([('valid', 'AF=.125;EUR=.125')]))
        self.assertEqual(self.audit('trait')['records_after_mhc'], 1)
        self.assertEqual(self.audit('trait')['excluded_missing_either'], 0)

    def test_missing_header_retries_twice_then_stops_pipeline(self):
        vcf = self.vcf([('missing', 'EUR=.125')], missing_header='AF')
        manifest = pd.DataFrame({'gwas_name': ['bad'], 'vcf_files': [vcf],
                                 'pop_prevalence': [float('nan')], 'sample_prevalence': [float('nan')]})
        with self.assertRaisesRegex(ValueError, 'VCF extraction jobs failed after 2 attempts'):
            extraction.run_vcf_to_table(1, manifest, str(self.root), self.filters)
        attempts = pd.read_csv(self.root/'LDSC_Extraction_Worker_Attempts.csv')
        self.assertEqual(attempts.Attempt.tolist(), [1, 2])
        self.assertEqual(attempts.Success.tolist(), [False, False])
        self.assertTrue(attempts.Error.str.contains('No such INFO field: AF').all())
        self.assertFalse((self.root/'LDSC_Trait_Prevalence_Metadata.csv').exists())


if __name__ == '__main__':
    unittest.main()
