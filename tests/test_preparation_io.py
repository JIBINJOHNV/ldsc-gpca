"""Shared VCF scans preserve original records and independent failure audits."""
import csv
import gzip
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import polars as pl
from ldsc_gpca.vcf_common import write_original_issue_reports


class OriginalIssueReportsTests(unittest.TestCase):
    def test_one_scan_preserves_text_overlapping_reasons_and_empty_reports(self):
        header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t001'
        rows = [f'1\t00{i}\trs{i}\tA\tG\t.\tPASS\tNOTE=a,b\tES:SE\t0.0100:0e0' for i in range(4)]
        text = '##fileformat=VCFv4.2\r\n'+header+'\r\n'+'\r\n'.join(rows)
        groups = [([(0, 'bad SE; bad Z', 'removed'), (3, 'P floor', 'p_adjusted')], 'NA'),
                  ([(1, 'INFO missing', 'removed'), (3, 'not HapMap', 'removed')], '001'), ([], 'empty')]
        for compressed in (False, True):
            with self.subTest(compressed=compressed), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                source = root/'original.vcf'
                source.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
                before = source.read_bytes()
                reports = []
                for i, (issues, trait) in enumerate(groups):
                    frame = pl.DataFrame(issues, schema={'_row':pl.UInt32, 'QC_reason':pl.String, 'QC_action':pl.String}, orient='row')
                    reports.append((root/f'{i}.csv', trait, frame))
                with patch('ldsc_gpca.vcf_common.gzip.open', wraps=gzip.open) as opened:
                    write_original_issue_reports(source, reports)
                self.assertEqual(opened.call_count, int(compressed))
                for i, (issues, trait) in enumerate(groups):
                    expected = io.StringIO(newline='')
                    writer = csv.writer(expected)
                    writer.writerow(header.split('\t')+['traitname','QC_action','QC_reason'])
                    writer.writerows(rows[row].split('\t')+[trait,action,reason] for row,reason,action in issues)
                    self.assertEqual((root/f'{i}.csv').read_bytes(), expected.getvalue().encode())
                self.assertEqual(source.read_bytes(), before)

    def test_incomplete_original_fails_and_retains_available_reports(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root/'original.vcf'
            source.write_text('#CHROM\tPOS\n1\t1\n')
            issues = lambda row: pl.DataFrame({'_row':[row], 'QC_reason':['bad'], 'QC_action':['removed']})
            with self.assertRaisesRegex(ValueError, 'ended before all affected'):
                write_original_issue_reports(source, [(root/'a.csv', 'A', issues(0)), (root/'b.csv', 'B', issues(2))])
            self.assertIn('1,1,A,removed,bad', (root/'a.csv').read_text())
            self.assertEqual(len((root/'b.csv').read_text().splitlines()), 1)

    def test_header_collision_remains_fatal_for_all_reports(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root/'original.vcf'
            source.write_text('#CHROM\ttraitname\n1\tA\n')
            empty = pl.DataFrame(schema={'_row':pl.UInt32, 'QC_reason':pl.String, 'QC_action':pl.String})
            with self.assertRaisesRegex(ValueError, 'column names must be unique'):
                write_original_issue_reports(source, [(root/'a.csv', 'A', empty), (root/'b.csv', 'B', empty)])


if __name__ == '__main__':
    unittest.main()
