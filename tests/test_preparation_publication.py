"""Fault injection at each preparation publication boundary, without a munger runtime."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ldsc_gpca import prepare_ldsc
from ldsc_gpca.interfaces import read_manifest
from ldsc_gpca.preparation_provenance import validate_completed_bundle


class PublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='publication test ')
        self.addCleanup(temporary.cleanup)
        self.out = Path(temporary.name).resolve()
        self.settings = self.out/'Preparation_Settings.json'
        self.settings.write_text('{}')
        (self.out/'munge_inputs').mkdir()
        for name in ('A', 'B'):
            (self.out/'munge_inputs'/f'{name}_munge_inputs.tsv').write_text('raw data retained\n')
        self.manifest = self.out/'Prepared_LDSC_Manifest.csv'
        metadata = {name: {'N': None, 'sample_prevalence': None, 'population_prevalence': None}
                    for name in ('A', 'B')}
        prepare_ldsc.write_ldsc_manifest(self.manifest, [('A', None), ('B', None)], metadata, self.out)
        self.original_manifest = self.manifest.read_bytes()
        self.opts = SimpleNamespace(outdir=self.out, munge_backend='python', info_filter=.9,
                                    maf_filter=.01, hapmap_file=self.out/'reference', prepare_workers=1)
        script = self.out/'munger.py'
        script.write_text('import gzip, sys\n'
                          'path = sys.argv[sys.argv.index("--out")+1]+".sumstats.gz"\n'
                          'with gzip.open(path, "wt") as f: f.write("SNP\\tN\\tZ\\tA1\\tA2\\nrs1\\t1000\\t2\\tA\\tG\\n")\n')
        self.command = [sys.executable, str(script)]

    def run_munging(self):
        with contextlib.redirect_stdout(io.StringIO()):
            prepare_ldsc.munge_prepared(self.opts, self.command)

    def assert_failed(self):
        self.assertEqual(json.loads(self.settings.read_text())['munging_status'], 'failed')
        for name in ('A', 'B'):
            self.assertEqual((self.out/'munge_inputs'/f'{name}_munge_inputs.tsv').read_text(), 'raw data retained\n')
            self.assertTrue((self.out/'munge_logs'/f'{name}.console.log').exists())
        self.assertFalse(list(self.out.glob('.manifest-*')))
        self.assertFalse(list(self.out.glob('.ldsc-checkpoint-*')))
        if (self.out/'munged').exists():
            sidecar = self.out/'munged/A.prevalence.json'
            with self.assertRaisesRegex(ValueError, 'not a completed'):
                validate_completed_bundle(self.out/'munged/A.sumstats.gz', sidecar, json.loads(sidecar.read_text()))

    def test_success_records_completion_after_files_and_manifest(self):
        original = prepare_ldsc.write_checkpoint
        def check(path, settings):
            if settings['munging_status'] == 'completed':
                _, rows = read_manifest(self.manifest)
                self.assertEqual(len(rows), 2)
                self.assertTrue(all(Path(row['munged_file']).is_file() for row in rows))
            return original(path, settings)
        with patch.object(prepare_ldsc, 'write_checkpoint', side_effect=check):
            self.run_munging()
        self.assertEqual(json.loads(self.settings.read_text())['munging_status'], 'completed')
        sidecar = self.out/'munged/A.prevalence.json'
        validate_completed_bundle(self.out/'munged/A.sumstats.gz', sidecar, json.loads(sidecar.read_text()))

    def test_directory_publication_failure_never_marks_completed(self):
        original = Path.rename
        def fail(source, target):
            if Path(target).resolve() == self.out/'munged':
                self.assertEqual(json.loads(self.settings.read_text())['munging_status'], 'running')
                raise OSError('directory publication failed')
            return original(source, target)
        with patch.object(Path, 'rename', fail), self.assertRaisesRegex(OSError, 'directory publication failed'):
            self.run_munging()
        self.assert_failed()
        self.assertFalse((self.out/'munged').exists())
        self.assertEqual(self.manifest.read_bytes(), self.original_manifest)

    def test_manifest_publication_failure_keeps_original_manifest(self):
        original = prepare_ldsc.os.replace
        def fail(source, target):
            if Path(target).resolve() == self.manifest:
                raise OSError('manifest publication failed')
            return original(source, target)
        with patch.object(prepare_ldsc.os, 'replace', side_effect=fail), self.assertRaisesRegex(OSError, 'manifest publication failed'):
            self.run_munging()
        self.assert_failed()
        self.assertTrue((self.out/'munged').is_dir())
        self.assertEqual(self.manifest.read_bytes(), self.original_manifest)

    def test_completion_publication_failure_is_a_failed_run(self):
        original = prepare_ldsc.os.replace
        def fail(source, target):
            if Path(target).resolve() == self.settings and json.loads(Path(source).read_text())['munging_status'] == 'completed':
                raise OSError('completion publication failed')
            return original(source, target)
        with patch.object(prepare_ldsc.os, 'replace', side_effect=fail), self.assertRaisesRegex(OSError, 'completion publication failed'):
            self.run_munging()
        self.assert_failed()
        _, rows = read_manifest(self.manifest)
        self.assertTrue(all(Path(row['munged_file']).is_file() for row in rows))


if __name__ == '__main__':
    unittest.main()
