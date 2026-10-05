"""Scientific input compatibility and immutable preparation-bundle boundaries."""
import contextlib
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from genomicsem_fixtures import make_bundle, seal_bundle
from ldsc_gpca import genomicsem_ldsc
from ldsc_gpca.genomicsem_inputs import InputCompatibilityError, REFERENCE_FILE, read_allele_reference, validate_alleles
from ldsc_gpca.restart import file_digest, write_checkpoint


class GenomicSEMInputTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='GenomicSEM reuse ')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        make_bundle(self.root)
        self.manifest = self.root/'Prepared_LDSC_Manifest.csv'
        self.path = self.root/'munged/B.sumstats.gz'
        self.sidecar = self.root/'munged/B.prevalence.json'
        self.settings = self.root/'Preparation_Settings.json'

    def args(self):
        return ['--input', str(self.manifest), '--outdir', str(self.root/'out'), '--ld_ref', 'unused',
                '--munged_dir', str(self.root/'munged')]

    def resolve(self):
        return genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(self.args()))

    def change_rows(self, rows):
        with gzip.open(self.path, 'wt') as stream:
            stream.write('SNP\tN\tZ\tA1\tA2\n'+rows)
        seal_bundle(self.root)  # Deliberately pass provenance to exercise the independent SNP scan.

    def test_valid_preserves_names_order_bytes_and_signed_products(self):
        before = {p.name: file_digest(p) for p in (self.root/'munged').iterdir()}
        rows, mode = self.resolve()
        self.assertEqual(mode, 'existing')
        self.assertEqual([r['traitname'] for r in rows], ['A', 'B'])
        self.assertEqual([r['_reuse_audit']['n_snps'] for r in rows], [2, 2])
        self.assertEqual(rows[0]['_reuse_audit']['snp_order_sha256'], hashlib.sha256(b'rs1\nrs2\n').hexdigest())
        data = []
        for row in rows:
            with gzip.open(row['source_file'], 'rt') as stream:
                data.append(list(csv.DictReader(stream, delimiter='\t')))
        self.assertEqual([float(a['Z'])*float(b['Z']) for a, b in zip(*data)], [6, 4])
        self.assertEqual(before, {p.name: file_digest(p) for p in (self.root/'munged').iterdir()})

    def test_all_unexpected_allele_orientations_fail_without_edits(self):
        for pair, reason in [('G\tA', 'swapped allele'), ('T\tC', 'strand-complemented'),
                             ('C\tT', 'swapped strand-complemented'), ('A\tC', 'incompatible'),
                             ('a\tg', 'incompatible'), ('A\tA', 'incompatible')]:
            with self.subTest(pair=pair):
                self.change_rows(f'rs1\t1000\t3\t{pair}\n')
                before = self.path.read_bytes()
                with self.assertRaisesRegex(InputCompatibilityError, reason) as caught:
                    self.resolve()
                self.assertEqual(caught.exception.audit['SNP'], 'rs1')
                self.assertEqual(caught.exception.audit['row'], 2)
                self.assertEqual(self.path.read_bytes(), before)

    def test_duplicates_unknown_snp_and_missing_values_are_fatal(self):
        for rows, reason in [
            ('rs1\t1000\t3\tA\tG\n'*2, 'duplicate SNP'),
            ('unknown\t1000\t3\tA\tG\n', 'absent from allele reference'),
            ('', 'no usable SNP'), ('rs1\t1000\tNA\tA\tG\n', 'nonnumeric'),
            ('rs1\t1000\tinf\tA\tG\n', 'invalid N/Z'),
            ('rs1\t0\t3\tA\tG\n', 'invalid N/Z'),
            ('rs1\tnan\t3\tA\tG\n', 'invalid N/Z'),
            ('rs1\t1000\t3\tA\n', 'malformed row')]:
            with self.subTest(reason=reason):
                self.change_rows(rows)
                with self.assertRaisesRegex(InputCompatibilityError, reason):
                    self.resolve()

    def test_changed_file_sidecar_reference_and_manifest_checksums_fail(self):
        for path in (self.path, self.sidecar, self.root/'munged'/REFERENCE_FILE, self.manifest):
            with self.subTest(file=path.name):
                before = path.read_bytes()
                path.write_bytes(before.replace(b',yes,', b',no,') if path == self.manifest else before+b' ')
                try:
                    with self.assertRaisesRegex(InputCompatibilityError, 'checksum mismatch'):
                        self.resolve()
                finally:
                    path.write_bytes(before)

    def test_backend_prevalences_n_and_runtime_metadata_conflicts_fail(self):
        original = json.loads(self.sidecar.read_text())
        for field, value, reason in [('munge_backend', 'python', 'backend'),
            ('sample_prevalence', .2, 'sample_prevalence'), ('population_prevalence', .05, 'population_prevalence'),
            ('N_source', 'manifest_N', 'N convention'), ('N_convention', 'total', 'N convention'),
            ('sample_size_column', 'N_TOTAL', 'N convention'), ('N_override', -1, 'N_override'),
            ('sha256', 'corrupt', 'checksum'), ('runtime', {}, 'runtime'), ('traitname', 'C', 'traitname')]:
            with self.subTest(field=field):
                write_checkpoint(self.sidecar, {**original, field: value})
                with self.assertRaisesRegex(InputCompatibilityError, reason):
                    self.resolve()
        write_checkpoint(self.sidecar, original)

    def test_missing_or_incomplete_bundle_cannot_be_reused(self):
        for path in (self.sidecar, self.settings, self.root/'munged'/REFERENCE_FILE):
            with self.subTest(path=path):
                backup = path.with_name(path.name+'.backup')
                path.rename(backup)
                try:
                    with self.assertRaises(InputCompatibilityError): self.resolve()
                finally:
                    backup.rename(path)
        settings = json.loads(self.settings.read_text())
        for status in ('running', 'failed', None):
            write_checkpoint(self.settings, {**settings, 'munging_status': status})
            with self.assertRaisesRegex(InputCompatibilityError, 'not a completed'):
                self.resolve()

    def test_binary_metadata_and_actual_n_agree(self):
        make_bundle(self.root, binary_traits=('B',))
        self.resolve()
        self.change_rows('rs1\t999\t3\tA\tG\n')
        with self.assertRaisesRegex(InputCompatibilityError, 'N differs from preparation override'):
            self.resolve()
        make_bundle(self.root, binary_traits=('B',))
        # A separate reuse manifest may omit an override; it cannot replace the saved N.
        self.manifest = self.root/'reuse.csv'
        self.manifest.write_text('traitname,sample_prevalence,population_prevalence,N\nA,,,\nB,.2,.05,2000\n')
        with self.assertRaisesRegex(InputCompatibilityError, 'manifest N differs'):
            self.resolve()
        self.manifest.write_text(self.manifest.read_text().replace(',2000', ','))
        self.resolve()

    def test_two_bundles_require_the_same_reference(self):
        other = self.root/'other'
        make_bundle(other, names=('B',))
        reference = other/'munged'/REFERENCE_FILE
        reference.write_text(reference.read_text()+'rs3\tA\tC\n')
        seal_bundle(other)
        self.manifest = self.root/'reuse.csv'
        self.manifest.write_text('traitname,munged_file,sample_prevalence,population_prevalence\n'
            f'A,{self.root}/munged/A.sumstats.gz,,\nB,{other}/munged/B.sumstats.gz,,\n')
        args = self.args()[:-2]+['--munged_input']
        with self.assertRaisesRegex(InputCompatibilityError, 'different allele references'):
            genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(args))

    def test_invalid_input_writes_audit_before_any_regression(self):
        self.change_rows('rs1\t1000\t3\tT\tC\n')
        with patch.object(genomicsem_ldsc.subprocess, 'run') as run, contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit): genomicsem_ldsc.main(self.args())
            run.assert_not_called()
        with (self.root/'out/GenomicSEM_Input_QC.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(rows[0]['SNP'], 'rs1')
        self.assertIn('strand-complemented', rows[0]['reason'])

    def test_reference_duplicates_and_invalid_alleles_fail(self):
        reference = self.root/'reference.tsv'
        for text in ('SNP A1 A2\nrs1 A G\nrs1 A G\n', 'SNP A1 A2\nrs1 A A\n', 'SNP A1 A2\n'):
            reference.write_text(text)
            with self.assertRaises(ValueError): read_allele_reference(reference)

    def test_relocated_bundle_and_literal_trait_names(self):
        make_bundle(self.root, names=('NA', 'NaN', '001'))
        destination = self.root/'relocated'
        destination.mkdir()
        for name in ('munged', 'Preparation_Settings.json', 'Prepared_LDSC_Manifest.csv'):
            (self.root/name).rename(destination/name)
        self.manifest = self.root/'reuse.csv'
        self.manifest.write_text('traitname,sample_prevalence,population_prevalence\nNA,,\nNaN,,\n001,,\n')
        args = self.args()[:-2]+['--munged_dir', str(destination/'munged')]
        rows, _ = genomicsem_ldsc.resolve_manifest(genomicsem_ldsc.build_parser().parse_args(args))
        self.assertEqual([r['traitname'] for r in rows], ['NA', 'NaN', '001'])

    def test_malformed_json_and_missing_runtime_provenance_fail(self):
        original = self.sidecar.read_text()
        for text in ('{broken', '[]', 'null'):
            self.sidecar.write_text(text)
            with self.assertRaises(InputCompatibilityError): self.resolve()
        metadata = json.loads(original)
        for field in ('schema_version', 'N_override', 'runtime', 'population_prevalence'):
            altered = {key: value for key, value in metadata.items() if key != field}
            self.sidecar.write_text(json.dumps(altered))
            with self.assertRaises(InputCompatibilityError): self.resolve()
        self.sidecar.write_text(original)

    def test_plain_files_and_literal_ids_with_partial_availability(self):
        reference = {'NA': ('A', 'G'), 'NaN': ('C', 'T'), '001': ('A', 'T')}
        path = self.root/'literal.sumstats'
        path.write_text('SNP\tN\tZ\tA1\tA2\nNA\t100\t0\tA\tG\n001\t100\t-2\tA\tT\n')
        audit = validate_alleles(path, 'NA', reference)
        self.assertEqual(audit['n_snps'], 2)
        self.assertEqual(audit['snp_order_sha256'], hashlib.sha256(b'NA\n001\n').hexdigest())


if __name__ == '__main__': unittest.main()
