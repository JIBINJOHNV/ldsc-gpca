"""Synthetic provenance for boundary tests; these fixtures do not run GenomicSEM."""
import gzip
import json

from ldsc_gpca.genomicsem_inputs import GENOMICSEM_COMMIT, REFERENCE_FILE
from ldsc_gpca.prepare_ldsc import write_ldsc_manifest
from ldsc_gpca.restart import file_digest, write_checkpoint


def seal_bundle(root):
    for path in (root/'munged').glob('*.prevalence.json'):
        meta = json.loads(path.read_text())
        meta['sha256'] = file_digest(root/'munged'/f"{meta['traitname']}.sumstats.gz")
        meta['reference_sha256'] = file_digest(root/'munged'/REFERENCE_FILE)
        write_checkpoint(path, meta)
    write_checkpoint(root/'Preparation_Settings.json', {
        'munging_status': 'completed', 'munge_backend': 'genomicsem',
        'munging_manifest_sha256': file_digest(root/'Prepared_LDSC_Manifest.csv'),
        'munging_artifacts': {p.name: file_digest(p) for p in (root/'munged').iterdir()}})


def make_bundle(root, names=('A', 'B'), binary_traits=()):
    (root/'munged').mkdir(parents=True, exist_ok=True)
    (root/'munged'/REFERENCE_FILE).write_text('SNP\tA1\tA2\nrs1\tA\tG\nrs2\tC\tT\n')
    metadata = {}
    for index, name in enumerate(names):
        with gzip.open(root/'munged'/f'{name}.sumstats.gz', 'wt') as stream:
            stream.write(f'SNP\tN\tZ\tA1\tA2\nrs1\t1000\t{index+2}\tA\tG\nrs2\t1000\t-2\tC\tT\n')
        binary = name in binary_traits
        metadata[name] = dict(N=1000 if binary else None, sample_prevalence=.2 if binary else None,
                              population_prevalence=.05 if binary else None)
        write_checkpoint(root/'munged'/f'{name}.prevalence.json', {
            'schema_version': 1, 'preparation_method': 'shared', 'munge_backend': 'genomicsem',
            'traitname': name, 'sample_prevalence': metadata[name]['sample_prevalence'],
            'population_prevalence': metadata[name]['population_prevalence'],
            'N_override': metadata[name]['N'], 'N_source': 'manifest_N' if binary else 'FORMAT/NEF',
            'N_convention': 'total' if binary else 'NEF', 'sample_size_column': 'N_TOTAL' if binary else 'NEF',
            'reference_file': REFERENCE_FILE,
            'runtime': {'source_commit': GENOMICSEM_COMMIT, 'package_version': '0.0.5', 'r_version': 'fixture'}})
    write_ldsc_manifest(root/'Prepared_LDSC_Manifest.csv', [(name, None) for name in names], metadata, root, munged=True)
    seal_bundle(root)
