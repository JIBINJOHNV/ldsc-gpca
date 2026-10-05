"""Completion and content checks shared by consumers of new preparation bundles."""
import json
from pathlib import Path

from .restart import file_digest


def read_json_object(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f'{path}: metadata must be a JSON object')
    return value


def validate_completed_bundle(path, sidecar, metadata, *, cache=None):
    """Return verified sumstats hash and artifact hashes; never edit the bundle."""
    path, sidecar = Path(path), Path(sidecar)
    cache = {} if cache is None else cache
    root = path.parent.parent
    if root not in cache:
        settings = read_json_object(root/'Preparation_Settings.json')
        if settings.get('munging_status') != 'completed':
            raise ValueError('preparation is not a completed bundle')
        if settings.get('munging_manifest_sha256') != file_digest(root/'Prepared_LDSC_Manifest.csv'):
            raise ValueError('prepared manifest checksum mismatch')
        cache[root] = settings
    settings = cache[root]
    if settings.get('munge_backend') != metadata.get('munge_backend'):
        raise ValueError('preparation backend differs from completion record')
    artifacts = settings.get('munging_artifacts', {})
    if not isinstance(artifacts, dict):
        raise ValueError('invalid completed preparation artifact checksums')
    actual = file_digest(path)
    if metadata.get('sha256') != actual or artifacts.get(path.name) != actual:
        raise ValueError('munged checksum mismatch')
    if artifacts.get(sidecar.name) != file_digest(sidecar):
        raise ValueError('sidecar checksum mismatch')
    return actual, artifacts
