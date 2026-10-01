"""Content-verified completion checkpoints for pairwise LDSC batches."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile

from .ldsc_runtime import fingerprint_runtime
from .results import _read_numerical_csv, _trait_name

CHECKPOINT_SUFFIX = '.checkpoint.json'


def _stamp(path):
    stat = os.stat(path)
    return stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def file_digest(path, *, decompress=False):
    """Hash actual contents; gzip timestamps do not change effective sumstats."""
    before = _stamp(path)
    digest = hashlib.sha256()
    opener = gzip.open if decompress else open
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    if before != _stamp(path):
        raise RuntimeError(f'Input changed while hashing: {path}')
    return digest.hexdigest()


def object_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                    separators=(',', ':')).encode()).hexdigest()


def write_checkpoint(path, value):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                dir=os.path.dirname(os.path.abspath(path)), prefix='.ldsc-checkpoint-',
                suffix='.tmp', delete=False) as handle:
            temporary = handle.name
            json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


class RestartContext:
    """Hash shared inputs once; each checkpoint covers only its own trait batch."""
    def __init__(self, ld_ref_dir, ld_weights_dir, sumstats_paths, runtime, parameters):
        reference_files = set()
        for directory in (ld_ref_dir, ld_weights_dir):
            files = set(Path(directory).glob('*.l2.ldscore*'))
            if not files:
                raise RuntimeError(f'No LD-score files available for restart provenance: {directory}')
            reference_files.update(files)
            reference_files.update(Path(directory).glob('*.l2.M*'))
        reference_files = sorted(str(p.absolute()) for p in reference_files if p.is_file())
        self.stamps = {p: _stamp(p) for p in reference_files}
        references = {p: file_digest(p) for p in reference_files}
        self.shared = {'schema': 1, 'references': object_digest(references),
                       'runtime': fingerprint_runtime(**runtime), 'parameters': parameters}
        self.sumstats = {}
        for path in sorted(set(sumstats_paths)):
            self.stamps[path] = _stamp(path)
            self.sumstats[path] = file_digest(path, decompress=True)

    def request(self, command, sumstats_paths):
        return {**self.shared, 'command': command,
                'sumstats': {p: self.sumstats[p] for p in sumstats_paths}}

    def check_unchanged(self, sumstats_paths):
        paths = set(self.stamps) - set(self.sumstats) | set(sumstats_paths)
        if any(self.stamps[p] != _stamp(p) for p in paths):
            raise RuntimeError('LDSC inputs changed during execution; restart checkpoint was not saved')


def validate_batch(path, reference, targets):
    """Validate structure/coverage, keeping numerical QC failures reportable."""
    frame = _read_numerical_csv(path, allow_failed=True)
    expected = [(_trait_name(reference), _trait_name(target)) for target in targets]
    actual = [(_trait_name(a), _trait_name(b)) for a, b in zip(frame.p1, frame.p2)]
    if len(actual) != len(expected) or sorted(actual) != sorted(expected):
        raise RuntimeError(f'Incomplete or unexpected LDSC batch comparisons in {path}')


def reusable_result(path, request, reference, targets):
    """Return an auditable reason; unverified legacy outputs are recomputed."""
    checkpoint = path + CHECKPOINT_SUFFIX
    try:
        with open(checkpoint) as handle:
            saved = json.load(handle)
        if not isinstance(saved, dict) or saved.get('status') != 'completed':
            return False, 'completion_not_verified'
        if saved.get('request') != request:
            return False, 'inputs_parameters_or_runtime_changed'
        if saved.get('result_sha256') != file_digest(path):
            return False, 'result_changed'
        validate_batch(path, reference, targets)
    except FileNotFoundError:
        return False, 'result_or_checkpoint_missing'
    except (OSError, ValueError, RuntimeError, EOFError):
        return False, 'result_or_checkpoint_invalid'
    return True, 'verified_completed_batch'
