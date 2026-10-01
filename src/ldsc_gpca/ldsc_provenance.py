"""Fingerprint the isolated LDSC runtime; also executable without this package."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys


def runtime_identity():
    import ldscore

    native = shutil.which('ldsc.py')
    if native is None:
        raise RuntimeError('Native ldsc.py is unavailable for restart provenance')
    root = Path(ldscore.__file__).resolve().parent
    sources = {'ldsc.py': Path(native), 'ldsc_export.py': Path(__file__).with_name('ldsc_export.py')}
    sources.update({str(p.relative_to(root)): p for p in root.rglob('*.py')})
    return {
        'python': sys.version, 'executable': sys.executable,
        'platform': platform.platform(),
        'packages': {name: importlib.metadata.version(name)
                     for name in ('numpy', 'pandas', 'scipy', 'bitarray')},
        'sources': {name: hashlib.sha256(path.read_bytes()).hexdigest()
                    for name, path in sorted(sources.items())},
    }


if __name__ == '__main__':
    print(json.dumps(runtime_identity(), sort_keys=True))
