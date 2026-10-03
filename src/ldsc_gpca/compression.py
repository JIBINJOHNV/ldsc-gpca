"""Shared gzip streams with optional, bounded pigz compression."""
from contextlib import contextmanager
import gzip
import io
import os
import shutil
import subprocess
import tempfile

from .threads import export_workers


def compression_settings(n_cores=0, gzip_level=1):
    n_cores = export_workers(n_cores)
    if not isinstance(gzip_level, int) or isinstance(gzip_level, bool) or not 1 <= gzip_level <= 9:
        raise ValueError('--gzip_level must be an integer from 1 to 9')
    pigz = shutil.which('pigz')
    return {'backend': 'pigz' if pigz else 'python_gzip', 'executable': pigz,
            'level': gzip_level, 'workers': n_cores if pigz else 1}


def _check_compressor(process, errors):
    code = process.wait()
    if code:
        errors.seek(0)
        detail = errors.read(4096).decode('utf-8', errors='replace').strip()
        raise OSError(f'pigz failed (exit {code}): {detail}')


@contextmanager
def compressed_writer(destination, settings, *, text=False):
    """Close/check the compressor before callers publish their staged output."""
    with open(destination, 'wb') as output:
        process = None
        with tempfile.TemporaryFile() as errors:
            if settings['backend'] == 'python_gzip':
                stream = gzip.GzipFile(filename='', fileobj=output, mode='wb',
                                       compresslevel=settings['level'], mtime=0)
            else:
                # Environment options must not change the format or worker budget.
                env = {k: v for k, v in os.environ.items() if k not in ('GZIP', 'PIGZ')}
                process = subprocess.Popen(
                    [settings['executable'], '-n', '-c', f"-{settings['level']}",
                     '-p', str(settings['workers'])],
                    stdin=subprocess.PIPE, stdout=output, stderr=errors, env=env)
                stream = process.stdin
            try:
                with io.TextIOWrapper(stream, encoding='utf-8', newline='') if text else stream as handle:
                    yield handle
                if process is not None:
                    _check_compressor(process, errors)
            except BrokenPipeError:
                if process is not None:
                    _check_compressor(process, errors)
                raise
            finally:
                if process is not None and process.poll() is None:
                    process.terminate()
                    process.wait()
