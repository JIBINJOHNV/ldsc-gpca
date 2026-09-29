"""Numerical thread defaults for the command-line process and its children."""
import os
import sys


THREAD_ENV_VARS = (
    'OPENBLAS_NUM_THREADS',
    'OMP_NUM_THREADS',
    'MKL_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS',
    'NUMEXPR_NUM_THREADS',
)


def configure_numerical_threads(*, report=False):
    """Set missing defaults before numerical imports; respect caller overrides.

    Child processes inherit these values. This does not change worker counts or
    reconfigure numerical libraries already initialized by a Python API caller.
    """
    settings = {name: os.environ.setdefault(name, '1') for name in THREAD_ENV_VARS}
    if report:
        print('[ldsc-gpca] Numerical thread environment: ' +
              ', '.join(f'{name}={value}' for name, value in settings.items()),
              file=sys.stderr, flush=True)
    return settings
