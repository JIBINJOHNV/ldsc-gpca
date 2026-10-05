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


def export_workers(n_cores=0):
    """Resolve the shared read/compression budget; 0 selects up to 22 CPUs."""
    if not isinstance(n_cores, int) or isinstance(n_cores, bool) or n_cores < 0:
        raise ValueError('--n_cores must be an integer >= 0 (0 = auto)')
    if n_cores:
        return n_cores
    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count()
    return min(22, available or 1)


def configure_table_threads(argv, *, default=0):
    """Set Polars' pool before preparation imports it; leave CLI errors to its parser.

    A preliminary parse is necessary because the full workflow parser imports the
    preparation module, which imports Polars. Respect an explicit environment
    limit. An already initialized Polars pool in a Python API caller cannot be
    resized; the exporter records its actual size.
    """
    import argparse
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False, exit_on_error=False)
    parser.add_argument('--n_cores', type=int, default=default)
    try:
        n_cores = parser.parse_known_args(argv)[0].n_cores
        workers = export_workers(n_cores)
    except (argparse.ArgumentError, ValueError):
        return
    os.environ.setdefault('POLARS_MAX_THREADS', str(workers))


def configure_gpca_table_threads(argv):
    """Compatibility wrapper for GPCA's automatic worker default."""
    configure_table_threads(argv)
