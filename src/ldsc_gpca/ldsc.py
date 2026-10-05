"""Select an LDSC backend while preserving its parser and analysis workflow."""
import copy
import sys

from .helptext import HelpParser


def ldsc_parser(backend):
    """Return the existing backend parser, also used by pipeline."""
    if backend == 'python':
        from .ldsc_cli import parser
        return parser
    if backend == 'genomicsem':
        from .genomicsem_ldsc import build_parser
        return build_parser()
    raise ValueError(f'Unknown LDSC backend: {backend}')


def add_backend_option(parser):
    parser.add_argument('--ldsc_backend', choices=['python', 'genomicsem'], default='python',
        help='Munging and regression implementation. Default: python. Reused files skip munging; '
             'their preparation must match the selected backend. --munge_backend belongs to prepare only.')


def build_parser(backend='python'):
    # Do not mutate the Python parser shared by direct callers and pipeline.
    parser = copy.deepcopy(ldsc_parser(backend))
    parser.prog = 'ldsc-gpca ldsc'
    parser.usage = '%(prog)s [--ldsc_backend {python,genomicsem}] --input MANIFEST.csv --outdir DIRECTORY --ld_ref DIRECTORY [options]'
    add_backend_option(parser.add_argument_group('Backend selection'))
    parser.epilog = (f'SELECTED BACKEND: {backend}\n'
        '  Raw inputs use this backend for both munging and regression.\n'
        '  Select --ldsc_backend genomicsem --help for GenomicSEM inputs and options.\n'
        '  Select --ldsc_backend python --help for Python inputs and options.\n'
        '  The existing genomicsem ldsc command remains available.\n\n' + parser.epilog)
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    selector = HelpParser(prog='ldsc-gpca ldsc', add_help=False)
    add_backend_option(selector)
    selected, remaining = selector.parse_known_args(argv)
    parser = build_parser(selected.ldsc_backend)
    if not remaining:
        parser.print_help()
        return 0
    # Validate with the selected parser before entering either analysis workflow.
    # Only the routing option is removed; backend options retain their exact values.
    parser.parse_args(argv)
    if selected.ldsc_backend == 'python':
        from .ldsc_cli import main as run
    else:
        from .genomicsem_ldsc import main as run
    return run(remaining)
