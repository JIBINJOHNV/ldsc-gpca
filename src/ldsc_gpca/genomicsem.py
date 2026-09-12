"""Native GenomicSEM input for the shared GPCA workflow."""
import sys
from .helptext import HelpParser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = HelpParser(prog='ldsc-gpca genomicsem',
        description='Native GenomicSEM munging/LDSC and GPCA/GWAMA.',
        epilog='ldsc runs munge then LDSC, or accepts existing munged files.\nUse ldsc-gpca genomicsem <command> --help for columns, separators and defaults.')
    parser.add_argument('command', choices=['ldsc','gpca'], nargs='?', help='ldsc: munge/LDSC; gpca: PCA/GWAMA from LDSCoutput RData.')
    if argv and argv[0] == 'ldsc':
        from .genomicsem_ldsc import main as run
        return run(argv[1:])
    if argv and argv[0] == 'gpca':
        from .gpca import main as run
        return run(argv[1:], r_script='gpsca_gwama_v2.r')
    if argv:
        parser.parse_args(argv)
    parser.print_help()
    return 0
