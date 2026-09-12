"""Native GenomicSEM input for the shared GPCA workflow."""
import sys
from .helptext import HelpParser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = HelpParser(prog='ldsc-gpca genomicsem',
        description='GenomicSEM backend: GPCA/GWAMA from an existing LDSCoutput RData.',
        epilog='Only gpca is integrated. GenomicSEM munge and ldsc must be run separately.\nUse ldsc-gpca genomicsem gpca --help for the manifest, RData and GWAMA file contracts.')
    parser.add_argument('command', choices=['gpca'], nargs='?', help='Run native GenomicSEM PCA/GWAMA.')
    if argv and argv[0] == 'gpca':
        from .gpca import main as run
        return run(argv[1:], r_script='gpsca_gwama_v2.r')
    if argv:
        parser.parse_args(argv)
    parser.print_help()
    return 0
