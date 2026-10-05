"""GenomicSEM input for the shared GPCA workflow."""
import sys
from .helptext import HelpParser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = HelpParser(prog='ldsc-gpca genomicsem',
        description='GenomicSEM munging/LDSC and GPCA/GWAMA.',
        epilog='''ldsc runs GenomicSEM munging then LDSC, or reuses munged files.
The same workflow is available as ldsc-gpca ldsc --ldsc_backend genomicsem.
For reuse, prepare files with --mode ldsc (or --mode both) and
--munge_backend genomicsem, using the same --hm3 allele reference for all traits.
Do not directly reuse Python-munged files; reuse does not verify allele orientation.
Use ldsc-gpca genomicsem <command> --help for columns, separators and defaults.''')
    parser.add_argument('command', metavar='COMMAND', choices=['ldsc','gpca'], nargs='?', help='ldsc: munge/LDSC; gpca: PCA/GWAMA from LDSCoutput RData.')
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
