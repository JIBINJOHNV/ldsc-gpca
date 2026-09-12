"""Dispatch subcommands without interpreting their existing arguments."""
import sys
from . import __version__
from .helptext import HelpParser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = HelpParser(
        prog="ldsc-gpca", description="Python LDSC and R genomicPCA/GWAMA workflows.",
        epilog="""AVAILABLE WORKFLOWS
  prepare          VCF -> GPCA inputs; optional LDSC munging input tables.
  ldsc             Run Python LDSC using Docker.
  gpca             Python-LDSC results -> PCA/GWAMA.
  genomicsem gpca  Existing GenomicSEM RData -> PCA/GWAMA.

NOT YET INTEGRATED
  genomicsem munge and genomicsem ldsc are not available commands.
  Run those upstream steps separately, then supply their RData to genomicsem gpca.
  GWAMA postprocessing runs automatically; postprocess is not a top-level command.

Use ldsc-gpca <command> --help for required columns, separators, defaults and choices.
An empty command also displays help; no input files are opened for help.""")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", choices=["prepare", "ldsc", "gpca", "genomicsem"], nargs="?",
                        help="prepare: VCF to input tables; ldsc: Python LDSC; gpca: Python-LDSC GPCA; genomicsem: native GenomicSEM workflow")
    if argv and argv[0] == "genomicsem":
        from .genomicsem import main as run
        return run(argv[1:])
    if argv and argv[0] in ("prepare", "ldsc", "gpca"):
        if argv[0] == "prepare":
            from .prepare import main as run
        elif argv[0] == "ldsc":
            from .ldsc_cli import main as run
        else:
            from .gpca import main as run
        return run(argv[1:])
    if argv:
        parser.parse_args(argv)
    parser.print_help()
    return 0
