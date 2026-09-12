"""Dispatch subcommands without interpreting their existing arguments."""
import argparse
import sys
from . import __version__


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="ldsc-gpca", description="Python LDSC and R genomicPCA/GWAMA workflows.",
        epilog="Use ldsc-gpca prepare --help, ldsc-gpca ldsc --help or ldsc-gpca gpca --help.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", choices=["prepare", "ldsc", "gpca"], nargs="?",
                        help="prepare: VCF to input tables; ldsc: pairwise LDSC; gpca: genomicPCA/GWAMA")
    if argv and argv[0] in ("prepare", "ldsc", "gpca"):
        if argv[0] == "prepare":
            from .prepare import main as run
        elif argv[0] == "ldsc":
            from .ldsc_cli import main as run
        else:
            from .gpca import main as run
        return run(argv[1:])
    parser.parse_args(argv)
    parser.print_help()
    return 0
