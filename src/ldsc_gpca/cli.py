"""Dispatch subcommands without interpreting their existing arguments."""
import argparse
import sys
from . import __version__


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="ldsc-gpca", description="Python LDSC and R genomicPCA/GWAMA workflows.",
        epilog="Use ldsc-gpca ldsc --help or ldsc-gpca gpca --help for analysis options.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", choices=["ldsc", "gpca"], nargs="?",
                        help="ldsc: VCF to pairwise LDSC; gpca: Python LDSC results to genomicPCA/GWAMA")
    if argv and argv[0] in ("ldsc", "gpca"):
        if argv[0] == "ldsc":
            from .ldsc_cli import main as run
        else:
            from .gpca import main as run
        return run(argv[1:])
    parser.parse_args(argv)
    parser.print_help()
    return 0
