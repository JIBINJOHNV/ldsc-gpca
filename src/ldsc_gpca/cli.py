"""Dispatch subcommands without interpreting their existing arguments."""
import sys
from . import __version__
from .helptext import HelpParser
from .threads import configure_numerical_threads, configure_table_threads


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Keep numerical imports inside the dispatch branches, after these defaults.
    help_flags = {'--help', '-h', '--version', '--prepare_help', '--postprocess_help'}
    has_workflow_args = len(argv) > (2 if argv[:1] == ['genomicsem'] else 1)
    configure_numerical_threads(report=has_workflow_args and not help_flags.intersection(argv))
    parser = HelpParser(
        prog="ldsc-gpca", description="Python LDSC and R genomicPCA/GWAMA workflows.",
        epilog="""AVAILABLE WORKFLOWS
  pipeline         GWAS VCFs -> Python/GenomicSEM LDSC -> GPCA/GWAMA/export.
  prepare          VCF -> GPCA inputs, shared raw LDSC tables and optional munging.
  ldsc             Munging/LDSC with --ldsc_backend python (default) or genomicsem.
  ldsc.py          Run the raw pinned CBIIT ldsc.py command.
  munge_sumstats.py Run the raw pinned CBIIT munge_sumstats.py command.
  gpca             Python-LDSC results -> PCA/GWAMA.
  genomicsem gpca  Existing GenomicSEM RData -> PCA/GWAMA.
  genomicsem ldsc  GenomicSEM munge -> LDSC, or start from existing munged files.

NOTES
  Numerical libraries default to one thread per worker; existing environment
  settings are preserved. Use --n_cores to choose analysis workers.
  Pipeline accepts VCF study inputs only and prepares all analysis inputs.
  Use standalone ldsc/gpca commands to reuse prepared inputs or completed fits.
  Preparation can stop after munging; ldsc also runs the selected regression.
  ldsc/pipeline use --ldsc_backend for both munging and regression.
  Standalone prepare uses --munge_backend; do not supply it to ldsc/pipeline.
  For GenomicSEM LDSC reuse, prepare files with --mode ldsc (or both) and
  --munge_backend genomicsem. Do not directly reuse Python-munged files.
  Raw .py passthrough commands bypass package validation, filtering and provenance.
  GWAMA postprocessing runs automatically; postprocess is not a top-level command.

Use ldsc-gpca <command> --help for required columns, separators, defaults and choices.
An empty command also displays help; no input files are opened for help.""")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", metavar="COMMAND", choices=["pipeline", "prepare", "ldsc", "ldsc.py", "munge_sumstats.py", "gpca", "genomicsem"], nargs="?",
                        help="Managed workflows or raw pinned CBIIT LDSC script passthrough")
    if argv and argv[0] in ("ldsc.py", "munge_sumstats.py"):
        from .ldsc_runtime import run_ldsc_script
        return run_ldsc_script(argv[0], argv[1:])
    if argv and argv[0] == "genomicsem":
        from .genomicsem import main as run
        return run(argv[1:])
    if argv and argv[0] in ("pipeline", "prepare", "ldsc", "gpca"):
        if argv[0] == "pipeline":
            from .pipeline import main as run
        elif argv[0] == "prepare":
            configure_table_threads(argv[1:], default=4)
            from .prepare import main as run
        elif argv[0] == "ldsc":
            from .ldsc import main as run
        else:
            from .gpca import main as run
        return run(argv[1:])
    if argv:
        parser.parse_args(argv)
    parser.print_help()
    return 0
