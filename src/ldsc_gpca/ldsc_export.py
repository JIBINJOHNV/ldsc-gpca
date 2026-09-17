"""Export native LDSC results in-process, before human-readable formatting.

This file runs as a standalone script in the isolated LDSC environment. It
does not import the main package or replace any regression/statistical code.
"""
import contextlib
import os
import runpy
import shutil
import sys
import tempfile


RESULT_SUFFIX = ".results.csv"
FLOAT_FORMAT = "%.17g"
BASE_RESULT_COLUMNS = (
    "p1", "p2", "rg", "se", "z", "p",
    "h2_int", "h2_int_se", "gcov_int", "gcov_int_se",
)
HERITABILITY_COLUMNS = (
    ("h2_obs", "h2_obs_se"), ("h2_liab", "h2_liab_se"),
)


def has_result_columns(columns):
    columns = set(columns)
    return set(BASE_RESULT_COLUMNS).issubset(columns) and any(
        set(pair).issubset(columns) for pair in HERITABILITY_COLUMNS
    )


def write_results_csv(frame, path):
    """Atomically save the original numerical values, including scientific notation."""
    directory = os.path.dirname(os.path.abspath(path))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=directory,
            prefix=".ldsc-results-", suffix=".tmp", delete=False,
        ) as handle:
            temporary = handle.name
            frame.to_csv(handle, index=False, float_format=FLOAT_FORMAT)
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


@contextlib.contextmanager
def export_rg_tables(sumstats, pandas):
    """Capture the exact DataFrame built by native _get_rg_table, without rebuilding it.

    The to_string hook is active only inside the native table builder and is
    always restored. Each worker is an isolated process; no global dependency
    files or scientific settings are modified.
    """
    original_table = getattr(sumstats, "_get_rg_table", None)
    if not callable(original_table):
        raise RuntimeError("Unsupported LDSC runtime: native _get_rg_table is unavailable.")
    exported = []

    def table_with_export(rg_paths, estimates, args):
        original_to_string = pandas.DataFrame.to_string
        captures = []

        def capture(frame, *positional, **keywords):
            if not has_result_columns(frame.columns):
                raise RuntimeError("Unsupported LDSC result-table schema; no estimates were exported.")
            if captures:
                raise RuntimeError("Native LDSC rendered more than one result table per batch.")
            # Preserve native log output; export this same frame, not its rendered text.
            rendered = original_to_string(frame, *positional, **keywords)
            path = os.fspath(args.out) + RESULT_SUFFIX
            write_results_csv(frame, path)
            captures.append(path)
            return rendered

        pandas.DataFrame.to_string = capture
        try:
            rendered = original_table(rg_paths, estimates, args)
        finally:
            pandas.DataFrame.to_string = original_to_string
        if len(captures) != 1:
            raise RuntimeError("Native LDSC did not expose exactly one numerical result table.")
        exported.extend(captures)
        return rendered

    sumstats._get_rg_table = table_with_export
    try:
        yield exported
    finally:
        sumstats._get_rg_table = original_table


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    native_script = shutil.which("ldsc.py")
    if native_script is None:
        raise RuntimeError("Native ldsc.py is not installed in the isolated LDSC environment.")
    import pandas
    import ldscore.sumstats as sumstats

    previous_argv = sys.argv
    try:
        sys.argv = [native_script, *arguments]
        with export_rg_tables(sumstats, pandas) as exported:
            runpy.run_path(native_script, run_name="__main__")
            if any(arg == "--rg" or arg.startswith("--rg=") for arg in arguments):
                if len(exported) != 1:
                    raise RuntimeError("LDSC finished without exactly one fresh numerical CSV.")
    finally:
        sys.argv = previous_argv
    return 0


if __name__ == "__main__":
    sys.exit(main())
