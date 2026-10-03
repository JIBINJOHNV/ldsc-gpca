"""Run native munging, then write its in-memory values without decimal rounding.

Executed in the isolated LDSC environment; no installed upstream files change.
"""
import os
import runpy
import shlex
import shutil
import sys
import tempfile


FLOAT_FORMAT = '%.17g'


def write_munged_sumstats(frame, path, *, keep_maf=False):
    """Preserve native column order/missing values and publish only a complete gzip."""
    if not {'SNP', 'N', 'Z'}.issubset(frame.columns):
        raise RuntimeError('Unsupported native munging output: expected SNP, N and Z')
    columns = [c for c in frame.columns if c in ('SNP', 'N', 'Z', 'A1', 'A2')]
    if keep_maf and 'FRQ' in frame.columns:
        columns.append('FRQ')
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=os.path.dirname(os.path.abspath(path)),
                prefix='.ldsc-munge-', suffix='.tmp', delete=False) as handle:
            temporary = handle.name
        frame.to_csv(temporary, sep='\t', index=False, columns=columns,
                     float_format=FLOAT_FORMAT, compression='gzip')
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    native_script = shutil.which('munge_sumstats.py')
    if native_script is None:
        raise RuntimeError('Native munge_sumstats.py is not installed in the isolated LDSC environment')
    # run_path does not add the native script directory for sibling ldsc imports.
    previous_path = sys.path[:]
    try:
        sys.path.insert(0, os.path.dirname(os.path.realpath(native_script)))
        native = runpy.run_path(native_script)
        args = native['parser'].parse_args(arguments)
        frame = native['munge_sumstats'](args, p=False)
    finally:
        sys.path[:] = previous_path
    write_munged_sumstats(frame, args.out + '.sumstats.gz', keep_maf=args.keep_maf)
    message = f'Managed munging wrote native in-memory values with float_format={FLOAT_FORMAT}.'
    with open(args.out + '.log', 'a') as handle:
        handle.write('\nManaged call: ' + shlex.join(arguments) + '\n' + message + '\n')
    print(message)
    return 0


if __name__ == '__main__':
    sys.exit(main())
