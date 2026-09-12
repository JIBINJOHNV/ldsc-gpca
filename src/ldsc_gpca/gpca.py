"""Run genomicPCA in R, then automatically combine its successful outputs."""
import argparse
from pathlib import Path
from importlib.resources import as_file, files
import shutil
import subprocess
import sys


def postprocess_parser():
    parser = argparse.ArgumentParser(prog='ldsc-gpca gpca', add_help=False, allow_abbrev=False,
                                     description='Automatic post-processing after successful GPCA/GWAMA. The selected-column summary is always saved in <outdir>/harmonised/.')
    parser.add_argument('--dataset-id', help='Dataset identifier used as the output filename prefix. Default: name of the --outdir folder.')
    parser.add_argument('--gwama-output-n-eff', dest='n_eff', type=float, help='Override N_eff only in the exported GWAMA selected-column summary (e.g. 330000), not LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--gwama-output-info', dest='info_value', type=float, help='Override INFO only in the exported GWAMA selected-column summary (e.g. 0.9), not variant filtering or LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--archive-chromosomes', action='store_true', help='Move current-run source results/logs to outdir/chromosome_wise after saving outputs. Default: keep originals.')
    parser.add_argument('--postprocess-help', action='store_true', help='Show these options without requiring R.')
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = postprocess_parser()
    opts, r_args = parser.parse_known_args(argv)
    help_requested = any(a in ('--help', '-h') for a in r_args)
    if opts.postprocess_help or help_requested:
        parser.print_help()
        print('\nR analysis options:')
        if opts.postprocess_help:
            return 0
    if not help_requested:
        from .postprocess import validate_overrides, filename_component
        try:
            validate_overrides(opts.n_eff, opts.info_value)
            if opts.dataset_id is not None:
                filename_component(opts.dataset_id)
        except ValueError as error:
            parser.error(str(error))
    rscript = shutil.which('Rscript')
    if rscript is None:
        print('ERROR: Rscript was not found on PATH. Install R and the R packages argparse, data.table and glue.', file=sys.stderr)
        return 127
    previous = None
    if not help_requested:
        out_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        out_parser.add_argument('--outdir', required=True)
        outdir = out_parser.parse_known_args(r_args)[0].outdir
        harmonised_output = Path(outdir) / 'harmonised'
        from .postprocess import snapshot_outputs
        previous = snapshot_outputs(outdir)
    script = files('ldsc_gpca').joinpath('r/gpsca_gwama_python_ldsc.r')
    with as_file(script) as path:
        code = subprocess.run([rscript, str(path), *r_args], check=False).returncode
    if code != 0 or help_requested:
        return code
    if '--validate_only' in r_args:
        print('Validation-only run: GWAMA post-processing skipped.')
        return 0
    from .postprocess import process_gwama_results
    try:
        process_gwama_results(outdir, harmonised_output, name=opts.dataset_id,
                             n_eff=opts.n_eff, info_value=opts.info_value,
                             archive=opts.archive_chromosomes, previous_files=previous)
    except (OSError, ValueError) as error:
        print(f'ERROR: GWAMA post-processing failed: {error}', file=sys.stderr)
        return 1
    return 0
