"""Run genomicPCA in R, then automatically combine its successful outputs."""
import argparse
from pathlib import Path
from importlib.resources import as_file, files
import shutil
import subprocess
import sys


def postprocess_parser():
    from .prepare import add_prepare_options
    parser = argparse.ArgumentParser(prog='ldsc-gpca gpca', add_help=False, allow_abbrev=False,
                                     description='Automatic post-processing after successful GPCA/GWAMA. The selected-column summary is always saved in <outdir>/harmonisation_input/.')
    parser.add_argument('--dataset-id', help='Dataset identifier used as the output filename prefix. Default: name of the --outdir folder.')
    parser.add_argument('--gwama-output-n-eff', dest='n_eff', type=float, help='Override N_eff only in the exported GWAMA selected-column summary (e.g. 330000), not LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--gwama-output-info', dest='info_value', type=float, help='Override INFO only in the exported GWAMA selected-column summary (e.g. 0.9), not variant filtering or LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--archive-chromosomes', action='store_true', help='Move current-run source results/logs to outdir/chromosome_wise after saving outputs. Default: keep originals.')
    parser.add_argument('--postprocess-help', action='store_true', help='Show these options without requiring R.')
    add_prepare_options(parser)
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
        out_parser.add_argument('--input')
        out_parser.add_argument('--gpca_input_folder')
        out_parser.add_argument('--splitby_chr', choices=['split','nosplit'], default='split')
        out_parser.add_argument('--validate_only', action='store_true')
        inputs = out_parser.parse_known_args(r_args)[0]
        outdir = inputs.outdir
        if inputs.gpca_input_folder is not None and not inputs.validate_only:
            if (opts.write_munge_inputs or opts.hapmap_file or opts.gpca_id_source != 'chr_pos_ref_alt'
                    or opts.munge_id_source != 'vcf_id' or opts.prepare_workers != 4 or opts.bcftools != 'bcftools'):
                parser.error('VCF preparation settings require omitting --gpca_input_folder; use ldsc-gpca prepare to create new tables')
        if inputs.gpca_input_folder is None and not inputs.validate_only:
            if not inputs.input:
                parser.error('Automatic VCF preparation requires --input with traitname and vcf_files columns')
            from .prepare import prepare_inputs, preparation_kwargs
            from polars.exceptions import PolarsError
            try:
                folder = prepare_inputs(inputs.input, outdir, splitby_chr=inputs.splitby_chr,
                                        **preparation_kwargs(opts))
            except (OSError, ValueError, PolarsError) as error:
                print(f'ERROR: VCF preparation failed: {error}', file=sys.stderr)
                return 1
            r_args += ['--gpca_input_folder', str(folder)]
        harmonised_output = Path(outdir) / 'harmonisation_input'
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
