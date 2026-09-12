"""Run genomicPCA in R, then automatically combine its successful outputs."""
from pathlib import Path
from importlib.resources import as_file, files
import shutil
import subprocess
import sys
from .helptext import HelpParser, print_analysis_help, PREPARE_INPUT_HELP


def postprocess_parser(r_script='gpsca_gwama_python_ldsc.r'):
    from .prepare import add_prepare_options
    prog = 'ldsc-gpca genomicsem gpca' if r_script == 'gpsca_gwama_v2.r' else 'ldsc-gpca gpca'
    parser = HelpParser(prog=prog, add_help=False,
                        description='Preparation and automatic GWAMA postprocessing options. Summary output: <outdir>/harmonisation_input/.',
                        epilog=PREPARE_INPUT_HELP + '\nFor automatic VCF preparation, omit --gpca_input_folder and include vcf_files in the manifest.\nWith an existing GPCA folder, only traitname is required in the manifest.')
    parser.add_argument('--dataset-id', help='Dataset identifier used as the output filename prefix. Default: name of the --outdir folder.')
    parser.add_argument('--gwama-output-n-eff', dest='n_eff', type=float, help='Override N_eff only in the exported GWAMA selected-column summary (e.g. 330000), not LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--gwama-output-info', dest='info_value', type=float, help='Override INFO only in the exported GWAMA selected-column summary (e.g. 0.9), not variant filtering or LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--archive-chromosomes', action='store_true', help='Move current-run source results/logs to outdir/chromosome_wise after saving outputs. Default: keep originals.')
    parser.add_argument('--postprocess-help', action='store_true', help='Show these options without requiring R.')
    add_prepare_options(parser)
    return parser


def main(argv=None, *, r_script='gpsca_gwama_python_ldsc.r'):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = postprocess_parser(r_script)
    def argument_error(message):
        print_analysis_help(r_script, file=sys.stderr)
        parser.print_help(sys.stderr)
        parser.exit(2, f'\nERROR: {message}\n')
    parser.error = argument_error
    if not argv or any(a in ('--help', '-h') for a in argv):
        print_analysis_help(r_script, argv)
        print('\nOPTIONAL VCF PREPARATION AND GWAMA EXPORT\n')
        parser.print_help()
        return 0
    opts, r_args = parser.parse_known_args(argv)
    if opts.postprocess_help:
        parser.print_help()
        return 0
    from .postprocess import validate_overrides, filename_component
    try:
        validate_overrides(opts.n_eff, opts.info_value)
        if opts.dataset_id is not None:
            filename_component(opts.dataset_id)
    except ValueError as error:
        parser.error(str(error))
    rscript = shutil.which('Rscript')
    if rscript is None:
        print_analysis_help(r_script, file=sys.stderr)
        parser.print_help(sys.stderr)
        print('ERROR: Rscript was not found on PATH. Install R and the R packages argparse, data.table and glue.', file=sys.stderr)
        return 127
    previous = None
    out_parser = HelpParser(prog=parser.prog, add_help=False)
    out_parser.error = argument_error
    out_parser.add_argument('--outdir', required=True)
    out_parser.add_argument('--input')
    out_parser.add_argument('--gpca_input_folder')
    out_parser.add_argument('--splitby_chr', choices=['split','nosplit'], default='split')
    out_parser.add_argument('--validate_only', action='store_true')
    inputs = out_parser.parse_known_args(r_args)[0]
    outdir = inputs.outdir
    from .prepare import preparation_kwargs
    prep_settings = preparation_kwargs(opts)
    if inputs.gpca_input_folder is not None and not inputs.validate_only:
        if any(value != parser.get_default(name) for name,value in prep_settings.items()):
            parser.error('VCF preparation settings require omitting --gpca_input_folder; use ldsc-gpca prepare to create new tables')
    if inputs.gpca_input_folder is None and not inputs.validate_only:
        if not inputs.input:
            parser.error('Automatic VCF preparation requires --input with traitname and vcf_files columns')
        from .prepare import prepare_inputs
        from polars.exceptions import PolarsError
        try:
            folder = prepare_inputs(inputs.input, outdir, splitby_chr=inputs.splitby_chr,
                                    **prep_settings)
        except (OSError, ValueError, PolarsError) as error:
            print(f'ERROR: VCF preparation failed: {error}', file=sys.stderr)
            return 1
        r_args += ['--gpca_input_folder', str(folder)]
    harmonised_output = Path(outdir) / 'harmonisation_input'
    from .postprocess import snapshot_outputs
    previous = snapshot_outputs(outdir)
    script = files('ldsc_gpca').joinpath('r', r_script)
    with as_file(script) as path:
        code = subprocess.run([rscript, str(path), *r_args], check=False).returncode
    if code != 0:
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
