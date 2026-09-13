"""Run genomicPCA in R, then automatically combine its successful outputs."""
from pathlib import Path
import argparse
from importlib.resources import as_file, files
import shutil
import subprocess
import sys
from .helptext import HelpParser, print_analysis_help, PREPARE_INPUT_HELP


def postprocess_parser(r_script='gpsca_gwama_python_ldsc.r', *, include_prepare=False, usage=None):
    from .prepare import add_prepare_options
    prog = 'ldsc-gpca genomicsem gpca' if r_script == 'gpsca_gwama_v2.r' else 'ldsc-gpca gpca'
    parser = HelpParser(prog=prog, add_help=False, usage=usage,
                        description='Automatic GWAMA export options. Summary output: <outdir>/harmonisation_input/.',
                        epilog='These options apply after a successful GWAMA run; --validate_only skips export.')
    parser.add_argument('--dataset-id', help='Dataset identifier used as the output filename prefix. Default: name of the --outdir folder.')
    parser.add_argument('--gwama-output-n-eff', dest='n_eff', type=float, help='Override N_eff only in the exported GWAMA selected-column summary (e.g. 330000), not LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--gwama-output-info', dest='info_value', type=float, help='Override INFO only in the exported GWAMA selected-column summary (e.g. 0.9), not variant filtering or LDSC/GPCA calculations. Default: preserve reported values.')
    parser.add_argument('--archive-chromosomes', action='store_true', help='Move current-run source results/logs to outdir/chromosome_wise after saving outputs. Default: keep originals.')
    parser.add_argument('--postprocess-help', action='store_true', help='Show these options without requiring R.')
    if include_prepare:
        add_prepare_options(parser)
    return parser


def preparation_help_parser(r_script):
    from .prepare import add_prepare_options
    prog = 'ldsc-gpca genomicsem gpca' if r_script == 'gpsca_gwama_v2.r' else 'ldsc-gpca gpca'
    parser = HelpParser(prog=prog, add_help=False,
        description='Optional automatic VCF preparation. Used only when --gpca_input_folder is omitted; skipped with --validate_only.',
        epilog=PREPARE_INPUT_HELP + '\nAnalysis uses --input MANIFEST.csv, --outdir DIRECTORY and --splitby_chr {split,nosplit} (default: split).\nSee --help for analysis inputs and options. With an existing GPCA folder these preparation options do not apply.')
    add_prepare_options(parser)
    return parser


def show_gpca_help(r_script, argv, *, section='main', file=None):
    file = sys.stdout if file is None else file
    if section == 'prepare':
        preparation_help_parser(r_script).print_help(file)
    elif section == 'postprocess':
        postprocess_parser(r_script).print_help(file)
    else:
        print_analysis_help(r_script, argv, file=file)
        print('\nGWAMA EXPORT\n', file=file)
        postprocess_parser(r_script, usage=argparse.SUPPRESS).print_help(file)
        print('\nOptional VCF preparation: use --prepare-help for its options and file columns.\n'
              'It runs only when --gpca_input_folder is omitted and --validate_only is not set.', file=file)


def main(argv=None, *, r_script='gpsca_gwama_python_ldsc.r'):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = postprocess_parser(r_script, include_prepare=True)
    def argument_error(message):
        prep_flags = {flag for action in preparation_help_parser(r_script)._actions for flag in action.option_strings}
        section = 'prepare' if any(flag in message for flag in prep_flags) else 'main'
        show_gpca_help(r_script, argv, section=section, file=sys.stderr)
        parser.exit(2, f'\nERROR: {message}\n')
    parser.error = argument_error
    if '--prepare-help' in argv:
        show_gpca_help(r_script, argv, section='prepare')
        return 0
    if '--postprocess-help' in argv:
        show_gpca_help(r_script, argv, section='postprocess')
        return 0
    if not argv or any(a in ('--help', '-h') for a in argv):
        show_gpca_help(r_script, argv)
        return 0
    opts, r_args = parser.parse_known_args(argv)
    from .postprocess import validate_overrides, filename_component
    try:
        validate_overrides(opts.n_eff, opts.info_value)
        if opts.dataset_id is not None:
            filename_component(opts.dataset_id)
    except ValueError as error:
        parser.error(str(error))
    rscript = shutil.which('Rscript')
    if rscript is None:
        show_gpca_help(r_script, argv, file=sys.stderr)
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
