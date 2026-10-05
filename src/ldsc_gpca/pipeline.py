"""Orchestrate existing preparation, LDSC and GPCA commands without changing them."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from . import __version__
from .helptext import HelpParser
from .interfaces import read_manifest
from .ldsc import ldsc_parser


COMMON_GPCA = {
    'gpca_input_folder': (str, None, None, 'Existing nine-column GWAMA input directory. Default: prepare from vcf_files unless --validate_only.'),
    'source_path': (str, None, None, 'Optional already-modified GWAMA source. Default: bundled GWAMA function.'),
    'splitby_chr': (str, 'split', ['split', 'nosplit'], 'GWAMA file layout.'),
    'failed_ldsc_action': (str, 'error', ['error', 'drop_traits'], 'GPCA failed-estimate policy; dropping traits changes PC1.'),
    'h2_z_warn_threshold': (float, 2, None, 'GPCA diagnostic h2/SE warning threshold; 0 disables.'),
    'rg_out_of_range_action': (str, 'warn', ['warn', 'error'], 'GPCA action for finite abs(rg)>1.'),
    'pca_matrix': (str, 'correlation', ['correlation', 'covariance'], 'Matrix used for PCA.'),
    'pc1_orientation': (str, 'tutorial', ['tutorial', 'as_computed'], 'Whole-vector PC1 sign convention.'),
    'negative_eigen_action': (str, 'warn', ['warn', 'error'], 'GPCA action for substantive negative eigenvalues.'),
    'matrix_eigen_tolerance': (float, 1e-8, None, 'GPCA relative negative-eigenvalue tolerance.'),
}
PYTHON_GPCA = {
    'rg_normalization': (str, 'pair', ['pair', 'trait_wide'], 'Correlation source: original Python rg or precomputed rg_trait_wide; CTI is unchanged.'),
    'heritability_scale': (str, 'auto', ['auto', 'observed', 'liability', 'mixed'], 'Python-table h2 scale policy; no conversion.'),
    'duplicate_tolerance': (float, .001, None, 'Absolute tolerance for duplicate LDSC estimates.'),
    'duplicate_z_tolerance': (float, .01, None, 'Absolute tolerance for duplicate LDSC z.'),
    'self_rg_tolerance': (float, .01, None, 'Allowed absolute self-rg difference from 1.'),
    'comparison_epsilon': (float, 1e-12, None, 'Non-negative numerical allowance in GPCA comparisons.'),
    'z_consistency_tolerance': (float, .01, None, 'Relative tolerance for reported z versus rg/SE.'),
    'z_consistency_action': (str, 'warn', ['warn', 'error'], 'GPCA action for inconsistent z.'),
    'ldsc_chunk_size': (int, 250000, None, 'Rows per Python LDSC input read chunk.'),
}


def build_parser(backend='python'):
    from .gpca import postprocess_parser
    from .prepare import add_prepare_options
    parser = HelpParser(prog='ldsc-gpca pipeline', add_help=False,
        parents=[copy.deepcopy(ldsc_parser(backend)), postprocess_parser()],
        usage='%(prog)s --input MANIFEST.csv --outdir DIRECTORY --ld_ref DIRECTORY [options]',
        description='One command: preparation -> selected LDSC backend -> GPCA/GWAMA -> export.',
        epilog='''INPUT AND EXECUTION
  --input: comma-separated CSV, >=2 unique traitname values, sample_prevalence,
    population_prevalence. Relative file paths resolve beside this manifest.
  Python VCF: vcf_files, with the existing Python ldsc VCF schema. ref may be
    omitted; pipeline sets every trait to yes for complete pair/self coverage.
    If ref is supplied, every value must already be yes.
  Python reuse: --ldsc_only --munged_dir; keep required prevalence sidecars.
  GenomicSEM: sumstats_file, or --munged_dir, or --munged_input + munged_file.
    Reused munged files must come from prepare --mode ldsc (or --mode both) with
    --munge_backend genomicsem and the same --hm3 allele reference for all traits.
    Do not directly reuse Python-munged files; HapMap matching can retain strand
    complements. Reuse checks do not verify the backend or allele orientation.
    --vcf_input selects the GenomicSEM INFO-preserving VCF adapter; binary traits
    require an appropriate explicit manifest N and both prevalences.
    Without --vcf_input or those LDSC tables, the older quantitative VCF fallback
    uses preparation's raw tables (N=NEF, no INFO). That fallback refuses binary
    traits. GenomicSEM munging aligns retained alleles to --hm3; no ancestry/build
    harmonisation or additional allele harmonisation of reused files is performed.
  GWAMA: --gpca_input_folder with nine-column TSVs, or vcf_files for preparation.
    Bundled GWAMA requires a justified --gwama_output_info for final export.
  --outdir must not exist. Writes prepare/, ldsc/, gpca/, manifests/ and
    Pipeline_Run_Status.json. Existing module outputs/headers are unchanged.
  --validate_only still RUNS LDSC, then PCA, but skips GWAMA/export. GenomicSEM VCF
    conversion still needs preparation to provide unmunged LDSC inputs.
  --chisq_max: Python accepts a positive integer or auto; GenomicSEM a positive
    number. Keep Z^2 <= cutoff independently per trait for LDSC only.
    Omitted: Python filter disabled; GenomicSEM automatic rule.
    It does not remove variants from GWAMA files or harmonize backend estimates.
  --n_cores is inherited from LDSC (Python 5; GenomicSEM 1) and also passed to
    GPCA/export. GenomicSEM LDSC stays sequential; --prepare_workers defaults to 4.
  Stops on failure; no pipeline resume/overwrite. Existing separate commands
    can reuse completed outputs. Original GWAMA implementation is unchanged.
  Use --ldsc_backend genomicsem --help for GenomicSEM-specific LDSC options.
''')
    execution = parser.add_argument_group('Pipeline execution')
    execution.add_argument('--ldsc_backend', choices=['python', 'genomicsem'], default=backend,
                        help='LDSC implementation; input contracts follow the selected backend.')
    execution.add_argument('--validate_only', action='store_true',
                        help='Run upstream LDSC and GPCA QC/PCA, without GWAMA/export.')
    execution.add_argument('--allow_missing_traits', action='store_true',
                        help='Allow audited GPCA removal of traits absent from LDSC results.')
    execution.add_argument('--prepare_help', action='store_true',
                        help='Show shared preparation options without running an analysis.')
    settings = {**COMMON_GPCA, **(PYTHON_GPCA if backend == 'python' else {})}
    group = parser.add_argument_group('GPCA settings')
    for name, (kind, default, choices, help_text) in settings.items():
        group.add_argument('--' + name, type=kind, default=default, choices=choices, help=help_text)
    prep = HelpParser(add_help=False)
    add_prepare_options(prep)
    group = parser.add_argument_group('VCF preparation settings')
    for action in prep._actions:
        # hm3 and Python bcftools already belong to the inherited LDSC parser.
        if not any(flag in parser._option_string_actions for flag in action.option_strings):
            group._add_action(copy.deepcopy(action))
    # Resume is a standalone LDSC operation; pipeline always requires a fresh outdir.
    if backend == 'python':
        action = parser._option_string_actions.pop('--restart')
        parser._remove_action(action)
        for group in parser._action_groups:
            if action in group._group_actions:
                group._group_actions.remove(action)
                group.title = 'LDSC execution and failed results'
        parser._option_string_actions['--munged_dir'].help = 'Existing munged files; required together with --ldsc_only. Default: unset; extract and munge VCFs.'
    parser._option_string_actions['--outdir'].help = 'Fresh pipeline output directory; must not already exist.'
    parser._option_string_actions['--n_cores'].help = 'LDSC and GPCA/export workers; integer >=1. GenomicSEM LDSC regression remains sequential.'
    parser._option_string_actions['--dataset_id'].help = 'Final export filename prefix. Default: pipeline --outdir directory name.'
    parser._option_string_actions['--hm3'].help = 'Whitespace reference with SNP,A1,A2 headers. Default: unset; required whenever munging or munging-table preparation runs.'
    parser._option_string_actions['--bcftools'].help = 'Local executable/path; required for VCF extraction or preparation. Default: bcftools on PATH.'
    parser._option_string_actions['--write_munge_inputs'].help = 'Also request HapMap-selected raw tables. Default: off; enabled automatically for the older GenomicSEM VCF fallback without --vcf_input.'
    parser._option_string_actions['--p_min'].help = 'VCF preparation P floor; also used by GenomicSEM --vcf_input. Default: 1e-300; range (0,1). Does not change GWAMA Z, but can affect Z reconstructed by GenomicSEM munging.'
    return parser


def optional_actions(parser, opts, *, omit=(), overrides=None):
    """Forward existing parser options using their canonical flags and destinations."""
    overrides = overrides or {}
    result = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction) or action.dest in omit:
            continue
        flags = [flag for flag in action.option_strings if flag.startswith('--')]
        if not flags:
            continue
        value = overrides.get(action.dest, getattr(opts, action.dest, None))
        if value is None or value is False:
            continue
        result.append(flags[0])
        if action.nargs != 0:
            result.append(str(value))
    return result


def write_manifest(path, columns, rows):
    with path.open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def resolve_inputs(opts):
    from .postprocess import filename_component, validate_overrides
    from .utils import optional_prevalence
    manifest = Path(opts.input).resolve()
    columns, rows = read_manifest(manifest)
    required = {'traitname', 'sample_prevalence', 'population_prevalence'}
    if not required.issubset(columns) or len(rows) < 2:
        raise ValueError('Pipeline manifest requires >=2 rows and traitname,sample_prevalence,population_prevalence headers.')
    if any(isinstance(value, float) and not math.isfinite(value) for value in vars(opts).values()):
        raise ValueError('Pipeline numeric options must be finite.')
    names = [row['traitname'] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate traitname in pipeline manifest.')
    for row in rows:
        name = row['traitname']
        filename_component(name)
        if any(c.isspace() for c in name):
            raise ValueError('traitname cannot contain whitespace: ' + repr(name))
        sample = optional_prevalence(row['sample_prevalence'], name + ': sample_prevalence')
        population = optional_prevalence(row['population_prevalence'], name + ': population_prevalence')
        if opts.ldsc_backend == 'genomicsem' and (sample is None) != (population is None):
            raise ValueError(name + ': GenomicSEM requires both prevalences or neither.')
        for column in ('vcf_files', 'sumstats_file', 'munged_file'):
            if row.get(column):
                path = Path(row[column])
                row[column] = str((manifest.parent / path).resolve() if not path.is_absolute() else path.resolve())
    if opts.n_cores < 1 or opts.prepare_workers < 1:
        raise ValueError('Pipeline --n_cores and --prepare_workers must be >=1.')
    for name in ('chisq_max', 'p_min'):
        value = getattr(opts, name)
        if name == 'chisq_max' and opts.ldsc_backend == 'python' and value == 'auto':
            continue
        if value is not None and (not math.isfinite(value) or value <= 0 or (name == 'p_min' and value >= 1)):
            raise ValueError('--' + name + ' is outside its valid positive finite range.')
    validate_overrides(opts.n_eff, opts.info_value)
    if opts.dataset_id is not None:
        filename_component(opts.dataset_id)
    if not opts.validate_only and opts.source_path is None and opts.info_value is None:
        raise ValueError('Bundled GWAMA omits INFO; full pipeline export requires a scientifically justified --gwama_output_info.')
    for name in ('ld_ref', 'ld_weights', 'hm3', 'munged_dir', 'gpca_input_folder', 'source_path'):
        value = getattr(opts, name, None)
        if value:
            path = Path(value).resolve()
            if not path.exists():
                raise ValueError('--' + name + ' does not exist: ' + str(path))
            setattr(opts, name, str(path))
    reuse = bool(getattr(opts, 'ldsc_only', False) or opts.munged_dir or getattr(opts, 'munged_input', False))
    if opts.ldsc_backend == 'python':
        if bool(opts.ldsc_only) != bool(opts.munged_dir):
            raise ValueError('Pipeline Python reuse requires --ldsc_only and --munged_dir together.')
        if 'ref' in columns and any(row['ref'] != 'yes' for row in rows):
            raise ValueError('Pipeline needs complete GPCA pair coverage: every supplied ref must be yes.')
        if 'ref' not in columns:
            columns.append('ref')
        for row in rows:
            row['ref'] = 'yes'
    direct_native_vcf = opts.ldsc_backend == 'genomicsem' and opts.vcf_input
    native_vcf = opts.ldsc_backend == 'genomicsem' and not reuse and not direct_native_vcf and not any(row.get('sumstats_file') for row in rows)
    if native_vcf and any(optional_prevalence(row['population_prevalence'], row['traitname']) is not None for row in rows):
        raise ValueError('The older GenomicSEM VCF fallback cannot infer binary N. Use --vcf_input with an appropriate manifest N and both prevalences, or supply raw/munged GenomicSEM tables.')
    needs_prepare = native_vcf or opts.write_munge_inputs or (not opts.validate_only and not opts.gpca_input_folder)
    needed_columns = []
    if needs_prepare or direct_native_vcf or (opts.ldsc_backend == 'python' and not reuse):
        needed_columns.append('vcf_files')
    if opts.ldsc_backend == 'genomicsem' and not opts.munged_dir and not native_vcf and not direct_native_vcf:
        needed_columns.append('munged_file' if opts.munged_input else 'sumstats_file')
    for column in needed_columns:
        for row in rows:
            if not row.get(column) or not Path(row[column]).is_file():
                raise ValueError(f'{row["traitname"]}: missing {column} file.')
    if direct_native_vcf:
        from .genomicsem_ldsc import resolve_manifest
        resolve_manifest(opts)  # Apply the same sample-size/prevalence contract before any stage runs.
    if (not reuse or opts.write_munge_inputs or native_vcf) and not opts.hm3:
        raise ValueError('--hm3 is required for munging/optional munging-table preparation.')
    for name in ('ld_ref', 'ld_weights', 'munged_dir', 'gpca_input_folder'):
        if getattr(opts, name, None) and not Path(getattr(opts, name)).is_dir():
            raise ValueError('--' + name + ' must be a directory.')
    chromosomes = getattr(opts, 'chromosomes', 22)
    if not 1 <= chromosomes <= 22:
        raise ValueError('--chromosomes must be between 1 and 22.')
    ld_ref = Path(opts.ld_ref)
    weights = Path(opts.ld_weights or opts.ld_ref)
    for chrom in range(1, chromosomes + 1):
        for path in (ld_ref/f'{chrom}.l2.ldscore.gz', ld_ref/f'{chrom}.l2.M_5_50', weights/f'{chrom}.l2.ldscore.gz'):
            if not path.is_file():
                raise ValueError('Missing LD reference/weight file: ' + str(path))
    if opts.hm3:
        with Path(opts.hm3).open() as stream:
            if not {'SNP', 'A1', 'A2'}.issubset(stream.readline().split()):
                raise ValueError('Pipeline --hm3 needs whitespace-separated SNP,A1,A2 headers.')
    if native_vcf:
        print('GenomicSEM VCF route: preparation exports N=NEF and no INFO. Use prior-QC quantitative VCFs; GenomicSEM INFO filtering cannot act on an absent field.', flush=True)
    out = Path(opts.outdir).resolve()
    if out.exists():
        raise ValueError('Pipeline --outdir must not exist; choose a fresh directory (no automatic resume/overwrite).')
    return out, columns, rows, reuse, native_vcf, needs_prepare


def prepare_hm3(source, destination):
    """Reuse preparation's SNP-only selection with a normalized one-column TSV."""
    with Path(source).open() as stream, destination.open('x', newline='') as target:
        columns = stream.readline().split()
        index = columns.index('SNP')
        writer = csv.writer(target, delimiter='\t')
        writer.writerow(['SNP'])
        for line in stream:
            if line.strip():
                fields = line.split()
                if len(fields) != len(columns):
                    raise ValueError('Malformed HapMap row; column count differs from header.')
                writer.writerow([fields[index]])


def retained_manifest(rows, audit=None):
    """Read GenomicSEM audit columns as data, then preserve original manifest order."""
    if audit is None:
        return [{'traitname': row['traitname']} for row in rows]
    with audit.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        if 'traitname' not in (reader.fieldnames or []):
            raise ValueError('LDSC retained-trait audit lacks traitname.')
        selected = [row['traitname'] for row in reader]
    original = {row['traitname'] for row in rows}
    if len(selected) < 2 or len(set(selected)) != len(selected) or not set(selected).issubset(original):
        raise ValueError('LDSC retained-trait audit has too few, duplicate or unknown traits.')
    return [{'traitname': row['traitname']} for row in rows if row['traitname'] in set(selected)]


def save_status(path, status):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(status, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def run_stage(name, arguments, required, status, status_path):
    command = [sys.executable, '-m', 'ldsc_gpca', *arguments]
    entry = {'stage': name, 'command': command, 'status': 'running',
             'started_at': datetime.now(timezone.utc).isoformat()}
    status['stages'].append(entry)
    save_status(status_path, status)
    print(f'[pipeline] {name}', flush=True)
    start = time.perf_counter()
    try:
        result = subprocess.run(command, check=False)
        entry['returncode'] = result.returncode
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, command)
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise ValueError(name + ' returned success without expected outputs: ' + ', '.join(missing))
        entry['status'] = 'completed'
    except BaseException as error:
        entry.update(status='failed', error=str(error))
        raise
    finally:
        entry['elapsed_seconds'] = time.perf_counter() - start
        save_status(status_path, status)


def run_pipeline(opts):
    from .gpca import postprocess_parser
    out, columns, rows, reuse, native_vcf, needs_prepare = resolve_inputs(opts)
    out.mkdir(parents=True)
    manifests = out / 'manifests'
    manifests.mkdir()
    status_path = out / 'Pipeline_Run_Status.json'
    status = {'version': __version__, 'status': 'running', 'backend': opts.ldsc_backend,
              'chisq_max': opts.chisq_max, 'chisq_policy': 'automatic_per_trait' if opts.chisq_max == 'auto' else
              'explicit_per_trait' if opts.chisq_max is not None else
              ('disabled' if opts.ldsc_backend == 'python' else 'native_automatic'),
              'input': str(Path(opts.input).resolve()), 'trait_order': [row['traitname'] for row in rows],
              'settings': vars(opts).copy(), 'stages': []}
    save_status(status_path, status)
    try:
        prepare_manifest = manifests / 'input_traits.csv'
        write_manifest(prepare_manifest, columns, rows)
        folder = Path(opts.gpca_input_folder) if opts.gpca_input_folder else out/'prepare'/'gpca_inputs'
        if needs_prepare:
            args = ['prepare', '--input', str(prepare_manifest), '--outdir', str(out/'prepare'),
                    '--n_cores', str(opts.prepare_workers), '--splitby_chr', opts.splitby_chr,
                    '--p_min', str(opts.p_min), '--gpca_id_source', opts.gpca_id_source, '--bcftools', opts.bcftools]
            if native_vcf or opts.write_munge_inputs:
                hm3 = manifests/'hm3_snps.tsv'
                prepare_hm3(opts.hm3, hm3)
                args += ['--write_munge_inputs', '--hm3', str(hm3), '--munge_id_source', opts.munge_id_source]
            run_stage('prepare', args, [out/'prepare'/'gpca_inputs'], status, status_path)
        if native_vcf:
            if 'sumstats_file' not in columns:
                columns.append('sumstats_file')
            for row in rows:
                row['sumstats_file'] = str(out/'prepare'/'munge_inputs'/f'{row["traitname"]}_munge_inputs.txt')
        ldsc_manifest = manifests/'ldsc_traits.csv'
        write_manifest(ldsc_manifest, columns, rows)
        ldsc_dir, gpca_dir = out/'ldsc', out/'gpca'
        args = ['ldsc'] if opts.ldsc_backend == 'python' else ['genomicsem', 'ldsc']
        args += optional_actions(ldsc_parser(opts.ldsc_backend), opts,
            overrides={'input': ldsc_manifest, 'outdir': ldsc_dir, 'hm3': None if reuse else opts.hm3})
        result = ldsc_dir/('ldsc_results.csv' if opts.ldsc_backend == 'python' else 'genomicsem_LDSC.RData')
        audit = (ldsc_dir/'Selected_Traits.csv' if opts.ldsc_backend == 'genomicsem' else
                 ldsc_dir/'LDSC_Retained_Traits.csv' if opts.result_failure_action == 'drop_traits' else None)
        run_stage('ldsc', args, [result] + ([audit] if audit else []), status, status_path)
        selected = retained_manifest(rows, audit)
        selected_manifest = manifests/'gpca_traits.csv'
        write_manifest(selected_manifest, ['traitname'], selected)
        status['gpca_input_trait_order'] = [row['traitname'] for row in selected]
        status['removed_before_gpca'] = [row['traitname'] for row in rows if row['traitname'] not in status['gpca_input_trait_order']]
        args = ['gpca'] if opts.ldsc_backend == 'python' else ['genomicsem', 'gpca']
        args += ['--input', str(selected_manifest), '--ldsc_results', str(result), '--outdir', str(gpca_dir),
                 '--n_cores', str(opts.n_cores)]
        for name in (*COMMON_GPCA, *(PYTHON_GPCA if opts.ldsc_backend == 'python' else {})):
            if name == 'gpca_input_folder':
                continue
            value = getattr(opts, name)
            if value is not None:
                args += ['--' + name, str(value)]
        if opts.validate_only:
            args.append('--validate_only')
        else:
            args += ['--gpca_input_folder', str(folder)]
        if opts.allow_missing_traits:
            args.append('--allow_missing_traits')
        dataset = opts.dataset_id or out.name
        args += optional_actions(postprocess_parser(), opts, omit={'postprocess_help'}, overrides={'dataset_id': dataset})
        required = [gpca_dir/'GenomicPCA_PC1_Weights_Used.csv']
        if not opts.validate_only:
            required += [gpca_dir/'GWAMA_Run_Status.csv', gpca_dir/f'{dataset}_postprocess.json',
                         gpca_dir/f'{dataset}_GWAMA_combined_results.txt.gz',
                         gpca_dir/'harmonisation_input'/f'{dataset}_GPCA_inputs.txt.gz']
        run_stage('gpca', args, required, status, status_path)
        status['status'] = 'completed'
        print(f'[pipeline] Completed. Results: {gpca_dir}', flush=True)
    except BaseException as error:
        status.update(status='failed', error=str(error))
        raise
    finally:
        save_status(status_path, status)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    selector = HelpParser(add_help=False)
    selector.add_argument('--ldsc_backend', choices=['python', 'genomicsem'], default='python')
    backend = selector.parse_known_args(argv)[0].ldsc_backend
    if '--prepare_help' in argv or '--postprocess_help' in argv:
        from .gpca import show_gpca_help
        show_gpca_help('gpsca_gwama_v2.r' if backend == 'genomicsem' else 'gpsca_gwama_python_ldsc.r', argv,
                      section='prepare' if '--prepare_help' in argv else 'postprocess')
        return 0
    parser = build_parser(backend)
    if not argv:
        parser.print_help()
        return 0
    opts = parser.parse_args(argv)
    try:
        return run_pipeline(opts)
    except KeyboardInterrupt:
        print('Pipeline interrupted; completed outputs and status are retained.', file=sys.stderr)
        return 130
    except subprocess.CalledProcessError as error:
        print('Pipeline stopped after a failed stage; see Pipeline_Run_Status.json.', file=sys.stderr)
        return error.returncode if error.returncode > 0 else 1
    except (OSError, ValueError) as error:
        print('Pipeline error: ' + str(error), file=sys.stderr)
        return 1
