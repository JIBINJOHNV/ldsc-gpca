"""CLI and orchestration for the Python LDSC workflow."""
import os
import math
import json
import hashlib
import sys
from .helptext import HelpParser, LDSC_INPUT_HELP
import pandas as pd
from .utils import optional_prevalence, is_valid_gz
from .extraction import run_vcf_to_table
from .munging import parallel_munge_sumstats
from .pairwise import parallel_ldsc_analysis
from .results import compile_results, check_saved_filters
from .ldsc_runtime import check_runtime

# Define command-line arguments
parser = HelpParser(prog="ldsc-gpca ldsc", description="Pairwise CBIIT Python LDSC using an isolated Conda environment; no Docker.", epilog=LDSC_INPUT_HELP)

parser.add_argument('-output_folder', '--output_folder', metavar='DIRECTORY', help="Output directory; LDSC tables and logs are written below it.", required=True)
parser.add_argument('-ld_ref_snp_file', '--ld_ref_snp_file', metavar='ALLELES.tsv', help="Whitespace-separated HapMap allele table with SNP,A1,A2 headers; passed to LDSC --merge-alleles.", required=True)
parser.add_argument('-input_file', '--input_file', metavar='MANIFEST.csv', help="Comma-separated trait manifest; required headers and prevalence rules below.", required=True)
parser.add_argument('-ld_ref', '--ld_ref', metavar='DIRECTORY', help="Chromosome LD-score reference directory (not an individual file).", required=True)
parser.add_argument('-n_cores', '--n_cores', help="Number of parallel local workers", default=5, type=int)
runtime_group = parser.add_argument_group('Local tools and isolated LDSC environment')
runtime_group.add_argument('--conda-executable', default=os.environ.get('CONDA_EXE', 'conda'), help='Conda executable name/path; default: CONDA_EXE when set, otherwise conda on PATH.')
target = runtime_group.add_mutually_exclusive_group()
target.add_argument('--ldsc-env', help='Child environment name. Default: ldsc-cbiit unless the setup script configured a prefix.')
target.add_argument('--ldsc-env-prefix', help='Child environment directory; overrides LDSC_GPCA_LDSC_PREFIX saved by setup. Default: saved prefix, otherwise use --ldsc-env.')
runtime_group.add_argument('--bcftools', default='bcftools', help='Local bcftools executable/path. Default: bcftools on PATH; not required with --ldsc_only.')

filter_group = parser.add_argument_group('Variant filters')
mhc = filter_group.add_mutually_exclusive_group()
mhc.add_argument('--exclude-mhc', action='store_true',
                 help='Exclude the MHC interval. Default: disabled; MHC variants are kept.')
mhc.add_argument('--no-mhc-exclude', action='store_false', dest='exclude_mhc',
                 help='Keep MHC variants. Default: enabled (no MHC exclusion).')
parser.set_defaults(exclude_mhc=False)
filter_group.add_argument('--mhc_chr', default='6', help='MHC chromosome. Default: 6.')
filter_group.add_argument('--mhc_start', type=int, default=25000000,
                          help='MHC interval start (inclusive). Default: 25000000.')
filter_group.add_argument('--mhc_end', type=int, default=35000000,
                          help='MHC interval end (inclusive). Default: 35000000.')
filter_group.add_argument('--info-min', type=float, default=0.7,
                          help='Minimum FORMAT/SI imputation score. Default: 0.7.')
filter_group.add_argument('--maf-min', type=float, default=0.01,
                          help='Extraction MAF threshold (inclusive). Default: 0.01. Munging also applies --munge-maf-min.')
filter_group.add_argument('--munge-maf-min', type=float, default=0.005,
                          help='Additional LDSC munging MAF threshold (strictly greater than). Default: 0.005.')
filter_group.add_argument('--max-af-difference', type=float, default=0.2,
                          help='Maximum abs(INFO/AF - INFO/EUR). Default: 0.2.')
filter_group.add_argument('--remove-palindrome', action='store_true',
                          help='Remove A/T and C/G variants in the AF interval during extraction. Default: disabled. Standard LDSC munging subsequently removes ALL palindromic SNPs regardless of this option.')
filter_group.add_argument('--paliandromaf_lower', type=float, default=0.45,
                          help='Lower AF bound for palindromic removal. Default: 0.45.')
filter_group.add_argument('--paliandromaf_upper', type=float, default=0.55,
                          help='Upper AF bound for palindromic removal. Default: 0.55.')

# Flags for LDSC-only execution
parser.add_argument('--ldsc_only', action='store_true', help="Reuse existing munged files; filters are NOT reapplied. Recorded filter settings must match, otherwise rerun without this flag. Total-N traits require .prevalence.json; legacy NEF files are accepted with a provenance warning.")
parser.add_argument('-ldsc_input_folder', '--ldsc_input_folder', help="Path to pre-munged sumstats.gz files.", default=None)
parser.add_argument('--ldsc-retries', type=int, default=1,
                    help='Additional attempts per failed pairwise LDSC command. Default: 1 (2 total attempts); 0 disables retries. Successful batches are not repeated. Result-validation failures are not retried.')
def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    output_folder = os.path.abspath(args.output_folder)
    snp_include_file = os.path.abspath(args.ld_ref_snp_file)
    ld_ref_dir = os.path.abspath(args.ld_ref) + os.sep
    n_parallel = args.n_cores
    if n_parallel < 1:
        parser.error('--n_cores must be positive')
    if args.ldsc_retries < 0:
        parser.error('--ldsc-retries must be >= 0')
    if args.mhc_start < 1 or args.mhc_end < args.mhc_start:
        parser.error('--mhc_start must be positive and --mhc_end must be >= --mhc_start')
    if not 0 <= args.info_min <= 1 or not 0 < args.maf_min < 0.5:
        parser.error('--info-min must be in [0,1] and --maf-min must be in (0,0.5)')
    if not 0 <= args.max_af_difference <= 1:
        parser.error('--max-af-difference must be in [0,1]')
    if not 0 <= args.munge_maf_min < 0.5:
        parser.error('--munge-maf-min must be in [0,0.5)')
    if not 0 <= args.paliandromaf_lower <= args.paliandromaf_upper <= 1:
        parser.error('paliandromaf bounds must satisfy 0 <= lower <= upper <= 1')
    filters = {
        'exclude_mhc': args.exclude_mhc, 'mhc_chr': str(args.mhc_chr),
        'mhc_start': args.mhc_start, 'mhc_end': args.mhc_end,
        'info_min': args.info_min, 'maf_min': args.maf_min, 'munge_maf_min': args.munge_maf_min,
        'max_af_difference': args.max_af_difference,
        'remove_palindrome': args.remove_palindrome,
        'pal_lower': args.paliandromaf_lower, 'pal_upper': args.paliandromaf_upper,
    }
    input_df = pd.read_csv(args.input_file)
    for column in ['gwas_name', 'vcf_files', 'ref', 'pop_prevalence', 'sample_prevalence']:
        if column not in input_df:
            parser.error(f'Manifest is missing {column}')
    if input_df['gwas_name'].isna().any() or input_df['gwas_name'].duplicated().any():
        parser.error('gwas_name must be non-empty and unique')
    if any(not str(name).strip() or '/' in str(name) or '\\' in str(name) for name in input_df['gwas_name']):
        parser.error('gwas_name must be a non-empty filename component')
    for column in ['pop_prevalence', 'sample_prevalence']:
        input_df[column] = [optional_prevalence(value, f'{name}: {column}')
                            for name, value in zip(input_df['gwas_name'], input_df[column])]
        input_df[column] = pd.to_numeric(input_df[column])
    no_population = input_df['pop_prevalence'].isna()
    input_df.loc[no_population, 'sample_prevalence'] = float('nan')
    input_df['sample_size_column'] = ['NEF' if missing else 'N_TOTAL' for missing in no_population]
    input_df['sample_prevalence_source'] = ['not_used' if missing else 'provided' for missing in no_population]
    munge_input_folder = os.path.join(output_folder, 'munge_input')
    ldsc_results_dir = os.path.join(output_folder, 'ldsc_results')
    ldsc_input_folder = os.path.abspath(args.ldsc_input_folder) if args.ldsc_input_folder else os.path.join(output_folder, 'ldsc_input')
    num_refs = len(input_df[input_df['ref'] == 'yes'])
    if not num_refs:
        parser.error('At least one trait must have ref=yes')
    runtime = dict(conda=args.conda_executable, environment=args.ldsc_env or 'ldsc-cbiit',
                   prefix=None if args.ldsc_env else (args.ldsc_env_prefix or os.environ.get('LDSC_GPCA_LDSC_PREFIX')))
    if runtime['prefix']:
        runtime['prefix'] = os.path.abspath(runtime['prefix'])
    try:
        check_runtime(**runtime, bcftools=None if args.ldsc_only else args.bcftools)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    os.makedirs(output_folder, exist_ok=True)
    with open(os.path.join(output_folder, 'LDSC_Runtime.json'), 'w') as handle:
        json.dump({**runtime, 'bcftools': args.bcftools, 'backend': 'CBIIT/ldsc'}, handle, indent=2)
    total_comparisons = num_refs * len(input_df)
    active_parallel = min(n_parallel, total_comparisons)
    batch_size = 1 if total_comparisons <= n_parallel else max(5, min(100, math.ceil(len(input_df) / (n_parallel / num_refs))))
    print("-" * 30)
    print("LDSC PIPELINE START")
    print("-" * 30)

    if not args.ldsc_only:
        print("\n[1/4] Converting VCF to TSV...")
        input_df = run_vcf_to_table(n_parallel, input_df, output_folder, filters, args.bcftools)

        print("\n[2/4] Munging Summary Statistics...")
        valid_names = parallel_munge_sumstats(n_parallel, input_df, output_folder, snp_include_file, filters,
                                             ldsc_input_folder=ldsc_input_folder, munge_input_folder=munge_input_folder, runtime=runtime)

        removed_count = len(input_df) - len(valid_names)
        input_df = input_df[input_df['gwas_name'].isin(valid_names)].reset_index(drop=True)

        print(f"\n--- PRE-FLIGHT SUMMARY ---")
        print(f"Traits Removed:  {removed_count}")
        print(f"Traits Retained: {len(input_df)}")
        print(f"Reference Traits: {len(input_df[input_df['ref'] == 'yes'])}")
    else:
        print("\n[!] Reusing munged files: extraction/munging filters are NOT reapplied.")
        for index, row in input_df.iterrows():
            prefix = os.path.join(ldsc_input_folder, row['gwas_name'])
            if not is_valid_gz(prefix + '.sumstats.gz'):
                raise ValueError(f"{row['gwas_name']}: missing or invalid munged file")
            sidecar = prefix + '.prevalence.json'
            if not os.path.exists(sidecar):
                if row['sample_size_column'] == 'N_TOTAL':
                    raise ValueError(f"{row['gwas_name']}: rerun without --ldsc_only to establish total-N provenance")
                check_saved_filters({}, filters, row['gwas_name'])
                continue  # Legacy NEF files remain supported.
            with open(sidecar) as handle:
                metadata = json.load(handle)
            check_saved_filters(metadata, filters, row['gwas_name'])
            with open(prefix + '.sumstats.gz', 'rb') as handle:
                digest = hashlib.sha256()
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            if metadata['sample_size_column'] != row['sample_size_column'] or metadata['sha256'] != digest.hexdigest():
                raise ValueError(f"{row['gwas_name']}: munged file or N convention changed; rerun without --ldsc_only")
            if row['sample_size_column'] == 'N_TOTAL' and pd.isna(row['sample_prevalence']):
                value = optional_prevalence(metadata['sample_prevalence'], row['gwas_name'])
                if value is None:
                    raise ValueError(f"{row['gwas_name']}: saved sample prevalence is missing")
                input_df.loc[index, 'sample_prevalence'] = value
                input_df.loc[index, 'sample_prevalence_source'] = 'saved_metadata'
        os.makedirs(output_folder, exist_ok=True)
        input_df.to_csv(os.path.join(output_folder, 'LDSC_Trait_Prevalence_Metadata.csv'), index=False)

    if not input_df.empty:
        print("\n[3/4] Running LDSC Genetic Correlation...")
        log_files = parallel_ldsc_analysis(active_parallel, batch_size, ldsc_results_dir, ld_ref_dir, input_df, ldsc_input_folder, retries=args.ldsc_retries, runtime=runtime)

        print("\n[4/4] Compiling final results...")
        compile_results(output_folder, log_files, input_df)

        print(f"\nSUCCESS: Results saved to {output_folder}/ldsc_results.csv")
    else:
        raise RuntimeError('No valid traits left to analyze')
