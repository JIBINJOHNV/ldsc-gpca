# Run Python LDSC

[README](../README.md) · [Inputs](INPUTS.md) · [Pipeline](PIPELINE.md) · [Next: GPCA/GWAMA](PYTHON_GPCA.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Files and references you need](#files-and-references-you-need)
- [Basic command](#basic-command)
- [What happens in order](#what-happens-in-order)
- [Change variant filters](#change-variant-filters)
- [Sample size and binary traits](#sample-size-and-binary-traits)
- [Limit large chi-square values](#limit-large-chi-square-values)
- [Reuse munged files or restart](#reuse-munged-files-or-restart)
- [Failed estimates and execution failures](#failed-estimates-and-execution-failures)
- [Results and next step](#results-and-next-step)
- [All options](#all-options)

## When to use this command

`ldsc-gpca ldsc` estimates trait heritabilities, genetic correlations and
intercepts using Python LDSC. Start with supported GWAS VCFs, or explicitly
reuse munged files. It stops after collecting LDSC results; it does not prepare
GWAMA tables or calculate PCA.

## Files and references you need

For VCF input, use the four-trait
[Python manifest](INPUTS.md#vcf-manifest-for-pipeline) saved as
`/data/traits_python.csv`. Required headers are
`traitname,vcf_files,ref,population_prevalence,sample_prevalence`.
Use `ref=yes` for **every trait** if these results will feed GPCA. This produces
all self-pairs and unordered pairs: four traits need ten unique rows.

VCFs need INFO `AF,EUR` and FORMAT `SI,AF,EZ,LP,NEF`; binary extraction also
needs `NC,NCO` when population prevalence is supplied. These requirements
differ from GWAMA preparation. Use absolute VCF paths: this standalone command
resolves relative VCF paths from the working directory.

Supply a chromosome LD-score directory and a HapMap allele table containing
`SNP,A1,A2`. Separate regression weights are optional. The package runs
regression in its configured child LDSC environment; keep the main
`ldsc-gpca` environment active. See [reference files](INPUTS.md#reference-files)
and [runtime configuration](REFERENCE.md#all-python-ldsc-options).

## Basic command

```bash
ldsc-gpca ldsc \
  --input /data/traits_python.csv \
  --outdir /results/python_ldsc \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist
```

The main successful output is `/results/python_ldsc/ldsc_results.csv`. Preserve
the output directory, munged files and provenance sidecars if you want to reuse
or restart the analysis later.

## What happens in order

1. Read traits, reference/target roles and prevalence settings.
2. Extract VCF association statistics, applying the extraction filters below.
3. Munge each trait in the Python LDSC environment, matching HapMap alleles and
   applying native munging filters.
4. If requested, filter high chi-square variants into separate input copies.
5. Run enough LDSC reference batches to cover the requested pairs, retrying
   failed commands according to `--ldsc_retries`.
6. Validate and collect numeric estimates. Add covariance and trait-wide
   normalization columns while writing `ldsc_results.csv`; preserve original
   Python `rg`, SE, Z and P.

With some `ref=no` traits, reference-versus-target output may not contain all
pairs required by GPCA. A single upstream `--rg A,B,C` command calculates A–B
and A–C, not B–C; complete coverage matters.

## Change variant filters

For example, make the extraction INFO threshold stricter, exclude MHC and
supply separate regression weights:

```bash
ldsc-gpca ldsc \
  --input /data/traits_python.csv \
  --outdir /results/python_ldsc_filtered \
  --ld_ref /references/eur_w_ld_chr \
  --ld_weights /references/weights_hm3_no_hla \
  --hm3 /references/w_hm3.snplist \
  --info_min 0.9 \
  --maf_min 0.01 \
  --munge_maf_min 0.005 \
  --max_af_difference 0.2 \
  --exclude_mhc \
  --chisq_max auto \
  --n_cores 8
```

These are analysis choices, not universally preferred thresholds. The default
extraction INFO threshold is `0.7`, using **FORMAT/SI**. Extraction MAF is
inclusive at `0.01`; native munging separately requires MAF strictly greater
than `0.005` by default. Frequency-difference filtering compares **INFO/AF with
INFO/EUR**. Missing values in either field are excluded and counted; FORMAT/AF
is not a fallback for that comparison.

MHC exclusion is off by default. When enabled, it removes the inclusive interval
on chromosome 6 from 25,000,000 through 35,000,000 unless overridden. Choose
coordinates appropriate for your inputs' build.

`--remove_palindrome` optionally removes A/T and C/G records in the requested
AF interval during extraction. Standard LDSC munging subsequently removes all
palindromic SNPs regardless of that extraction setting. These LDSC filters do
not filter the separate GWAMA tables.

## Sample size and binary traits

With population prevalence absent, extraction uses FORMAT/NEF for N and ignores
sample prevalence. With population prevalence supplied, it uses NC+NCO;
missing sample prevalence is inferred from the median per-SNP case fraction.
There is no constant-N CLI override for this extraction route. See
[N and prevalence](INPUTS.md#sample-size-and-prevalence) before choosing the
manifest values. Preserve generated provenance sidecars for reuse.

## Limit large chi-square values

The default is **no additional managed chi-square cutoff**. Use
`--chisq_max 80` for a fixed positive integer cutoff, or `--chisq_max auto` for:

```text
cutoff for each trait = max(80, 0.001 × maximum matched N)
retain SNPs with Z² <= cutoff
```

Automatic N is calculated after matching both reference and weight LD-score
SNPs and dropping incomplete matched rows. For maximum matched N of 50,000,
the cutoff is 80; for 200,000, it is 200. This runs separately for each trait,
including reused munged inputs, and writes filtered copies and audits. Original
munged files stay available. Rows with missing Z are left for native LDSC's
missing-value handling.

This preprocessing is not the native Python LDSC cross-product filter, nor a
guarantee of identical estimates to GenomicSEM. Do not pass fractional values
or zero to this managed option.

## Reuse munged files or restart

Save the [Python reuse manifest](INPUTS.md#reusing-munged-files) as
`/data/traits_python_munged.csv`, and keep files such as
`/data/python_munged/Trait_A.sumstats.gz` together with their sidecars:

```bash
ldsc-gpca ldsc \
  --input /data/traits_python_munged.csv \
  --outdir /results/python_ldsc_reused \
  --ld_ref /references/eur_w_ld_chr \
  --ldsc_only \
  --munged_dir /data/python_munged \
  --chisq_max auto
```

`--ldsc_only` selects reuse. `--munged_dir` alone does not: without the flag it
changes the destination of newly munged files. If omitted with `--ldsc_only`,
the default directory is `<outdir>/ldsc_input`.

Reuse skips extraction and munging, so their filters are not reapplied. Recorded
filter settings and input hashes must match the requested provenance. Total-N
traits require their `.prevalence.json` sidecars; legacy NEF files can be used
with a provenance warning. Previously rounded Z values cannot be recovered by
reusing a munged file. See the [full reuse/restart rules](REFERENCE.md#restarting-interrupted-ldsc-runs).

To continue an interrupted run, repeat the same command and output path with
`--restart`. Only completed batches with matching input contents, settings,
runtime and intact results are reused. Changed or unverified batches are rerun.
`--restart` alone does not skip extraction/munging; also use `--ldsc_only` and
the correct munged directory when you want that behavior.

## Failed estimates and execution failures

Each failed pairwise command or incomplete export gets one additional attempt
by default: **two attempts in total**. `--ldsc_retries 2` allows three total;
`0` disables these retries. Exhausted execution failures stop analysis.
Extraction/munging workers also have a second attempt, separately from this flag.

Numerically invalid estimates are handled by `--result_failure_action`:

- `error` (default): save diagnostics and stop.
- `report`: save diagnostic results; a failed run may have no usable final CSV.
- `drop_traits`: audit removals and write a complete retained subset. Requires
  `ref=yes` for all traits and at least two retained traits.

These policies do not excuse structural errors or conflicting estimates.
Changing the retained trait set changes the downstream PCA. Read
`LDSC_Retained_Traits.csv` and `LDSC_Dropped_Traits.csv` when dropping is enabled.

## Results and next step

`ldsc_results.csv` includes original Python `rg,se,z,p`, heritabilities,
intercepts and their SEs. Additional columns include `gcov_native_obs`,
`h2_p1_pair_obs`, `h2_p2_pair_obs`, `gcov_pair`, `rg_trait_wide` and
`normalization_status`. Native covariance, reconstructed covariance and a
correlation are different quantities; see
[normalization and covariance definitions](REFERENCE.md#optional-trait-wide-correlation-normalization).

Check numeric-QC reports and pair coverage before passing the CSV to
[`ldsc-gpca gpca`](PYTHON_GPCA.md). Four selected traits need ten unique
self/pair combinations. GPCA defaults to original `rg`; merely adding
`rg_trait_wide` does not select it for PCA. See the
[complete output catalog](REFERENCE.md#python-ldsc-outputs).

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Output directory. Use a fresh directory for a separate analysis. |
| `--ld_ref` | Required | Directory of chromosome LD scores and M reference files. |
| `--ld_weights` | Use `--ld_ref` | Optional separate directory of regression-weight LD scores. |
| `--hm3` | Unset | Whitespace SNP/A1/A2 allele reference, required unless reusing munged files. |
| `--ldsc_only` | Off | Reuse Python munged files; skip extraction/munging. Their INFO/MAF and other extraction filters are not reapplied. |
| `--munged_dir` | <outdir>/ldsc_input | With --ldsc_only, read existing {trait}.sumstats.gz files. Otherwise choose the destination for new munging; this option alone does not enable reuse. |

### LDSC filtering and estimation

| Option | Default | What it changes |
| --- | --- | --- |
| `--info_min` | `0.7` | Python extraction: inclusive minimum FORMAT/SI, from 0 to 1. |
| `--maf_min` | `0.01` | Python extraction: retain FORMAT/AF between MAF and 1−MAF, inclusive. Threshold strictly between 0 and 0.5. |
| `--munge_maf_min` | `0.005` | Python munging: require MAF strictly greater than this value; threshold from 0 up to, but not including, 0.5. |
| `--max_af_difference` | `0.2` | Python extraction: maximum absolute INFO/AF minus INFO/EUR difference, from 0 to 1. Missing fields are excluded and counted. |
| `--exclude_mhc` | Off | Python extraction: remove the inclusive interval defined by the three MHC options below. |
| `--mhc_chr` | `6` | Chromosome for MHC exclusion; relevant when `--exclude_mhc` is enabled. |
| `--mhc_start` | `25000000` | Inclusive MHC start, positive integer; check your genome build. |
| `--mhc_end` | `35000000` | Inclusive MHC end, at least the start position. |
| `--remove_palindrome` | Off | Python extraction: remove A/T and C/G variants within the specified AF interval. Native Python munging later removes all palindromic SNPs regardless. |
| `--paliandromaf_lower` | `0.45` | Lower AF bound for extraction palindrome removal; use this exact option spelling. Bounds must satisfy 0 ≤ lower ≤ upper ≤ 1. |
| `--paliandromaf_upper` | `0.55` | Upper AF bound for extraction palindrome removal; use this exact option spelling. |
| `--chisq_max` | Disabled | Positive integer or auto. Auto uses max(80, 0.001 × maximum matched N) per trait; keep Z² ≤ cutoff in separate input copies. Applies to reused inputs too. |

### Failures and reuse

| Option | Default | What it changes |
| --- | --- | --- |
| `--result_failure_action` | `error` | Python LDSC collection: `error`, `report`, `drop_traits`. All save diagnostics. Report may leave no final CSV; drop needs all ref=yes and at least two retained traits. Structural conflicts remain fatal. |
| `--restart` | Off | Standalone Python LDSC only: reuse verified matching completed batches in the same outdir. Does not itself skip extraction/munging. |
| `--ldsc_retries` | `1` | Additional attempts per failed Python LDSC command/incomplete export; integer ≥0. Default means two total attempts. Numerical failures follow QC policies. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--n_cores` | `5` | Positive integer; concurrent extraction, munging, filtering and LDSC workers. |
| `--bcftools` | `bcftools` on PATH | Executable name/path for VCF extraction or preparation; unnecessary when neither runs. |
| `--conda_executable` | `CONDA_EXE`, otherwise `conda` | Executable name/path used to launch the child Python LDSC environment. |
| `--ldsc_env` | Configured prefix, otherwise `ldsc-cbiit` | Explicit child environment name, mutually exclusive with `--ldsc_env_prefix`. The supplied installer configures a prefix, so normally omit both. |
| `--ldsc_env_prefix` | `LDSC_GPCA_LDSC_PREFIX` when set | Explicit child environment directory overrides the saved prefix; mutually exclusive with `--ldsc_env`. Without a prefix, use the named environment. |

### Help

| Option | Default | What it changes |
| --- | --- | --- |
| `--help` | Off | Show command usage, defaults and choices; exit without analysis. |
