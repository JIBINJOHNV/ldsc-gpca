# Run Python LDSC

[README](../README.md) · [Input formats](INPUTS.md) · [Output file reference](PYTHON_LDSC_OUTPUTS.md) · [Next: GPCA/GWAMA](PYTHON_GPCA.md)

`ldsc-gpca ldsc` defaults to **Python LDSC**; `--ldsc_backend python` selects it
explicitly. It estimates heritabilities, genetic correlations and LDSC
intercepts, then collects them in **`<outdir>/ldsc_results.csv`**. Choose the
instructions that match the files you already have:

| Your starting files | Input mode | Work performed |
| --- | --- | --- |
| One supported GWAS summary-statistics VCF per trait | [Start from VCF files](#start-from-vcf-files); omit `--ldsc_only` | Extract and filter VCF records → munge → optionally filter Z² → run LDSC → collect results. |
| One Python LDSC `{traitname}.sumstats.gz` file per trait | [Start from munged files](#start-from-munged-files); include `--ldsc_only` | Check existing inputs/provenance → optionally filter Z² → run LDSC → collect results. |

Both routes produce the same final result format. This command ends at LDSC;
use [`gpca`](PYTHON_GPCA.md) for PCA/GWAMA, or [`pipeline`](PIPELINE.md) to run
the stages together. Munged LDSC files are not GWAMA input tables.

For **GenomicSEM LDSC**, use [`ldsc-gpca ldsc --ldsc_backend genomicsem`](GENOMICSEM_LDSC.md),
which writes `genomicsem_LDSC.RData`. Reused munged files for that command must
come from `prepare --mode ldsc` (or `--mode both`) with
`--munge_backend genomicsem`, using the same `--hm3` allele reference for every
trait. Do not pass Python-munged files directly to GenomicSEM LDSC: HapMap
matching can retain strand-complement coding without normalizing its orientation.
See the [preparation and reuse example](GENOMICSEM_LDSC.md#reuse-munged-files).

`--ldsc_backend` selects the matching munging and regression implementation.
Do not supply `--munge_backend` to `ldsc`; it belongs to standalone `prepare`.
Use `ldsc-gpca ldsc --ldsc_backend python --help` or
`ldsc-gpca ldsc --ldsc_backend genomicsem --help` for the selected backend's
options. Backend-specific flags are not interchangeable.

## In this guide

- [Start from VCF files](#start-from-vcf-files)
- [Start from munged files](#start-from-munged-files)
- [Which options apply to each input mode?](#which-options-apply-to-each-input-mode)
- [Limit large chi-square values](#limit-large-chi-square-values)
- [Failed estimates and execution failures](#failed-estimates-and-execution-failures)
- [Results and next step: final columns and a real example](#results-and-next-step)
- [All options and defaults](#all-options)

## Start from VCF files

This route runs Python munging and Python LDSC. To run GenomicSEM LDSC from
VCFs, follow the [GenomicSEM VCF guide](GENOMICSEM_LDSC.md#start-from-gwas-vcfs).

### 1. Prepare the manifest and references

Use GWAS **summary-statistics VCFs**, not individual-level genotype VCFs. Save
this example as `/data/traits_python.csv`, replacing the trait names and paths
with your own. The command examples use illustrative paths; data and references
are not bundled.

```csv
traitname,vcf_files,ref,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,yes,,
Trait_B,/data/Trait_B.vcf.gz,yes,,
Trait_C,/data/Trait_C.vcf.gz,yes,,
Trait_D,/data/Trait_D.vcf.gz,yes,,
```

| Manifest column | What to supply |
| --- | --- |
| `traitname` | A unique, non-empty trait name suitable for a filename. This determines output filenames and matrix ordering. |
| `vcf_files` | Path to that trait's VCF. Use absolute paths: standalone `ldsc` resolves relative paths from the working directory. |
| `ref` | Use `yes` for **every trait** when preparing results for GPCA. LDSC runs each reference against all listed traits. |
| `population_prevalence` | Leave empty for a quantitative trait. For a binary trait, supply the population prevalence when requesting liability conversion. |
| `sample_prevalence` | For a binary trait with population prevalence, supply its case fraction or leave empty to infer it from NC and NCO. Ignored without population prevalence. |

The two prevalence headers are required even when their cells are empty.
With all four traits marked `ref=yes`, a complete managed run produces 16
directed rows, including self-pairs: these represent 10 unique unordered/self
combinations. Some `ref=no` traits can leave pairs missing for downstream GPCA.

Each VCF needs INFO `AF,EUR` and FORMAT `SI,AF,EZ,LP,NEF`, along with variant
IDs, coordinates and alleles. Binary extraction with population prevalence also
needs FORMAT `NC,NCO`. The two AF fields have different uses; see
[filter settings](#change-variant-filters).

Supply an ancestry-appropriate chromosome LD-score directory for `--ld_ref`
and a whitespace-separated HapMap allele table with `SNP,A1,A2` headers for
`--hm3`. Separate regression weights can be supplied with `--ld_weights`;
otherwise the LD-reference directory supplies them too. Check builds and allele
conventions against the [reference requirements](INPUTS.md#reference-files).
Keep the main `ldsc-gpca` environment active; it launches the configured child
Python LDSC environment automatically.

### 2. Run with the default filters

```bash
ldsc-gpca ldsc \
  --input /data/traits_python.csv \
  --outdir /results/python_ldsc_vcf \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist
```

### Change variant filters

For example, use a stricter INFO threshold, exclude the MHC interval, supply
separate weights and enable automatic chi-square filtering:

```bash
ldsc-gpca ldsc \
  --input /data/traits_python.csv \
  --outdir /results/python_ldsc_vcf_filtered \
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

These are optional analysis choices. By default:

- Extraction keeps **FORMAT/SI ≥ 0.7** and FORMAT/AF between `0.01` and `0.99`,
  inclusive. `--info_min` and `--maf_min` control these filters.
- Frequency comparison keeps records with an absolute **INFO/AF − INFO/EUR**
  difference ≤ `0.2`. Missing either value excludes the record; FORMAT/AF is
  not a substitute. `--max_af_difference` sets this limit.
- Python LDSC munging applies its own QC and requires MAF **strictly greater than
  `0.005`**, controlled by `--munge_maf_min`.
- MHC exclusion is off. `--exclude_mhc` removes the inclusive chromosome 6
  interval 25,000,000–35,000,000 by default; choose coordinates for your build.
- `--remove_palindrome` optionally removes A/T and C/G records in the chosen
  AF interval during extraction. Python LDSC munging later removes all
  palindromic SNPs regardless of that extraction flag.

These filters apply to LDSC inputs. GWAMA preparation has its own options.

### Sample size and binary traits

| Manifest setting | N passed to munging | Prevalence handling |
| --- | --- | --- |
| No population prevalence | FORMAT/NEF | Sample prevalence is ignored; no liability conversion. |
| Population prevalence supplied | FORMAT/NC + FORMAT/NCO | Use the supplied sample prevalence, or infer the median per-SNP case fraction. |

There is no constant-N override in this command. The generated `.prevalence.json`
files record the N convention for later reuse. See
[sample size and prevalence](INPUTS.md#sample-size-and-prevalence) before choosing
these settings.

### 3. Understand the processing steps and files

1. Read the manifest and verify the required tools and LDSC runtime.
2. Extract each VCF into `munge_input/{trait}_mungeinput.tsv`, applying the
   extraction filters and recording missing allele-frequency counts.
3. Munge each extracted table: match HapMap alleles, apply Python LDSC QC and write
   `ldsc_input/{trait}.sumstats.gz`, its log and provenance sidecar.
4. If `--chisq_max` is supplied, create separately filtered munged copies and
   exclusion reports. Original munged files remain available.
5. Run LDSC reference/target batches, saving numerical results and readable logs.
6. Check estimates and collect `ldsc_results.csv`. Covariance and trait-wide
   correlation annotations are added during collection without extra regressions.

A successful run with default folder settings has this layout. `{trait}` and
`{batch}` stand for generated names, not literal filenames:

```text
/results/python_ldsc_vcf/
├── ldsc_results.csv                       ← final estimates for GPCA
├── ldsc_results_diagnostic.csv            ← all collected estimates before trait removal
├── LDSC_Compilation_Status.csv            ← whether collection produced usable results
├── LDSC_Pair_Status.csv                   ← per-pair QC and source locations
├── LDSC_Trait_Status.csv                  ← self-pair QC and failed-pair counts
├── LDSC_Trait_Prevalence_Metadata.csv     ← resolved trait/N/prevalence settings
├── LDSC_Runtime.json                      ← tools, runtime and run settings
├── LDSC_Extraction_Worker_Attempts.csv    ← extraction successes/retries
├── LDSC_Munging_Worker_Attempts.csv       ← munging successes/retries
├── munge_input/
│   ├── {trait}_mungeinput.tsv             ← extracted association statistics
│   └── {trait}_AF_Filter_QC.csv            ← missing INFO/AF and INFO/EUR counts
├── ldsc_input/
│   ├── {trait}.sumstats.gz                ← munged SNP/N/Z/allele table
│   ├── {trait}.log                        ← Python LDSC munging log
│   └── {trait}.prevalence.json            ← file hash, filters and N convention
└── ldsc_results/
    ├── {batch}.results.csv                ← numerical LDSC batch estimates
    ├── {batch}.log                        ← readable LDSC regression log
    ├── {batch}.results.csv.checkpoint.json ← restart verification record
    └── LDSC_Batch_Status.csv              ← batch completion and attempts
```

`--munged_dir /other/path` moves the **new munged outputs** there; it does not
select reuse. Optional chi-square and trait-removal files are listed
[below](#additional-files-created-by-optional-settings). For every intermediate
and audit schema, use the [output file reference](PYTHON_LDSC_OUTPUTS.md).

## Start from munged files

### 1. Prepare the existing files and manifest

Use one Python LDSC file named **`{traitname}.sumstats.gz`** per trait. These
are gzip-compressed, whitespace-delimited summary statistics containing
`SNP,N,Z,A1,A2`; the managed munging output uses tabs and that column order.
They must already have been munged with suitable SNP/allele references.
A raw GWAS table, a VCF, or a GWAMA input table is not interchangeable with this
input. Preserve each file's `.prevalence.json` sidecar when available.

For example, the input directory contains:

```text
/data/python_munged/
├── Trait_A.sumstats.gz
├── Trait_A.prevalence.json
├── Trait_B.sumstats.gz
├── Trait_B.prevalence.json
├── Trait_C.sumstats.gz
├── Trait_C.prevalence.json
├── Trait_D.sumstats.gz
└── Trait_D.prevalence.json
```

Save `/data/traits_python_munged.csv`. `vcf_files` is not required in this mode:

```csv
traitname,ref,population_prevalence,sample_prevalence
Trait_A,yes,,
Trait_B,yes,,
Trait_C,yes,,
Trait_D,yes,,
```

Names must match the munged filenames exactly. Use `ref=yes` for every trait
for complete GPCA coverage. Prevalence values must agree with how the existing
N values were produced: adding population prevalence does not convert an
existing NEF column into total N. Total-N traits require their sidecars; legacy
NEF inputs without sidecars are accepted with a provenance warning.

Supply the same LD-reference inputs required for regression in the VCF route.
**`--hm3` and `bcftools` are not needed** because extraction and munging are
skipped. The configured child LDSC runtime is still required.

### 2. Run LDSC on the existing inputs

For inputs produced with the default filter settings:

```bash
ldsc-gpca ldsc \
  --input /data/traits_python_munged.csv \
  --outdir /results/python_ldsc_reused \
  --ld_ref /references/eur_w_ld_chr \
  --ldsc_only \
  --munged_dir /data/python_munged
```

`--ldsc_only` explicitly selects this route. If `--munged_dir` is omitted,
files are read from `<outdir>/ldsc_input`. The program does not detect this
mode from filenames or from `--munged_dir` alone.

### Reuse munged files or restart

You can create these files with [`prepare --mode ldsc`](PREPARE.md#produce-munged-files).
Use its `Prepared_LDSC_Manifest.csv` with `--ldsc_only --munged_dir DIR/munged`,
and keep the provenance sidecars. Shared preparation records its own munging
filters; Python VCF extraction filters are not applied to those files.

Reuse checks saved file hashes, the N convention and recorded filter settings.
**Extraction and munging filters are not run again.** When a sidecar records
non-default settings, repeat those settings so the provenance check matches.
For example, to reuse files originally made with `--info_min 0.9 --exclude_mhc`
and otherwise default extraction/munging settings:

```bash
ldsc-gpca ldsc \
  --input /data/traits_python_munged.csv \
  --outdir /results/python_ldsc_reused_filtered \
  --ld_ref /references/eur_w_ld_chr \
  --ldsc_only \
  --munged_dir /data/python_munged_filtered \
  --info_min 0.9 \
  --exclude_mhc \
  --chisq_max auto
```

Here, `--info_min 0.9 --exclude_mhc` describe the existing inputs;
`--chisq_max auto` **does perform a new filtering step** on separate copies.
Match every recorded extraction/munging setting, including stored MHC and
palindrome bounds. To change those filters, return to the original VCFs and
munge again. Missing legacy filter metadata produces a warning rather than
proof that filters matched. Reuse cannot restore Z precision lost in an older
munged file.

To resume an interrupted run, repeat its command and output directory with
`--restart`. Only completed batches whose input contents, references, settings,
runtime and result hashes still match are reused. Changed or unverified batches
are rerun. `--restart` alone does not skip VCF extraction/munging; include
`--ldsc_only` and the correct munged directory to skip those steps. See
[full restart rules](REFERENCE.md#restarting-interrupted-ldsc-runs).

### 3. Understand the processing steps and files

1. Validate the manifest, runtime, munged files and available provenance.
2. If requested, apply `--chisq_max` to separate copies of the existing files.
3. Run LDSC batches, or reuse verified completed batches with `--restart`.
4. Check estimates and collect the same final columns as the VCF route.

A fresh reuse output directory has this layout after a successful run:

```text
/results/python_ldsc_reused/
├── ldsc_results.csv                       ← final estimates for GPCA
├── ldsc_results_diagnostic.csv            ← all collected estimates before trait removal
├── LDSC_Compilation_Status.csv            ← collection outcome
├── LDSC_Pair_Status.csv                   ← per-pair QC and source locations
├── LDSC_Trait_Status.csv                  ← self-pair QC and failed-pair counts
├── LDSC_Trait_Prevalence_Metadata.csv     ← resolved trait/N/prevalence settings
├── LDSC_Runtime.json                      ← tools, runtime and run settings
└── ldsc_results/
    ├── {batch}.results.csv                ← numerical LDSC batch estimates
    ├── {batch}.log                        ← readable regression log
    ├── {batch}.results.csv.checkpoint.json ← restart verification record
    └── LDSC_Batch_Status.csv              ← batch completion and attempts
```

The original munged inputs stay in `--munged_dir`; this route does not create
new extraction tables, munging logs, sidecars or extraction/munging attempt
reports. Files left from an earlier VCF run can still be present if you reuse
its output directory. Optional outputs are described next.

## Which options apply to each input mode?

| Option or setting | Starting from VCFs | With `--ldsc_only` |
| --- | --- | --- |
| Manifest `vcf_files` | Required; read during extraction. | Not required; VCFs are not read. |
| `--hm3` | Required for munging. | Unused; existing files are already munged. |
| `--bcftools` | Used for extraction. | Not required or used for extraction. |
| `--munged_dir` | **Output directory** for new munged files. | **Input directory** holding existing munged files. |
| `--info_min`, `--maf_min`, `--max_af_difference`, `--exclude_mhc`, `--mhc_*`, `--remove_palindrome`, `--paliandromaf_*` | Apply during VCF extraction. | Do not filter again; must match recorded settings when sidecar filter metadata exists. |
| `--munge_maf_min` | Apply during Python LDSC munging. | Do not munge again; must match recorded settings when available. |
| Manifest prevalence | Controls extracted N and requested liability conversion. | Must match the existing N convention; also controls requested liability conversion. |
| `--chisq_max` | Optional filtering after munging. | Optional filtering of existing munged inputs. |
| `--ld_ref`, `--ld_weights`, runtime, workers, retries, result policy | Used for LDSC. | Used for LDSC. |
| `--restart` | Can reuse verified LDSC batches; extraction/munging still run. | Can reuse verified LDSC batches; extraction/munging stay skipped. |

## Limit large chi-square values

The default is **no additional managed chi-square cutoff**. Use
`--chisq_max 80` for a fixed positive integer cutoff, or `--chisq_max auto` for:

```text
cutoff for each trait = max(80, 0.001 × maximum matched N)
retain SNPs with Z² <= cutoff
```

Automatic N is calculated after matching both reference and weight LD-score
SNPs and dropping incomplete matched rows. A maximum matched N of 50,000 gives
a cutoff of 80; 200,000 gives 200. Filtering runs separately for each trait,
including reused munged inputs. Missing Z rows remain for Python LDSC's
missing-value handling; they are counted separately in the audit.

This step writes filtered copies and preserves the source munged files. It is
separate from Python LDSC's cross-product filtering and does not make
Python and GenomicSEM estimates identical. Fractions and zero are not accepted
as explicit managed cutoffs.

### Additional files created by optional settings

All paths below are relative to `--outdir`, for either input mode.

| Setting / event | Additional path | What it represents |
| --- | --- | --- |
| `--chisq_max` supplied | `ldsc_input_chisq_filtered/{trait}.sumstats.gz` | Munged copy actually passed to LDSC after the managed Z² filter. |
| `--chisq_max` supplied | `ldsc_input_chisq_filtered/{trait}.chisq_excluded.tsv.gz` | That trait's removed SNPs, Z values and Z² values. |
| `--chisq_max` supplied | `LDSC_ChiSquare_Filter_Summary.csv` | Per-trait cutoff, source paths and before/after/missing-Z counts. |
| `--chisq_max` supplied | `LDSC_ChiSquare_Excluded_Variants.tsv.gz` | Combined removed-SNP table for all traits. |
| `--chisq_max` supplied | `LDSC_ChiSquare_Worker_Attempts.csv` | Filter-worker successes and retries. |
| `--result_failure_action drop_traits` | `LDSC_Retained_Traits.csv` | Manifest of traits kept for downstream analysis, in original order. |
| `--result_failure_action drop_traits` | `LDSC_Dropped_Traits.csv` | Removal order and reasons; can be empty apart from its header. |
| A subprocess fails | `execution_errors.log` | Captured command errors, when produced; also inspect batch/worker reports. |
| Collection files already exist | `*.previous-<uuid>` | Backups of selected earlier results/audits, not the current run's results. |

The [output reference](PYTHON_LDSC_OUTPUTS.md) lists the columns and explains
how to read these files. A failure can stop the run before later files exist.

## Failed estimates and execution failures

Failed pairwise commands or incomplete numerical exports get one additional
attempt by default: **two attempts in total**. `--ldsc_retries 2` allows three;
`0` disables retries. Exhausted execution failures stop analysis before final
collection. Extraction, munging and chi-square workers also have a second
attempt, independently of this flag.

Numerically invalid estimates follow `--result_failure_action`:

- `error` (default): save diagnostics and stop without a usable final CSV.
- `report`: save diagnostics; the command can exit successfully even when
  `LDSC_Compilation_Status.csv` reports failed estimates and no final CSV exists.
- `drop_traits`: save removals and a complete retained subset. Requires all
  `ref=yes` and at least two retained traits to produce a usable final CSV.

Structural errors remain fatal. Trait removal also checks duplicate estimates
for conflicts before choosing a subset; GPCA checks duplicate consistency before
analysis. Read the compilation status, pair and trait reports even when the
command finishes. Trait removal changes the downstream PCA; use the retained
manifest with the retained results.

## Results and next step

### The main output: `ldsc_results.csv`

For either input route, start with **`<outdir>/ldsc_results.csv`**. Each row is
a directed trait comparison; `p1 == p2` identifies a self-pair. Managed all-reference
runs retain both orientations of each between-trait pair. Rows follow manifest
trait order. Four traits normally give 16 rows; GPCA checks the 10 unique
unordered/self combinations and checks duplicate orientations for consistency.

These are the columns from current managed native exports. Optional h2 scales
can be empty, and their position can vary with the scales exported by the
batches; read columns by name. Older files may lack the three native fields.

| Column | Meaning |
| --- | --- |
| `p1` | First/reference trait name. |
| `p2` | Second/target trait name. |
| `rg` | Original Python genetic correlation, normalized with the two pair-specific heritabilities. Default correlation used by Python GPCA. |
| `se` | Standard error of original `rg`. |
| `z` | Test statistic for original `rg`. |
| `p` | P value for original `rg`. |
| `h2_int` | Univariate LDSC intercept for `p2` in this fit. Its self-pair value supplies the CTI diagonal. |
| `h2_int_se` | SE of `h2_int`. |
| `gcov_int` | Bivariate LDSC intercept. Between-trait values supply the CTI off-diagonal entries. This is not genetic covariance. |
| `gcov_int_se` | SE of `gcov_int`. |
| `h2_obs` | Unconverted heritability for `p2` in this fit, when available. Use its self-pair row for a trait-wide estimate. |
| `h2_obs_se` | SE of `h2_obs`. |
| `gcov_native_obs` | Native fitted genetic covariance, before liability conversion. |
| `h2_p1_pair_obs` | First trait's native unconverted heritability in this pair fit. |
| `h2_p2_pair_obs` | Second trait's native unconverted heritability in this pair fit. |
| `h2_liab` | Liability-scale heritability for `p2` when that conversion is requested and available. Empty for the quantitative example below. |
| `h2_liab_se` | SE of `h2_liab`. |
| `h2_scale` | Collection label: `NEF_unconverted` without population prevalence; `liability` when it is supplied. Check N provenance as well as this label. |
| `gcov_pair` | Covariance reconstructed from original `rg` and pair heritabilities on their reported scales. It can have different units from `gcov_native_obs` after liability conversion. |
| `rg_trait_wide` | Genetic correlation using native covariance and the two traits' self-pair heritabilities; self-correlations are set to 1. |
| `normalization_status` | `calculated`, or the reason derived values could not be calculated. Original estimates are retained. |

“Unconverted” means before liability conversion. With effective sample size
NEF, the name `obs` does not establish the usual population observed-scale
interpretation. Read the manifest and N provenance together.

For native exports, the difference between the two correlations is:

```text
rg(A,B) = gcov_native_obs(A,B) / sqrt(h2_A_pair_obs × h2_B_pair_obs)
rg_trait_wide(A,B) = gcov_native_obs(A,B) / sqrt(h2_A_self_obs × h2_B_self_obs)
```

The self estimates come from A–A and B–B. Collection adds the derived values
in memory; it does not rerun LDSC or change original `rg,se,z,p`. In particular,
`se,z,p` are **not new uncertainty estimates for `rg_trait_wide`**. When needed
inputs or scales are invalid, derived values stay empty and the status explains
why. See [normalization details and older-file behavior](REFERENCE.md#optional-trait-wide-correlation-normalization).

### A real-data result row

This is an actual between-trait row from `T3_C0_C1_C2`, shown vertically so all
columns are readable. Numerical strings below are copied from the compiled CSV.
It comes from saved regression fits on real munged data, recompiled with the
current collector on 2026-10-04; it is not a fresh VCF run.

| Column | Value |
| --- | --- |
| `p1` | `AGR2_O95994_OID20896_v1_Neurology` |
| `p2` | `AGRP_O00253_OID20658_v1_Inflammation` |
| `rg` | `0.21759828728629294` |
| `se` | `0.1161828427773906` |
| `z` | `1.872895189035932` |
| `p` | `0.061082859669125543` |
| `h2_int` | `0.99908304377925494` |
| `h2_int_se` | `0.0085624830776470633` |
| `gcov_int` | `0.11372629519394008` |
| `gcov_int_se` | `0.00518269692027908` |
| `h2_obs` | `0.15421852646762016` |
| `h2_obs_se` | `0.024544434823656189` |
| `gcov_native_obs` | `0.021553153240404802` |
| `h2_p1_pair_obs` | `0.063617139695829716` |
| `h2_p2_pair_obs` | `0.15421852646762016` |
| `h2_liab` | Empty |
| `h2_liab_se` | Empty |
| `h2_scale` | `NEF_unconverted` |
| `gcov_pair` | `0.021553153240404802` |
| `rg_trait_wide` | `0.21786011704514013` |
| `normalization_status` | `calculated` |

Here original `rg` and trait-wide `rg` differ slightly. Native and reconstructed
covariance agree on this unconverted scale; neither is the intercept `gcov_int`.
The [real-data verification record](PYTHON_LDSC_OUTPUTS.md#real-data-verification)
includes both checked datasets, row counts and the limits of that validation.

### Check results before GPCA

Read `LDSC_Compilation_Status.csv`, then `LDSC_Pair_Status.csv` and
`LDSC_Trait_Status.csv`. The diagnostic CSV keeps all collected estimates before
trait selection, including failed estimates; it is written on successful runs
too. Its existence alone does not establish that results are ready for GPCA.
Self-pairs can legitimately have `se=0,z=Inf,p=0` together when `rg` is within
the allowed tolerance of 1; this narrow exception is audited. Other required
SEs still need to be finite and positive.

Pass the final CSV and your selected manifest to [`ldsc-gpca gpca`](PYTHON_GPCA.md).
If traits were removed during collection, use `LDSC_Retained_Traits.csv`.
PCA/GWAMA also needs the separate per-SNP GWAMA inputs described in that guide.

Python GPCA defaults to original `rg` (`--rg_normalization pair`). Select
`--rg_normalization trait_wide` **on `gpca` or the Python pipeline** to use the
additional column. `--rg_normalization` is not an `ldsc` option. Merely writing
the column does not select it for PCA, and P values do not filter or weight PCA.

## All options

Required options have no default. “Off” means omit the flag; include the flag
alone to enable it. These options apply to `ldsc-gpca ldsc --ldsc_backend python`;
omitting the selector also chooses Python. Preparation,
PCA and GWAMA export settings belong to their respective commands.

### Files and input modes

| Option | Default | Meaning and input-mode applicability |
| --- | --- | --- |
| `--ldsc_backend` | `python` | `python` or `genomicsem`; selects matching munging/regression and backend-specific options. This guide covers Python. |
| `--input` | Required | Manifest CSV: VCF route needs `vcf_files`; reuse does not. |
| `--outdir` | Required | Results directory for either route. Use a new directory for a separate analysis. |
| `--ld_ref` | Required | Chromosome LD-score/M-file reference directory for either route. |
| `--ld_weights` | Use `--ld_ref` | Separate regression-weight LD-score directory, either route. |
| `--hm3` | Unset | SNP/A1/A2 allele reference; required for VCF munging, unused with `--ldsc_only`. |
| `--ldsc_only` | Off | Read existing munged inputs; skip extraction and munging. |
| `--munged_dir` | `<outdir>/ldsc_input` | VCF route: destination for new munged files. Reuse: source of existing files. |

### Extraction and munging filters

These filters act on VCF inputs. With `--ldsc_only`, they only describe/check
recorded provenance and do not filter the files again.

| Option | Default | Meaning |
| --- | --- | --- |
| `--info_min` | `0.7` | Inclusive minimum FORMAT/SI; accepted range `[0,1]`. |
| `--maf_min` | `0.01` | Keep FORMAT/AF between this value and 1 minus it, inclusive. Threshold must be strictly between 0 and 0.5. |
| `--munge_maf_min` | `0.005` | Python LDSC munging keeps MAF strictly greater than this value; threshold range `[0,0.5)`. |
| `--max_af_difference` | `0.2` | Maximum absolute INFO/AF − INFO/EUR difference; range `[0,1]`. Missing either value excludes the record. |
| `--exclude_mhc` | Off | Remove the inclusive interval defined below. |
| `--mhc_chr` | `6` | Chromosome for MHC exclusion. |
| `--mhc_start` | `25000000` | Inclusive start; positive integer. Check genome build. |
| `--mhc_end` | `35000000` | Inclusive end; at least the start position. |
| `--remove_palindrome` | Off | Remove A/T and C/G records in the extraction AF interval below. Python LDSC munging later removes all palindromic SNPs. |
| `--paliandromaf_lower` | `0.45` | Lower AF bound; use this exact flag spelling. Bounds must satisfy 0 ≤ lower ≤ upper ≤ 1. |
| `--paliandromaf_upper` | `0.55` | Upper AF bound; use this exact flag spelling. |

### Settings shared by both routes

| Option | Default | Meaning |
| --- | --- | --- |
| `--chisq_max` | Disabled | Positive integer or `auto`; filter separate munged copies before regression. |
| `--result_failure_action` | `error` | `error`, `report` or `drop_traits`; see [failure handling](#failed-estimates-and-execution-failures). |
| `--restart` | Off | Reuse verified completed LDSC batches in the same output directory. Does not skip extraction/munging by itself. |
| `--ldsc_retries` | `1` | Additional attempts per failed LDSC command/incomplete export; integer ≥ 0. Default is two total attempts. |
| `--n_cores` | `5` | Concurrent workers; positive integer. Applies to the stages this input route runs. |
| `--conda_executable` | `CONDA_EXE`, otherwise `conda` | Executable used to launch the child LDSC environment. |
| `--ldsc_env` | Configured prefix, otherwise `ldsc-cbiit` | Explicit child environment name; mutually exclusive with `--ldsc_env_prefix`. Normally omit both after installation. |
| `--ldsc_env_prefix` | Saved `LDSC_GPCA_LDSC_PREFIX`, if set | Explicit child environment directory overrides the saved prefix. Without a prefix, use the named environment. |
| `--help` | Off | Show usage, defaults and choices; exit without analysis. |

### VCF extraction executable

| Option | Default | Meaning |
| --- | --- | --- |
| `--bcftools` | `bcftools` on PATH | Executable name/path used for VCF extraction; unnecessary with `--ldsc_only`. |
