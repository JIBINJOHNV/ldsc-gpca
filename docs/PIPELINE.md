# Run the complete pipeline

[README](../README.md) · [Input formats](INPUTS.md) · [Full reference](REFERENCE.md)

## In this guide

- [What pipeline does](#what-pipeline-does)
- [Choose the two inputs for your run](#choose-the-two-inputs-for-your-run)
- [Python from VCFs](#python-from-vcfs)
- [Python from munged files](#python-from-munged-files)
- [Native GenomicSEM from raw tables](#native-genomicsem-from-raw-tables)
- [Native GenomicSEM from munged files](#native-genomicsem-from-munged-files)
- [Native GenomicSEM with explicit VCF input](#native-genomicsem-with-explicit-vcf-input)
- [Native GenomicSEM from quantitative VCFs](#native-genomicsem-from-quantitative-vcfs)
- [Run LDSC and PCA without GWAMA](#run-ldsc-and-pca-without-gwama)
- [Change filters, PCA and export settings](#change-filters-pca-and-export-settings)
- [What runs and where results go](#what-runs-and-where-results-go)
- [Failures and continuing work](#failures-and-continuing-work)
- [All pipeline options](#all-pipeline-options)

## What pipeline does

`ldsc-gpca pipeline` coordinates preparation, LDSC, PCA, PC1 GWAMA and final
export. Python LDSC is the default; `--ldsc_backend genomicsem` selects native
GenomicSEM LDSC and its matching GPCA reader.

Use pipeline when you still need to **run LDSC**. If you already have completed
LDSC results, start with [Python GPCA](PYTHON_GPCA.md) or
[native GPCA](GENOMICSEM_GPCA.md) instead. Pipeline does not accept completed
results as a shortcut around its LDSC stage.

## Choose the two inputs for your run

First choose what LDSC will read:

- **Python:** supported VCFs, or existing Python munged files selected with
  both `--ldsc_only` and `--munged_dir`.
- **GenomicSEM:** raw tables listed in `sumstats_file`, existing munged files
  selected with `--munged_dir` or `--munged_input`, or VCFs selected with
  `--vcf_input`. The older quantitative-trait VCF fallback is also retained.

Then choose what GWAMA will read:

- Pass existing per-trait nine-column tables through `--gpca_input_folder`; or
- include suitable `vcf_files` in the manifest so pipeline prepares the tables.

Raw/munged LDSC tables do not supply the full GWAMA input. If you want to stop
after LDSC and PCA, use `--validate_only`; no GWAMA source is needed unless VCF
preparation is also needed to create native LDSC inputs or optional raw tables.

For every route, use at least two traits and retain
`traitname,population_prevalence,sample_prevalence` headers. Paths inside the
manifest resolve beside it. The Python pipeline supplies `ref=yes` when absent
and requires all supplied `ref` values to be `yes`. The same manifest row order
is carried through any audited trait removal.

**Choose an output directory that does not exist.** Pipeline has no overwrite
or whole-run resume mode. Use `Rscript` on PATH for GPCA and the configured
backend runtime for LDSC. See [reference requirements](INPUTS.md#reference-files).

## Python from VCFs

Use `/data/traits_vcf.csv`:

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

These are quantitative examples. Each VCF must meet **both** the Python LDSC
extraction and GWAMA preparation [field requirements](INPUTS.md#vcf-fields).

```bash
ldsc-gpca pipeline \
  --input /data/traits_vcf.csv \
  --outdir /results/four_traits_python \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Set `GWAMA_INFO` to your justified export constant first; the bundled GWAMA
function has no INFO output. This is summary metadata, not a filter or quality
estimate. Read [export metadata](INPUTS.md#export-metadata) before choosing it.
The Python pipeline's preparation and LDSC extraction use separate VCF paths
through the code and separate filtering rules. LDSC filters do not also filter
the GWAMA tables.

## Python from munged files

Use the [Python reuse manifest](INPUTS.md#reusing-munged-files), preserve its
required provenance sidecars, and supply the GWAMA tables separately:

```bash
ldsc-gpca pipeline \
  --input /data/traits_python_munged.csv \
  --outdir /results/python_reused_pipeline \
  --ld_ref /references/eur_w_ld_chr \
  --ldsc_only \
  --munged_dir /data/python_munged \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Both reuse flags are required in pipeline. Extraction and munging are skipped,
so their filters are not reapplied. `--chisq_max`, if supplied, still filters
separate LDSC input copies. Total-N traits need provenance sidecars; legacy NEF
inputs can produce a warning. See [Python reuse details](PYTHON_LDSC.md#reuse-munged-files-or-restart).

If you have VCFs instead of prepared GWAMA tables, add `vcf_files` to this
manifest and omit `--gpca_input_folder`. This prepares GWAMA inputs while still
reusing the selected munged files for LDSC.

## Native GenomicSEM from raw tables

Use `/data/traits_native_raw.csv` with
`traitname,sumstats_file,population_prevalence,sample_prevalence`; see the
[four-trait example](INPUTS.md#raw-table-manifest-for-genomicsem).
Provide prepared GWAMA tables separately:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_native_raw.csv \
  --outdir /results/native_raw_pipeline \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Native munging reads the raw tables, including suitable N/effect/P columns and
any recognized INFO/frequency fields. A positive manifest N replaces that
trait's raw-file N. Binary traits require both prevalence values; this wrapper
does not infer one. See [native input and sample-size rules](GENOMICSEM_LDSC.md).

To prepare GWAMA tables from VCFs in this route, add a `vcf_files` column alongside
`sumstats_file` and omit `--gpca_input_folder`. LDSC still uses the supplied raw
tables; the VCFs supply GWAMA statistics.

## Native GenomicSEM from munged files

For a single directory, use `/data/traits_native_munged.csv` with trait names
and the two prevalence columns:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_native_munged.csv \
  --outdir /results/native_munged_pipeline \
  --ld_ref /references/eur_w_ld_chr \
  --munged_dir /data/native_munged \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Alternatively, supply `munged_file` paths in the manifest and use
`--munged_input` in place of `--munged_dir /data/native_munged`. For example:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_native_munged_paths.csv \
  --outdir /results/native_paths_pipeline \
  --ld_ref /references/eur_w_ld_chr \
  --munged_input \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

See the [explicit-path manifest](GENOMICSEM_LDSC.md#reuse-munged-files).
Do not combine native reuse flags or add Python's `--ldsc_only`. Native reuse
skips munging; it uses the N values already in the files and does not apply raw
INFO/MAF filters again. To prepare GWAMA inputs instead, include VCF paths and
omit `--gpca_input_folder` as described above.

## Native GenomicSEM with explicit VCF input

Select `--vcf_input` to convert VCFs directly in the native LDSC stage, retaining
FORMAT/SI as INFO for native filtering:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --vcf_input \
  --input /data/traits_vcf.csv \
  --outdir /results/native_vcf_with_info \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

The native LDSC VCFs need `AF,ES,SE,LP,SI` and NEF unless a manifest N override
is supplied. Binary VCFs require an explicit appropriate N and both prevalences;
see the [native VCF contract](GENOMICSEM_LDSC.md#start-from-gwas-vcfs).
For GWAMA, supply prepared tables through `--gpca_input_folder`, or let pipeline
prepare them from the VCFs. That separate preparation still requires NEF and
uses it for GWAMA N, independently of a native LDSC N override.

With `--validate_only`, this route converts the VCFs for LDSC and runs LDSC/PCA
without preparing GWAMA inputs. The conversion writes its raw inputs and
`GenomicSEM_VCF_*` audits under `ldsc/`; native munging and final RData follow
the existing workflow. `--p_min` applies to VCF conversion in LDSC and, when
run, to the separate GWAMA preparation. Native extraction uses `--n_cores`;
GWAMA preparation uses `--prepare_workers`.

## Native GenomicSEM from quantitative VCFs

Without `--vcf_input`, if no raw `sumstats_file` paths and no munged reuse mode are supplied, the native
pipeline can turn quantitative-trait VCFs into raw tables for native munging:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_vcf.csv \
  --outdir /results/native_vcf_pipeline \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

VCFs need the **preparation** fields `AF,ES,SE,LP,NEF`. Pipeline automatically
requests the additional raw tables, then munges them with GenomicSEM. These
raw tables use N=NEF and contain **no INFO**, so `--info_filter` cannot provide
INFO-based filtering for this route. Existing GWAMA tables can be supplied,
but the VCFs are still needed to create the LDSC raw inputs.

This older fallback refuses binary traits because the appropriate native N
convention cannot be inferred. Use the explicit `--vcf_input` route with an
appropriate manifest N and both prevalences, or supply raw/munged native inputs.
The fallback is a pipeline feature; standalone `genomicsem ldsc` requires
`--vcf_input` to read VCFs.

## Run LDSC and PCA without GWAMA

```bash
ldsc-gpca pipeline \
  --input /data/traits_vcf.csv \
  --outdir /results/python_ldsc_pca_only \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --chisq_max auto \
  --rg_normalization trait_wide \
  --validate_only
```

This **runs LDSC**, then checks its results and calculates PCA. It skips GWAMA,
SNP-result export and ordinary GWAMA-only preparation. It does not check GWAMA
SNP tables. The older native VCF fallback still needs preparation to create LDSC raw
tables; requesting `--write_munge_inputs` also requires preparation. For PCA
from already completed results, use a standalone GPCA command instead.

## Change filters, PCA and export settings

This full Python example shows optional settings together:

```bash
ldsc-gpca pipeline \
  --input /data/traits_vcf.csv \
  --outdir /results/python_optional_settings \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --info_min 0.9 \
  --maf_min 0.01 \
  --exclude_mhc \
  --chisq_max auto \
  --rg_normalization trait_wide \
  --pca_matrix correlation \
  --pc1_orientation tutorial \
  --splitby_chr nosplit \
  --prepare_workers 4 \
  --n_cores 8 \
  --dataset_id four_traits_custom \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}" \
  --gwama_output_n_eff "${GWAMA_N_EFF:?Set a justified summary N_eff}" \
  --gzip_level 1 \
  --archive_chromosomes
```

These choices are illustrative, not a recommendation to override every default.
Omit `--gwama_output_n_eff` to retain GWAMA's reported N_eff. The override affects
only the selected-column summary; it does not change calculation N, Z, P, BETA,
SE or the full combined output.

For native input, use `--info_filter` and `--maf_filter` for raw munging. Omit
`--chisq_max` to use native automatic filtering, or supply a positive number;
native does not accept the string `auto`. Native GPCA uses its own `S_Stand` or
`S`, so Python's `--rg_normalization` and `--heritability_scale` do not apply.

`--pca_matrix covariance` is available for either backend, with different
matrix sources. Read [Python PCA choices](PYTHON_GPCA.md#choose-the-correlation-and-pca-method)
or [native PCA choices](GENOMICSEM_GPCA.md#choose-the-pca-method) first.

## What runs and where results go

The pipeline validates inputs, creates stage manifests, runs preparation when
needed, runs LDSC, selects the retained traits, then invokes the matching GPCA
command. A full run continues through GWAMA and export. Stage commands,
timings and failures are recorded.

```text
/results/four_traits_python/
  Pipeline_Run_Status.json
  manifests/
  prepare/                         # when preparation is needed
    gpca_inputs/
    munge_inputs/                  # only when requested/needed
  ldsc/
    ldsc_results.csv               # Python backend
    genomicPCA_LDSC.RData           # native backend instead
  gpca/
    GenomicPCA_PC1_Weights_Used.csv
    <dataset_id>_GWAMA_combined_results.txt.gz
    harmonisation_input/<dataset_id>_GPCA_inputs.txt.gz
    <dataset_id>_postprocess.json
```

Only the selected backend's LDSC result exists. Validation-only runs do not
produce the GWAMA/export files. Pipeline's default `dataset_id` is its output
directory name. `manifests/gpca_traits.csv` records the LDSC-retained trait list
passed to GPCA; inspect downstream QC for any further GPCA removals. See the
[full output catalog](REFERENCE.md#find-and-interpret-the-outputs).

## Failures and continuing work

Pipeline stops at an exhausted stage failure and records it in
`Pipeline_Run_Status.json`. Inspect the failed stage's log; do not treat the
presence of intermediate files as successful completion.

There is no `pipeline --restart`. Use individual commands to continue from
verified outputs: standalone [Python LDSC restart](PYTHON_LDSC.md#reuse-munged-files-or-restart)
can reuse checked batches; a completed LDSC result can be passed to its matching
GPCA command with a fresh GPCA output directory. Repeat your intended analysis
settings when continuing manually. Do not assume pipeline metadata supplies
them as defaults.

Execution retries and scientific QC are different. Python pairwise command
retries default to one additional attempt. Native munging and GWAMA workers
also get a second attempt, but numerical failures follow their stage's QC
policy. Native LDSC initially drops invalid-h2 traits by default; Python result
collection and downstream GPCA default to stopping on invalid estimates.
Inspect all removal reports because changing traits changes PC1.

## All pipeline options

Use `ldsc-gpca pipeline --help` for Python defaults or `ldsc-gpca pipeline --ldsc_backend genomicsem --help` for native defaults. Required options have no default. “Off” flags are omitted. Backend-specific flags must match the selected backend.

Preparation settings act only when that stage runs. Raw extraction/munging filters are not reapplied to existing munged files. Export settings do not change the analysis and are skipped by `--validate_only`.

### Files and input modes

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--input` | Both | **Required**. Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Both | **Required**. Pipeline root directory; must not already exist. No whole-run resume. |
| `--ldsc_backend` | Both | **`python`**. `python` or `genomicsem`; selects LDSC and the matching GPCA reader. Backend-specific flags differ. |
| `--ld_ref` | Both | **Required**. Directory of chromosome LD scores and M reference files. |
| `--ld_weights` | Both | **Use `--ld_ref`**. Optional separate directory of regression-weight LD scores. |
| `--hm3` | Both | **Unset**. Whitespace SNP/A1/A2 allele reference, required for munging or optional raw-table preparation. Pipeline derives the SNP-only preparation list. |
| `--ldsc_only` | Python | **Off**. Reuse Python munged files; skip extraction/munging. Their INFO/MAF and other extraction filters are not reapplied. |
| `--munged_dir` | Both | **Unset**. Python: requires --ldsc_only too. Native: selects directory reuse by itself, mutually exclusive with --munged_input. See backend file/provenance rules. |
| `--munged_input` | GenomicSEM | **Off**. Native reuse from manifest `munged_file` paths; mutually exclusive with `--munged_dir`. |
| `--vcf_input` | GenomicSEM | **Off**. Convert manifest `vcf_files` to raw tables with INFO during the native LDSC stage. Mutually exclusive with either munged mode; binary traits require an appropriate explicit manifest N. |
| `--gpca_input_folder` | Both | **Unset**. Directory of existing nine-column GWAMA tables. Otherwise prepare from manifest VCFs for a full run. Not needed for ordinary validation-only runs. |

### VCF preparation and file layout

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--splitby_chr` | Both | **`split`**. `split`: all chromosomes 1–22 per trait. `nosplit`: one autosomal table per trait. Preparation and GWAMA must use matching layouts. |
| `--gpca_id_source` | Both | **`chr_pos_ref_alt`**. `chr_pos_ref_alt` or `vcf_id` for prepared GWAMA SNPIDs. Choose compatible IDs across traits. |
| `--p_min` | Both | **`1e-300`**. VCF preparation P floor, also used by native `--vcf_input`. Does not change ES/SE-derived GWAMA Z; it can change extreme Z reconstructed by native munging. Finite and strictly between 0 and 1. |
| `--write_munge_inputs` | Both | **Off**. Also request HapMap-selected raw tables. Automatically enabled for the older native VCF fallback. Requires suitable VCFs and --hm3 even in validation-only mode. |
| `--munge_id_source` | Both | **`vcf_id`**. `vcf_id` or `chr_pos_ref_alt` for optional raw-table SNP IDs; must match HapMap IDs. |
| `--prepare_workers` | Both | **`4`**. Positive integer; workers for VCF preparation, separate from LDSC/GWAMA workers. |

### LDSC filtering and estimation

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--info_min` | Python | **`0.7`**. Python extraction: inclusive minimum FORMAT/SI, from 0 to 1. |
| `--maf_min` | Python | **`0.01`**. Python extraction: retain FORMAT/AF between MAF and 1−MAF, inclusive. Threshold strictly between 0 and 0.5. |
| `--munge_maf_min` | Python | **`0.005`**. Python munging: require MAF strictly greater than this value; threshold from 0 up to, but not including, 0.5. |
| `--max_af_difference` | Python | **`0.2`**. Python extraction: maximum absolute INFO/AF minus INFO/EUR difference, from 0 to 1. Missing fields are excluded and counted. |
| `--exclude_mhc` | Python | **Off**. Python extraction: remove the inclusive interval defined by the three MHC options below. |
| `--mhc_chr` | Python | **`6`**. Chromosome for MHC exclusion; relevant when `--exclude_mhc` is enabled. |
| `--mhc_start` | Python | **`25000000`**. Inclusive MHC start, positive integer; check your genome build. |
| `--mhc_end` | Python | **`35000000`**. Inclusive MHC end, at least the start position. |
| `--remove_palindrome` | Python | **Off**. Python extraction: remove A/T and C/G variants within the specified AF interval. Native Python munging later removes all palindromic SNPs regardless. |
| `--paliandromaf_lower` | Python | **`0.45`**. Lower AF bound for extraction palindrome removal; use this exact option spelling. Bounds must satisfy 0 ≤ lower ≤ upper ≤ 1. |
| `--paliandromaf_upper` | Python | **`0.55`**. Upper AF bound for extraction palindrome removal; use this exact option spelling. |
| `--info_filter` | GenomicSEM | **`0.9`**. Native raw munging: inclusive INFO threshold from 0 to 1, when a recognized INFO column exists. Inactive for munged reuse. |
| `--maf_filter` | GenomicSEM | **`0.01`**. Native raw munging: inclusive MAF threshold from 0 to 0.5, when recognized frequency data exist. Inactive for munged reuse. |
| `--chisq_max` | Both | **Python: disabled; native: automatic**. Python accepts a positive integer or auto; native accepts a positive number and uses automatic filtering when omitted. Both act on LDSC, with different implementations; see recipes above. |
| `--chromosomes` | GenomicSEM | **`22`**. Native LDSC: use reference chromosomes 1 through this integer (1–22). Does not change GWAMA split-file requirements. |
| `--n_blocks` | GenomicSEM | **`200`**. Native LDSC: requested jackknife blocks, integer ≥2. Pinned GenomicSEM overrides the count for more than 18 traits; inspect its log. |

### PCA method and quality checks

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--rg_normalization` | Python | **`pair`**. Python GPCA: `pair` uses original rg; `trait_wide` requires valid precomputed `rg_trait_wide` and `normalization_status`. CTI and original columns stay unchanged. |
| `--pca_matrix` | Both | **`correlation`**. `correlation` or `covariance`. Python derives covariance from selected rg and self h2; native uses S. CTI is unchanged. |
| `--heritability_scale` | Python | **`auto`**. Python GPCA: `auto`, `observed`, `liability`, `mixed`. Selects existing h2 columns, never converts scales; covariance auto needs a common complete scale. |
| `--pc1_orientation` | Both | **`tutorial`**. `tutorial` flips all PC1 loadings if their median is negative; `as_computed` keeps the eigenvector sign. |
| `--allow_missing_traits` | Both | **Off**. GPCA: permit audited removal of manifest traits absent from LDSC results. Otherwise stop. |
| `--failed_ldsc_action` | Both | **`error`**. GPCA: `error` or `drop_traits` for invalid/missing estimates. Dropping changes the analysed trait set; no imputation. |
| `--h2_z_warn_threshold` | Both | **`2`**. GPCA: warn below this retained-trait h2/SE ratio. Nonnegative; 0 disables. Diagnostic only, never a removal rule. |
| `--rg_out_of_range_action` | Both | **`warn`**. GPCA: `warn` or `error` for finite off-diagonal rg outside [−1,1]. Never clamps values. |
| `--negative_eigen_action` | Both | **`warn`**. GPCA: `warn` or `error` for substantive negative PCA eigenvalues. No silent matrix repair. |
| `--matrix_eigen_tolerance` | Both | **`1e-8`**. GPCA: positive relative tolerance separating substantive negative eigenvalues from floating-point noise. |
| `--duplicate_tolerance` | Python | **`0.001`**. Python GPCA: positive absolute tolerance for duplicate rg, SE, P and intercept estimates. |
| `--duplicate_z_tolerance` | Python | **`0.01`**. Python GPCA: positive absolute tolerance for duplicate-orientation Z values. |
| `--self_rg_tolerance` | Python | **`0.01`**. Python GPCA: positive tolerance around self rg=1. SE=0 is accepted only with Z=+Inf and P=0; other checks still apply. |
| `--comparison_epsilon` | Python | **`1e-12`**. Python GPCA: nonnegative floating-point allowance added to boundary comparisons. |
| `--z_consistency_tolerance` | Python | **`0.01`**. Python GPCA: positive relative tolerance when comparing reported Z with rg/SE after duplicate collapse. |
| `--z_consistency_action` | Python | **`warn`**. Python GPCA: `warn` or `error` when reported Z is inconsistent with rg/SE. |

### Failures and reuse

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--result_failure_action` | Python | **`error`**. Python LDSC collection: `error`, `report`, `drop_traits`. All save diagnostics. Report may leave no final CSV; drop needs all ref=yes and at least two retained traits. Structural conflicts remain fatal. |
| `--invalid_h2_action` | GenomicSEM | **`drop`**. Native first-pass h2 QC: `drop` audits/removes nonfinite or nonpositive h2 traits; `error` stops. At least two must remain. |
| `--ldsc_retries` | Python | **`1`**. Additional attempts per failed Python LDSC command/incomplete export; integer ≥0. Default means two total attempts. Numerical failures follow QC policies. |

### Execution

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--validate_only` | Both | **Off**. Run LDSC and then QC/PCA; skip GWAMA/export. Preparation still runs if needed for native VCF LDSC inputs or --write_munge_inputs. Does not inspect GWAMA tables. |
| `--source_path` | Both | **Bundled modified GWAMA**. Optional custom R source defining the expected modified GWAMA function and output interface. Unused in validation-only mode. |
| `--n_cores` | Both | **`1` (native), `5` (Python)**. Positive integer forwarded to LDSC and GPCA/export. Native LDSC regression remains sequential; --prepare_workers is separate. |
| `--bcftools` | Both | **`bcftools` on PATH**. Executable name/path for VCF extraction or preparation; unnecessary when neither runs. |
| `--conda_executable` | Python | **`CONDA_EXE`, otherwise `conda`**. Executable name/path used to launch the child Python LDSC environment. |
| `--ldsc_env` | Python | **Configured prefix, otherwise `ldsc-cbiit`**. Explicit child environment name, mutually exclusive with `--ldsc_env_prefix`. The supplied installer configures a prefix, so normally omit both. |
| `--ldsc_env_prefix` | Python | **`LDSC_GPCA_LDSC_PREFIX` when set**. Explicit child environment directory overrides the saved prefix; mutually exclusive with `--ldsc_env`. Without a prefix, use the named environment. |
| `--rscript` | GenomicSEM | **`Rscript`**. Native LDSC R executable/path. In pipeline this affects LDSC only; GPCA still needs Rscript on PATH. |
| `--ldsc_chunk_size` | Python | **`250000`**. Python GPCA: positive number of LDSC rows read per streaming chunk; controls memory/read size, not statistics. |

### Final export

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--dataset_id` | Both | **Output directory name**. Prefix for final combined/summary files and export audit. Pipeline uses the pipeline outdir name. |
| `--gwama_output_info` | Both | **Preserve reported INFO**. Optional finite constant in [0,1] for the selected-column summary only. Required for bundled GWAMA because it has no INFO. Does not filter variants. |
| `--gwama_output_n_eff` | Both | **Preserve reported N_eff**. Optional finite positive constant for the selected-column summary only. Does not change calculation N, weights, Z, P, BETA or SE. |
| `--gzip_level` | Both | **`1`**. Integer 1–9 for final export only: 1 faster, 9 smaller. Uses pigz when installed; otherwise single-worker Python gzip. |
| `--archive_chromosomes` | Both | **Off**. After successful export, move current-run source result/log files into outdir/chromosome_wise. |

### Help

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--help` | Both | **Off**. Show command usage, defaults and choices; exit without analysis. |
| `--prepare_help` | Both | **Off**. Show shared preparation options and exit without analysis. |
| `--postprocess_help` | Both | **Off**. Show automatic export options and exit without analysis. |
