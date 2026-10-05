# Run the complete LDSC → GPCA → GWAMA pipeline

[README](../README.md) · [Input formats](INPUTS.md) · [Full reference](REFERENCE.md)

## In this guide

- [What pipeline does](#what-pipeline-does)
- [VCF inputs and references](#vcf-inputs-and-references)
- [Python LDSC from VCFs](#python-ldsc-from-vcfs)
- [GenomicSEM LDSC from VCFs](#genomicsem-ldsc-from-vcfs)
- [Migrate existing pipeline commands](#migrate-existing-pipeline-commands)
- [Run LDSC and PCA without GWAMA](#run-ldsc-and-pca-without-gwama)
- [Change filters, PCA and export settings](#change-filters-pca-and-export-settings)
- [What runs and where results go](#what-runs-and-where-results-go)
- [Failures and continuing work](#failures-and-continuing-work)
- [All pipeline options](#all-pipeline-options)

## What pipeline does

`ldsc-gpca pipeline` accepts **GWAS summary-statistics VCFs as its only study-data
inputs**. It prepares GWAMA SNP tables and backend-specific LDSC inputs, runs
LDSC and PCA, then runs PC1 GWAMA and exports the results. Select
`--ldsc_backend python` (default) or `--ldsc_backend genomicsem`.

The backend selects both munging and regression. `--munge_backend` belongs to
standalone `prepare`. Pipeline prepares the GWAMA tables automatically; no
`--gpca_input_folder` is accepted. Existing raw/munged files and completed LDSC
results belong to the [standalone commands](#migrate-existing-pipeline-commands).

<a id="choose-the-two-inputs-for-your-run"></a>

## VCF inputs and references

Use `/data/traits_vcf.csv`, with at least two traits in the required order:

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

Empty prevalence cells denote quantitative traits. Each VCF must contain exactly
one GWAS sample. These are summary-statistics VCFs, not individual-level genotype
VCFs. Relative paths resolve beside the manifest. Populated `sumstats_file` or
`munged_file` columns are rejected, even when VCF paths are also present.

For a full run, each VCF must satisfy **both** its LDSC contract and the GWAMA
preparation contract. GWAMA needs FORMAT `AF,ES,SE,LP,NEF`.
Python LDSC additionally needs FORMAT `SI,EZ` and INFO `AF,EUR`; population-
prevalence traits also need FORMAT `NC,NCO`. GenomicSEM LDSC needs FORMAT `SI`
and valid variant IDs; binary traits require an appropriate explicit manifest
`N` and both prevalences. See [all VCF fields and sample-size rules](INPUTS.md#vcf-fields).
The LDSC and GWAMA stages retain their own filtering and N conventions.

Supply `--hm3` with `SNP,A1,A2` allele-reference headers and `--ld_ref` with
chromosome LD scores/M files. `--ld_weights` optionally selects a separate weights
directory. These references remain required even though the study input is VCF-only.
For Python LDSC, pipeline supplies `ref=yes` if absent; every supplied value must
already be `yes` to obtain complete pair/self coverage.

**Choose an output directory that does not exist.** Pipeline has no overwrite
or whole-run resume mode. Rscript, bcftools and the selected LDSC runtime must
be available. The bundled GWAMA output lacks INFO: set the shell variable
`GWAMA_INFO` to a justified constant for final export. It is export metadata,
not measured imputation quality or a filtering threshold.

<a id="python-from-vcfs"></a>

## Python LDSC from VCFs

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
  --input /data/traits_vcf.csv \
  --outdir /results/four_traits_python \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

The original Python extraction/munging filters remain in force. GWAMA preparation
runs independently from the same VCFs. LDSC filters do not filter the GWAMA tables.

<a id="genomicsem-with-explicit-vcf-input"></a>
<a id="genomicsem-ldsc-with-explicit-vcf-input"></a>

## GenomicSEM LDSC from VCFs

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_vcf.csv \
  --outdir /results/four_traits_genomicsem \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

GenomicSEM VCF conversion is automatic: omit `--vcf_input`. It retains FORMAT/SI
as INFO in the raw tables for GenomicSEM munging. Quantitative LDSC N uses NEF
unless manifest N overrides it; binary traits require an appropriate explicit
N and both prevalences. GWAMA always keeps its original NEF. GenomicSEM's existing
two-pass LDSC and retained-trait selection are unchanged.

<a id="python-from-munged-files"></a>
<a id="python-ldsc-from-munged-files"></a>
<a id="genomicsem-from-raw-tables"></a>
<a id="genomicsem-ldsc-from-raw-tables"></a>
<a id="genomicsem-from-munged-files"></a>
<a id="genomicsem-ldsc-from-munged-files"></a>
<a id="genomicsem-from-quantitative-vcfs"></a>
<a id="genomicsem-ldsc-from-quantitative-vcfs"></a>

## Migrate existing pipeline commands

Version 0.8.0 introduces the VCF-only contract as an intentional CLI change. Pipeline rejects
`--ldsc_only`, `--munged_dir`, `--munged_input`, `--gpca_input_folder`,
`--ldsc_results` and `--vcf_input`; it does not silently ignore them.

| Your existing input | Use instead |
| --- | --- |
| Python LDSC munged files | [Standalone Python LDSC](PYTHON_LDSC.md#start-from-munged-files), then GPCA |
| GenomicSEM raw or compatible munged files | [Standalone GenomicSEM LDSC](GENOMICSEM_LDSC.md), then GPCA |
| Completed Python LDSC estimates | [Python LDSC results to GPCA/GWAMA](PYTHON_GPCA.md) |
| Completed GenomicSEM LDSC RData | [GenomicSEM LDSC results to GPCA/GWAMA](GENOMICSEM_GPCA.md) |
| Prepared GWAMA tables | Pass `--gpca_input_folder` to the appropriate standalone GPCA command |

For an existing VCF-based GenomicSEM pipeline command with `--vcf_input`, simply
remove that switch. The same INFO-preserving adapter now runs automatically.
The older implicit quantitative-VCF fallback (N=NEF, no INFO) is removed from
pipeline. Runs that used it can retain different SNPs because INFO filtering
now applies; binary inputs also require explicit N. No estimator source changed.

## Run LDSC and PCA without GWAMA

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
  --input /data/traits_vcf.csv \
  --outdir /results/python_ldsc_pca_only \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --chisq_max auto \
  --rg_normalization trait_wide \
  --validate_only
```

This still extracts VCFs, munges, runs LDSC and calculates PCA. It skips GWAMA,
export and GWAMA preparation unless `--write_munge_inputs` is requested.
GenomicSEM likewise converts VCFs inside its LDSC stage. For PCA from saved
LDSC results, use a standalone GPCA command.

## Change filters, PCA and export settings

This full Python example shows optional settings together:

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
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

For GenomicSEM LDSC, use `--info_filter` and `--maf_filter` for raw munging. Omit
`--chisq_max` to use GenomicSEM automatic filtering, or supply a positive number;
GenomicSEM does not accept the string `auto`. GenomicSEM GPCA uses its own `S_Stand` or
`S`, so Python's `--rg_normalization` and `--heritability_scale` do not apply.

`--pca_matrix covariance` is available for either backend, with different
matrix sources. Read [Python PCA choices](PYTHON_GPCA.md#choose-the-correlation-and-pca-method)
or [GenomicSEM PCA choices](GENOMICSEM_GPCA.md#choose-the-pca-method) first.

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
    munge_inputs/                  # only with --write_munge_inputs
  ldsc/
    ldsc_results.csv               # Python backend
    genomicsem_LDSC.RData           # GenomicSEM backend instead
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
retries default to one additional attempt. GenomicSEM munging and GWAMA workers
also get a second attempt, but numerical failures follow their stage's QC
policy. GenomicSEM LDSC initially drops invalid-h2 traits by default; Python result
collection and downstream GPCA default to stopping on invalid estimates.
Inspect all removal reports because changing traits changes PC1.

## All pipeline options

Use `ldsc-gpca pipeline --ldsc_backend python --help` for Python defaults or `ldsc-gpca pipeline --ldsc_backend genomicsem --help` for GenomicSEM defaults. Required options have no default. “Off” flags are omitted. Backend-specific flags must match the selected backend.

Preparation settings apply to VCFs. Every pipeline run extracts and munges for LDSC; full runs also prepare GWAMA tables. Export settings do not change the analysis and are skipped by `--validate_only`.

<a id="files-and-input-modes"></a>

### VCF manifest and references

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--input` | Both | **Required**. Manifest CSV. Requires traitname,vcf_files,sample_prevalence,population_prevalence; preserve trait order. |
| `--outdir` | Both | **Required**. Pipeline root directory; must not already exist. No whole-run resume. |
| `--ldsc_backend` | Both | **`python`**. `python` or `genomicsem`; selects LDSC and the matching GPCA reader. Backend-specific flags differ. |
| `--ld_ref` | Both | **Required**. Directory of chromosome LD scores and M reference files. |
| `--ld_weights` | Both | **Use `--ld_ref`**. Optional separate directory of regression-weight LD scores. |
| `--hm3` | Both | **Required**. Whitespace SNP/A1/A2 allele reference. Pipeline derives a SNP-only list only for optional legacy raw-table export. |

### VCF preparation and file layout

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--splitby_chr` | Both | **`split`**. `split`: all chromosomes 1–22 per trait. `nosplit`: one autosomal table per trait. Preparation and GWAMA must use matching layouts. |
| `--gpca_id_source` | Both | **`chr_pos_ref_alt`**. `chr_pos_ref_alt` or `vcf_id` for prepared GWAMA SNPIDs. Choose compatible IDs across traits. |
| `--p_min` | Both | **`1e-300`**. VCF preparation P floor, also used by automatic GenomicSEM VCF conversion. Does not change ES/SE-derived GWAMA Z; it can change extreme Z reconstructed by GenomicSEM munging. Finite and strictly between 0 and 1. |
| `--write_munge_inputs` | Both | **Off**. Also export legacy HapMap-selected raw tables with GWAMA preparation. Does not select the LDSC input route; LDSC still processes the VCFs. |
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
| `--remove_palindrome` | Python | **Off**. Python extraction: remove A/T and C/G variants within the specified AF interval. Python munging later removes all palindromic SNPs regardless. |
| `--paliandromaf_lower` | Python | **`0.45`**. Lower AF bound for extraction palindrome removal; use this exact option spelling. Bounds must satisfy 0 ≤ lower ≤ upper ≤ 1. |
| `--paliandromaf_upper` | Python | **`0.55`**. Upper AF bound for extraction palindrome removal; use this exact option spelling. |
| `--info_filter` | GenomicSEM | **`0.9`**. GenomicSEM raw munging: inclusive INFO threshold from 0 to 1, when a recognized INFO column exists. |
| `--maf_filter` | GenomicSEM | **`0.01`**. GenomicSEM raw munging: inclusive MAF threshold from 0 to 0.5, when recognized frequency data exist. |
| `--chisq_max` | Both | **Python: disabled; GenomicSEM: automatic**. Python accepts a positive integer or auto; GenomicSEM accepts a positive number and uses automatic filtering when omitted. Both act on LDSC, with different implementations; see recipes above. |
| `--chromosomes` | GenomicSEM | **`22`**. GenomicSEM LDSC: use reference chromosomes 1 through this integer (1–22). Does not change GWAMA split-file requirements. |
| `--n_blocks` | GenomicSEM | **`200`**. GenomicSEM LDSC: requested jackknife blocks, integer ≥2. Pinned GenomicSEM overrides the count for more than 18 traits; inspect its log. |

### PCA method and quality checks

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--rg_normalization` | Python | **`pair`**. Python GPCA: `pair` uses original rg; `trait_wide` requires valid precomputed `rg_trait_wide` and `normalization_status`. CTI and original columns stay unchanged. |
| `--pca_matrix` | Both | **`correlation`**. `correlation` or `covariance`. Python derives covariance from selected rg and self h2; GenomicSEM uses S. CTI is unchanged. |
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
| `--invalid_h2_action` | GenomicSEM | **`drop`**. GenomicSEM first-pass h2 QC: `drop` audits/removes nonfinite or nonpositive h2 traits; `error` stops. At least two must remain. |
| `--ldsc_retries` | Python | **`1`**. Additional attempts per failed Python LDSC command/incomplete export; integer ≥0. Default means two total attempts. Numerical failures follow QC policies. |

### Execution

| Option | Backend | Default and effect |
| --- | --- | --- |
| `--validate_only` | Both | **Off**. Run LDSC and then QC/PCA; skip GWAMA/export. VCF extraction/munging still runs. GWAMA preparation is skipped unless --write_munge_inputs is requested. Does not inspect GWAMA tables. |
| `--source_path` | Both | **Bundled modified GWAMA**. Optional custom R source defining the expected modified GWAMA function and output interface. Unused in validation-only mode. |
| `--n_cores` | Both | **`1` (GenomicSEM), `5` (Python)**. Positive integer forwarded to LDSC and GPCA/export. GenomicSEM LDSC regression remains sequential; --prepare_workers is separate. |
| `--bcftools` | Both | **`bcftools` on PATH**. Executable name/path; required for VCF extraction and preparation. |
| `--conda_executable` | Python | **`CONDA_EXE`, otherwise `conda`**. Executable name/path used to launch the child Python LDSC environment. |
| `--ldsc_env` | Python | **Configured prefix, otherwise `ldsc-cbiit`**. Explicit child environment name, mutually exclusive with `--ldsc_env_prefix`. The supplied installer configures a prefix, so normally omit both. |
| `--ldsc_env_prefix` | Python | **`LDSC_GPCA_LDSC_PREFIX` when set**. Explicit child environment directory overrides the saved prefix; mutually exclusive with `--ldsc_env`. Without a prefix, use the named environment. |
| `--rscript` | GenomicSEM | **`Rscript`**. GenomicSEM LDSC R executable/path. In pipeline this affects LDSC only; GPCA still needs Rscript on PATH. |
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
