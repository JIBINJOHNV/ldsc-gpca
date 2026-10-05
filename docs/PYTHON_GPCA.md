# PCA and PC1 GWAMA from Python LDSC results

[README](../README.md) · [Inputs](INPUTS.md) · [Python LDSC](PYTHON_LDSC.md) · [GenomicSEM GPCA](GENOMICSEM_GPCA.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Files you need](#files-you-need)
- [Check results and calculate PCA](#check-results-and-calculate-pca)
- [What happens in order](#what-happens-in-order)
- [Run PC1 GWAMA](#run-pc1-gwama)
- [Choose the correlation and PCA method](#choose-the-correlation-and-pca-method)
- [Prepare GWAMA tables automatically](#prepare-gwama-tables-automatically)
- [QC choices and failure handling](#qc-choices-and-failure-handling)
- [Exports and what to inspect](#exports-and-what-to-inspect)
- [All options](#all-options)

## When to use this command

Use `ldsc-gpca gpca` after Python LDSC has finished. It reads completed pairwise
estimates, checks them and calculates PCA. A full run then combines per-SNP
statistics with the PC1 loadings using GWAMA and exports the results. It does
not rerun LDSC. Both GPCA routes require `Rscript` on PATH.

## Files you need

Always provide:

- `--input`: a CSV with unique `traitname` values in the desired order, such as
  the [four-trait selection manifest](INPUTS.md#selecting-traits-from-completed-ldsc-results).
- `--ldsc_results`: the completed Python LDSC table, not munged per-trait files.
- `--outdir`: a destination for QC, matrices, weights and any GWAMA results.
  Use a fresh directory for each analysis choice.

The LDSC table can be CSV, tab- or whitespace-delimited, optionally gzipped.
Required correlation/intercept columns are:

```text
p1 p2 rg se z p h2_int h2_int_se gcov_int gcov_int_se
```

It also needs the appropriate self-pair heritability estimate/SE columns:
`h2_obs,h2_obs_se` or `h2_liab,h2_liab_se`. For k selected traits there must be
`k × (k+1) / 2` unique self/pair combinations: ten for Trait_A through Trait_D.
A triangle or a consistent symmetric table is accepted. Conflicting duplicate
orientations are fatal. See the [full table contract](REFERENCE.md#python-ldsc-results-table).

For **GWAMA as well**, either supply the nine-column per-trait tables using
`--gpca_input_folder`, or add suitable `vcf_files` to the manifest so the command
can prepare them. The VCFs provide GWAMA statistics; completed LDSC results are
still required. Check [GWAMA columns, names and alleles](INPUTS.md#gwama-tables).

## Check results and calculate PCA

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_pca_check \
  --validate_only
```

This writes QC and PCA reports without preparing SNP tables, running GWAMA or
exporting SNP results. **It does not validate the per-SNP GWAMA files.** Review
warnings, retained traits and weights before running the full analysis.

## What happens in order

1. Read selected traits and LDSC rows; check pair coverage, numeric values,
   self-pairs and duplicate orientations.
2. Apply any explicitly requested trait-removal policies, preserving manifest
   order. Do not impute missing estimates.
3. Build the genetic matrix and intercept matrix. Correlation PCA uses diagonal
   1 and the selected pair correlations off diagonal. The intercept matrix
   (CTI) uses self `h2_int` and off-diagonal `gcov_int`.
4. Check symmetry, positive-definite CTI and matrix eigenvalues. Calculate PC1
   loadings using symmetric eigen decomposition and write audit files.
5. Unless `--validate_only`, run GWAMA with those loadings, CTI and the retained
   traits' SNP tables; then combine and export successful current-run results.

When automatic VCF preparation is requested, it runs before the R QC/trait
selection stage. A VCF preparation error can therefore stop a full command
even if that trait would later be dropped.

P, SE and Z are diagnostics, not PCA weights or significance filters. Substantive
negative genetic eigenvalues are reported; the matrix is not silently replaced with a positive-definite approximation. Nonpositive eigenvalues are guarded when forming loadings and
reported, and PC1 must have a finite positive eigenvalue.

## Run PC1 GWAMA

After checking the same trait set and analysis choices, run in a fresh directory:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/python_pc1 \
  --splitby_chr nosplit \
  --dataset_id four_traits_python \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

The prepared filenames must match `nosplit`. For files created in the default
split layout, use `--splitby_chr split` (or omit it) and provide all 22 chromosome
files per trait. The bundled GWAMA function lacks INFO, so set a justified
export constant before the full command. See [export metadata](INPUTS.md#export-metadata).

GWAMA uses the available traits at each SNP, with weights based on
`sqrt(N) × PC1 loading` and the LDSC intercept matrix accounting for correlated
errors. Keep input IDs, allele conventions and genome builds compatible. Export
N_eff overrides do not change this weighting or recalibrate BETA/SE.

## Choose the correlation and PCA method

Two independent settings control this choice:

| Setting | Default | Alternative and meaning |
| --- | --- | --- |
| `--rg_normalization` | `pair`: original Python rg | `trait_wide`: precomputed covariance normalized by each trait's self-pair observed h2. |
| `--pca_matrix` | `correlation`: unit diagonal | `covariance`: selected rg multiplied by the square roots of selected self heritabilities; self h2 on the diagonal. |

For a pair A–B, original Python rg uses heritabilities estimated in that pair's
regression. Trait-wide rg uses one observed self-pair h2 for A and one for B
across all their pairings. The managed LDSC stage writes `rg_trait_wide` and
`normalization_status` by default; GPCA uses them only when requested. All
selected rows must have a successfully calculated normalization status.

For example, use trait-wide correlations:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_trait_wide_check \
  --rg_normalization trait_wide \
  --pca_matrix correlation \
  --validate_only
```

Original rg/SE/Z/P stay unchanged. The original SE/Z/P are not uncertainty
estimates for the new rg. CTI also stays unchanged. This normalization does
not reproduce GenomicSEM's regression weighting and cannot ensure identical
GenomicSEM results. For older result files, see
[the annotation utility and covariance definitions](REFERENCE.md#optional-trait-wide-correlation-normalization).

To use covariance PCA with observed-scale self heritabilities:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_covariance_check \
  --rg_normalization pair \
  --pca_matrix covariance \
  --heritability_scale observed \
  --validate_only
```

This requires populated observed-scale h2 columns. `--heritability_scale`
selects existing columns; it does not convert scales. `auto` uses available
per-trait scales for correlation PCA, but requires one complete common scale
for covariance PCA. `mixed` explicitly permits trait-specific covariance scales
only when each trait has one available scale. Ambiguous choices stop the run.
Consider trait scales carefully when interpreting covariance PCA.

The covariance matrix is constructed from the selected rg and self h2. There
is no separate PCA mode that directly selects the `gcov_pair` column. With
trait-wide rg and observed self h2, this construction recovers its covariance
numerator, subject to numeric precision.

To run GWAMA with either alternative, repeat its choices in a full command,
adding the per-SNP input and export options shown above. QC-only settings are
not saved as automatic defaults for a later invocation.

## Prepare GWAMA tables automatically

Include `vcf_files` in the manifest and omit `--gpca_input_folder`:

```bash
ldsc-gpca gpca \
  --input /data/traits_vcf.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_pc1_auto_prepare \
  --splitby_chr nosplit \
  --gpca_id_source chr_pos_ref_alt \
  --p_min 1e-300 \
  --prepare_workers 4 \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Preparation needs FORMAT `AF,ES,SE,LP,NEF`; it does not apply LDSC INFO/MAF
filters. `--prepare_workers` controls this stage separately from `--n_cores`.
Do not request nondefault preparation options while supplying existing GWAMA
tables to this standalone command. See [preparation](PREPARE.md).

## QC choices and failure handling

Missing traits and invalid selected estimates stop by default.
`--allow_missing_traits` permits audited removal of absent traits;
`--failed_ldsc_action drop_traits` permits removal to obtain a valid complete
subset. Both change what is analysed; neither repairs conflicting duplicates.
Low h2/SE is a diagnostic warning at the default threshold 2, not automatic
trait exclusion.

A self-pair can have rg near 1, SE=0, Z=+Inf and P=0. This dedicated self-pair
case is accepted and retained in diagnostics; other required SEs and nonself
estimates must satisfy normal checks. See
[QC policies](REFERENCE.md#understand-qc-and-scale-choices).

The default PC1 sign convention flips the whole loading vector if its median
is negative. `--pc1_orientation as_computed` preserves R's eigenvector sign.
Signs are arbitrary; changing the whole-vector sign does not change association
strength or P values.

Failed chromosome GWAMA jobs retry once; failure after two attempts stops.
`--n_cores 0` selects workers automatically (Windows runs sequentially). Choose
fresh output directories for comparisons: completed export files are not
silently overwritten, and failed-run outputs should be inspected before reuse.

## Exports and what to inspect

For the full example, check `GWAMA_Run_Status.csv`, the matrices and
`GenomicPCA_PC1_Weights_Used.csv`, then:

- `four_traits_python_GWAMA_combined_results.txt.gz`: original combined output.
- `harmonisation_input/four_traits_python_GPCA_inputs.txt.gz`: selected-column
  result summary.
- `four_traits_python_postprocess.json`: export audit.

To override only the summary's N_eff, add
`--gwama_output_n_eff "${GWAMA_N_EFF:?Set a justified summary N_eff}"` to the
full command. Omit it to preserve reported values. `--gzip_level 1` is fastest;
`9` favours smaller final files. `--archive_chromosomes` moves current-run
source results and logs into `chromosome_wise` after export.

Python pairwise results do not contain GenomicSEM full sampling covariance matrices
V/V_Stand. They are sufficient for this PC1/GWAMA calculation, but cannot be
used to fabricate those matrices or run methods requiring them, such as paLDSC.
See the [output catalog and interpretation](REFERENCE.md#find-and-interpret-the-outputs).

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

Preparation settings apply only when VCF preparation runs. Export settings apply only after a full successful GWAMA run; both stages are skipped by `--validate_only`.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Output directory. Use a fresh directory for a separate analysis. |
| `--ldsc_results` | Required | Complete selected Python LDSC self/pair table, including heritabilities and intercepts; CSV/TSV/whitespace, optionally gzipped. |
| `--hm3` | Unset | Tab-separated file with SNP header; needed only for optional --write_munge_inputs. Preparation selects IDs without allele alignment. |
| `--gpca_input_folder` | Unset | Directory of existing nine-column GWAMA tables. Otherwise prepare from manifest VCFs for a full run. Not needed for ordinary validation-only runs. |

### VCF preparation and file layout

| Option | Default | What it changes |
| --- | --- | --- |
| `--splitby_chr` | `split` | `split`: all chromosomes 1–22 per trait. `nosplit`: one autosomal table per trait. Preparation and GWAMA must use matching layouts. |
| `--gpca_id_source` | `chr_pos_ref_alt` | `chr_pos_ref_alt` or `vcf_id` for prepared GWAMA SNPIDs. Choose compatible IDs across traits. |
| `--p_min` | `1e-300` | Preparation P floor, finite and strictly between 0 and 1. Does not change the ES/SE-derived GWAMA Z. |
| `--write_munge_inputs` | Off | Also prepare HapMap-selected raw tables, using NEF and no INFO column. Does not itself run munging. |
| `--munge_id_source` | `vcf_id` | `vcf_id` or `chr_pos_ref_alt` for optional raw-table SNP IDs; must match HapMap IDs. |
| `--prepare_workers` | `4` | Positive integer; workers for VCF preparation, separate from LDSC/GWAMA workers. |

### PCA method and quality checks

| Option | Default | What it changes |
| --- | --- | --- |
| `--rg_normalization` | `pair` | Python GPCA: `pair` uses original rg; `trait_wide` requires valid precomputed `rg_trait_wide` and `normalization_status`. CTI and original columns stay unchanged. |
| `--pca_matrix` | `correlation` | correlation uses selected rg with unit diagonal; covariance uses selected rg and self h2. CTI is unchanged. |
| `--heritability_scale` | `auto` | Python GPCA: `auto`, `observed`, `liability`, `mixed`. Selects existing h2 columns, never converts scales; covariance auto needs a common complete scale. |
| `--pc1_orientation` | `tutorial` | `tutorial` flips all PC1 loadings if their median is negative; `as_computed` keeps the eigenvector sign. |
| `--allow_missing_traits` | Off | GPCA: permit audited removal of manifest traits absent from LDSC results. Otherwise stop. |
| `--failed_ldsc_action` | `error` | GPCA: `error` or `drop_traits` for invalid/missing estimates. Dropping changes the analysed trait set; no imputation. |
| `--h2_z_warn_threshold` | `2` | GPCA: warn below this retained-trait h2/SE ratio. Nonnegative; 0 disables. Diagnostic only, never a removal rule. |
| `--rg_out_of_range_action` | `warn` | GPCA: `warn` or `error` for finite off-diagonal rg outside [−1,1]. Never clamps values. |
| `--negative_eigen_action` | `warn` | GPCA: `warn` or `error` for substantive negative PCA eigenvalues. No silent matrix repair. |
| `--matrix_eigen_tolerance` | `1e-8` | GPCA: positive relative tolerance separating substantive negative eigenvalues from floating-point noise. |
| `--duplicate_tolerance` | `0.001` | Python GPCA: positive absolute tolerance for duplicate rg, SE, P and intercept estimates. |
| `--duplicate_z_tolerance` | `0.01` | Python GPCA: positive absolute tolerance for duplicate-orientation Z values. |
| `--self_rg_tolerance` | `0.01` | Python GPCA: positive tolerance around self rg=1. SE=0 is accepted only with Z=+Inf and P=0; other checks still apply. |
| `--comparison_epsilon` | `1e-12` | Python GPCA: nonnegative floating-point allowance added to boundary comparisons. |
| `--z_consistency_tolerance` | `0.01` | Python GPCA: positive relative tolerance when comparing reported Z with rg/SE after duplicate collapse. |
| `--z_consistency_action` | `warn` | Python GPCA: `warn` or `error` when reported Z is inconsistent with rg/SE. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--validate_only` | Off | Stop after QC/PCA, skipping GWAMA preparation, GWAMA and export. Does not inspect per-SNP tables. |
| `--source_path` | Bundled modified GWAMA | Optional custom R source defining the expected modified GWAMA function and output interface. Unused in validation-only mode. |
| `--n_cores` | `0` (automatic) | Nonnegative integer. Automatic split-GWAMA workers use available physical cores, capped by 22 jobs; Windows is sequential. Export auto uses available logical/affinity cores up to 22. Nosplit GWAMA is one job. |
| `--bcftools` | `bcftools` on PATH | Executable name/path for VCF extraction or preparation; unnecessary when neither runs. |
| `--ldsc_chunk_size` | `250000` | Python GPCA: positive number of LDSC rows read per streaming chunk; controls memory/read size, not statistics. |
| `--color` | `auto` | Python GPCA help: `auto`, `always`, `never`. Auto respects NO_COLOR; no analysis effect. |

### Final export

| Option | Default | What it changes |
| --- | --- | --- |
| `--dataset_id` | Output directory name | Prefix for final combined/summary files and export audit. Pipeline uses the pipeline outdir name. |
| `--gwama_output_info` | Preserve reported INFO | Optional finite constant in [0,1] for the selected-column summary only. Required for bundled GWAMA because it has no INFO. Does not filter variants. |
| `--gwama_output_n_eff` | Preserve reported N_eff | Optional finite positive constant for the selected-column summary only. Does not change calculation N, weights, Z, P, BETA or SE. |
| `--gzip_level` | `1` | Integer 1–9 for final export only: 1 faster, 9 smaller. Uses pigz when installed; otherwise single-worker Python gzip. |
| `--archive_chromosomes` | Off | After successful export, move current-run source result/log files into outdir/chromosome_wise. |

### Help

| Option | Default | What it changes |
| --- | --- | --- |
| `--help` | Off | Show command usage, defaults and choices; exit without analysis. |
| `--prepare_help` | Off | Show shared preparation options and exit without analysis. |
| `--postprocess_help` | Off | Show automatic export options and exit without analysis. |
