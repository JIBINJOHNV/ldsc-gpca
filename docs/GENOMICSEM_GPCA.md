# PCA and PC1 GWAMA from GenomicSEM LDSC results

[README](../README.md) · [Inputs](INPUTS.md) · [GenomicSEM LDSC](GENOMICSEM_LDSC.md) · [Python GPCA](PYTHON_GPCA.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Files you need](#files-you-need)
- [Check results and calculate PCA](#check-results-and-calculate-pca)
- [What happens in order](#what-happens-in-order)
- [Run PC1 GWAMA](#run-pc1-gwama)
- [Choose the PCA method](#choose-the-pca-method)
- [Prepare GWAMA tables automatically](#prepare-gwama-tables-automatically)
- [QC and failure handling](#qc-and-failure-handling)
- [Results and export settings](#results-and-export-settings)
- [All options](#all-options)

## When to use this command

`ldsc-gpca genomicsem gpca` uses completed GenomicSEM LDSC results to
check the selected traits and calculate PCA. A full run then uses PC1 loadings
and the LDSC intercept matrix to combine per-SNP statistics with GWAMA.
It does not rerun LDSC. It requires `Rscript` on PATH.

To produce the GenomicSEM LDSC results first, use
[`ldsc-gpca ldsc --ldsc_backend genomicsem`](GENOMICSEM_LDSC.md).
`--ldsc_backend` applies to `ldsc` and `pipeline`; this GPCA command reads GenomicSEM
results and does not accept that selector.

## Files you need

Supply a trait-selection CSV, the final GenomicSEM RData and a fresh output folder.
The manifest needs only `traitname`, in the desired analysis order; use
[the four-trait example](INPUTS.md#selecting-traits-from-completed-ldsc-results)
or the `Selected_Traits.csv` produced by GenomicSEM LDSC.

`--ldsc_results` must load an object named `LDSCoutput` containing all five
matrices: `S,V,I,S_Stand,V_Stand`. Use `genomicsem_LDSC.RData` from the
[GenomicSEM LDSC command](GENOMICSEM_LDSC.md). Its first-pass `_raw.RData` is an
audit result, not the downstream GPCA input. A Python CSV is also not a valid
replacement for this GenomicSEM object.

A full GWAMA run additionally needs per-trait nine-column tables via
`--gpca_input_folder`, or suitable `vcf_files` in the manifest for automatic
preparation. See [filenames and allele requirements](INPUTS.md#gwama-tables).

## Check results and calculate PCA

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --outdir /results/genomicsem_pca_check \
  --validate_only
```

This checks GenomicSEM estimates and writes matrices, PC1 loadings and QC reports.
It does not prepare or inspect per-SNP GWAMA tables, run GWAMA or export SNP
results. Review retained traits and warnings before the full run.

## What happens in order

1. Load the GenomicSEM object and validate its matrices and trait names.
2. Select traits in manifest order, applying any explicitly requested missing
   or failed-estimate policies. Subset/reorder the real sampling covariance
   matrices consistently.
3. Use GenomicSEM `S_Stand` for correlation PCA (default), or GenomicSEM `S` for
   covariance PCA. Use GenomicSEM `I` for CTI in either case.
4. Check matrices, calculate PC1 with symmetric eigen decomposition and write
   QC, eigenvalue and loading reports.
5. Unless `--validate_only`, run GWAMA on the retained traits' per-SNP tables
   and automatically combine/export successful current-run results.

Automatic VCF preparation, when needed, occurs before R QC and trait removal.
Only PC1 goes into GWAMA even though eigenvalue reports describe all components.
Neither GenomicSEM nor Python GPCA runs paLDSC in this package.

## Run PC1 GWAMA

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/genomicsem_pc1 \
  --splitby_chr nosplit \
  --dataset_id four_traits_genomicsem \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Use the same layout as preparation: `nosplit` means one autosomal file per
trait; default `split` needs chromosomes 1–22 for each trait. Set `GWAMA_INFO`
to your justified downstream export constant before the command, because the
bundled GWAMA output has no INFO. This is not an imputation-quality estimate or
a variant filter. See [export metadata](INPUTS.md#export-metadata).

The bundled GWAMA uses available traits per SNP, compatible IDs and allele
swaps, weights based on `sqrt(N) × PC1 loading`, and CTI for correlated errors.
It does not convert genome builds or resolve strand flips. Export overrides do
not change GWAMA's inherited BETA, SE or N_eff calculation.

## Choose the PCA method

Correlation PCA uses GenomicSEM `S_Stand`, with traits on a standardized scale.
Covariance PCA uses GenomicSEM `S`, preserving its original trait scales:

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --outdir /results/genomicsem_covariance_check \
  --pca_matrix covariance \
  --pc1_orientation tutorial \
  --validate_only
```

Covariance PCA can give more influence to traits with larger genetic variance;
check scale compatibility before interpreting it. GenomicSEM GPCA has no Python
`--rg_normalization` or `--heritability_scale` setting. It uses the matrices
already estimated by GenomicSEM. CTI remains `I` for either PCA method.

The default sign convention flips all PC1 loadings together when their median
is negative. `--pc1_orientation as_computed` keeps the eigenvector sign.
Changing the whole-vector sign does not change association strength or P values.
To run GWAMA with a nondefault method, repeat those choices in your full command.

## Prepare GWAMA tables automatically

With suitable `vcf_files` in the manifest, omit `--gpca_input_folder`:

```bash
ldsc-gpca genomicsem gpca \
  --input /data/traits_vcf.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --outdir /results/genomicsem_pc1_auto_prepare \
  --splitby_chr nosplit \
  --gpca_id_source chr_pos_ref_alt \
  --p_min 1e-300 \
  --prepare_workers 4 \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

VCFs need FORMAT `AF,ES,SE,LP,NEF`. Their NEF supplies GWAMA N, separately from
N used in the completed GenomicSEM LDSC analysis. Preparation has no INFO threshold
or constant-N override. Do not request nondefault preparation options with an
existing `--gpca_input_folder` in this standalone command. See
[preparation details](PREPARE.md).

## QC and failure handling

The default is to stop for absent selected traits or invalid estimates.
`--allow_missing_traits` permits audited removal of absent traits.
`--failed_ldsc_action drop_traits` permits audited removal for failed estimates;
it does not impute missing values. At least two valid traits must remain.
These are downstream GPCA policies, separate from GenomicSEM LDSC's initial
`--invalid_h2_action drop` default.

Weak h2/SE is only a warning below 2 by default. Estimated correlations outside
`[-1,1]` and substantive negative genetic eigenvalues warn by default; their
options can require an error instead. Values are not clamped and matrices are
not silently repaired. CTI must be positive definite and PC1's eigenvalue must
be finite and positive. Any nonpositive-eigenvalue guard used in loadings is
reported. See [QC interpretation](REFERENCE.md#understand-qc-and-scale-choices).

Chromosome GWAMA jobs retry once before an exhausted failure stops analysis.
Workers default to automatic selection; Windows runs sequentially. Use fresh
output directories to avoid conflicting audit or completed export files.

After GWAMA writes its raw results, numerical QC checks every SNP before a
worker can report success. The run fails for zero available weight, nonfinite
BETA/Z/P, nonpositive or nonfinite SE/N, or P outside `[0,1]`. P underflow to
zero is allowed. Signed loadings and partial trait availability remain valid
when they produce a defined association; weights are not tested by their
signed sum. For example, a SNP available only in a trait with zero PC1 loading
has no available weight and fails even when its reported N_eff is positive.

Raw results remain available alongside `<output>.GWAMA_QC_Issues.csv` and
`<output>.GWAMA_QC_Summary.csv`. Issues identify result row, SNPID and reason;
summary files record the fixed `error` policy. Failed jobs remain failed in
`GWAMA_Run_Status.csv` and the existing worker-attempt audit. No SNP is silently
excluded and no Z is fabricated. The exporter independently checks numerical
values before overrides, combined output or archival, writing
`<output>.GWAMA_Export_QC_Issues.csv` and `..._Summary.csv`. An N_eff override
cannot hide an invalid N_eff already present in the raw output.

## Results and export settings

Inspect QC, `GWAMA_Run_Status.csv`, the matrices and
`GenomicPCA_PC1_Weights_Used.csv`. The full example also writes:

- `four_traits_genomicsem_GWAMA_combined_results.txt.gz`.
- `harmonisation_input/four_traits_genomicsem_GPCA_inputs.txt.gz`.
- `four_traits_genomicsem_postprocess.json`.

The first file retains original GWAMA columns; the second is a 12-column PC1
summary, not a nine-column input for another trait-level GWAMA analysis.

Optionally add `--gwama_output_n_eff "${GWAMA_N_EFF:?Set a justified summary N_eff}"`
to override only that summary's N_eff. Omit it to preserve reported values.
`--gzip_level` defaults to 1 for faster final export; 9 favours smaller files.
`--archive_chromosomes` moves current-run source results/logs after export.
See the [output catalog](REFERENCE.md#find-and-interpret-the-outputs).

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

Preparation settings apply only when VCF preparation runs. Export settings apply only after a full successful GWAMA run; both stages are skipped by `--validate_only`.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Output directory. Use a fresh directory for a separate analysis. |
| `--ldsc_results` | Required | Final RData containing LDSCoutput with S, V, I, S_Stand and V_Stand. |
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
| `--pca_matrix` | `correlation` | correlation uses GenomicSEM S_Stand; covariance uses GenomicSEM S on its original scales. CTI remains GenomicSEM I. |
| `--pc1_orientation` | `tutorial` | `tutorial` flips all PC1 loadings if their median is negative; `as_computed` keeps the eigenvector sign. |
| `--allow_missing_traits` | Off | GPCA: permit audited removal of manifest traits absent from LDSC results. Otherwise stop. |
| `--failed_ldsc_action` | `error` | GPCA: `error` or `drop_traits` for invalid/missing estimates. Dropping changes the analysed trait set; no imputation. |
| `--h2_z_warn_threshold` | `2` | GPCA: warn below this retained-trait h2/SE ratio. Nonnegative; 0 disables. Diagnostic only, never a removal rule. |
| `--rg_out_of_range_action` | `warn` | GPCA: `warn` or `error` for finite off-diagonal rg outside [−1,1]. Never clamps values. |
| `--negative_eigen_action` | `warn` | GPCA: `warn` or `error` for substantive negative PCA eigenvalues. No silent matrix repair. |
| `--matrix_eigen_tolerance` | `1e-8` | GPCA: positive relative tolerance separating substantive negative eigenvalues from floating-point noise. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--validate_only` | Off | Stop after QC/PCA, skipping GWAMA preparation, GWAMA and export. Does not inspect per-SNP tables. |
| `--source_path` | Bundled modified GWAMA | Optional custom R source defining the expected modified GWAMA function and output interface. Unused in validation-only mode. |
| `--n_cores` | `0` (automatic) | Nonnegative integer. Automatic split-GWAMA workers use available physical cores, capped by 22 jobs; Windows is sequential. Export auto uses available logical/affinity cores up to 22. Nosplit GWAMA is one job. |
| `--bcftools` | `bcftools` on PATH | Executable name/path for VCF extraction or preparation; unnecessary when neither runs. |

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
