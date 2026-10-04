# ldsc-gpca

Combine GWAS summary statistics from several traits into a SNP-level association
analysis of their first genetic principal component (PC1). Start from supported
GWAS files, or continue from completed LDSC estimates.

The package connects Python LDSC or native GenomicSEM LDSC to the
[Fürtjes genomicPCA procedure](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html)
and its modified GWAMA weighting. It reports all principal components, but
**uses only PC1 for SNP-level GWAMA**. The modified GWAMA function is bundled.

## Contents

- [Understand the workflow](#understand-the-workflow)
- [Install and activate](#install-and-activate)
- [Choose a workflow](#choose-a-workflow)
- [Run a complete analysis](#run-a-complete-analysis)
- [Use completed LDSC results](#use-completed-ldsc-results)
- [Prepare the right inputs](#prepare-the-right-inputs)
- [Choose analysis settings](#choose-analysis-settings)
- [Find and interpret results](#find-and-interpret-results)
- [Validation and reproducibility](#validation-and-reproducibility)
- [Troubleshooting](#troubleshooting)
- [Help, methods and attribution](#help-methods-and-attribution)

## Understand the workflow

A full analysis needs two data streams from the same traits:

```text
GWAS summary statistics
  ├─ Filter and munge → LDSC → genetic matrix → PCA → PC1 loadings
  │                         └→ intercept/error-covariance matrix (CTI)
  └─ Prepare per-trait SNP tables ────────────────────────────────┐
                   PC1 loadings + CTI + SNP tables → GWAMA → exports
```

*Munging* converts and filters summary statistics into a backend's LDSC input
format. LDSC estimates trait heritabilities, genetic relationships and
intercepts. PCA summarizes the genetic matrix; GWAMA combines per-SNP
statistics using PC1 loadings and the intercept matrix, which accounts for
correlation between trait-level errors, including sample overlap.

**Munged files, completed LDSC results and GWAMA inputs are different files.**
Completed LDSC results are enough for QC and PCA. SNP-level GWAMA also needs
per-trait tables containing alleles, frequency, N and association statistics.

## Install and activate

From a Bash terminal with Git and internet access:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
conda activate ldsc-gpca
ldsc-gpca --version
ldsc-gpca --help
```

Activate **`ldsc-gpca`** in each new terminal. The installer creates the main
Python/R environment and a separate `ldsc-gpca-ldsc` environment for Python
LDSC. Managed commands launch that isolated runtime automatically; keep the
main environment active. If Conda is absent, the script installs Miniforge
using a verified checksum; follow its printed shell-initialization instructions
before activating. This bootstrap also requires `curl`.

The supplied [main environment](environment.yml) requests Python 3.11, R 4.3,
bcftools and pigz. Package metadata allows Python **3.10–3.12**; this is distinct
from the installer's chosen version. The [LDSC environment](environment.ldsc.yml)
uses Python 3.10 and older pinned numerical dependencies. Setup pins
[CBIIT Python LDSC](https://github.com/CBIIT/ldsc/tree/6c673952cee74bd5c57aef1555a03b1c015399a0)
and [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM/tree/6b65ca5db39fdade08b0d811477be1cdd57b5039)
to source commits. The environment YAML files are not complete dependency locks.

The setup script recognizes Linux and macOS on x86_64 and ARM64; successful
installation still depends on Conda packages being available for that platform.
It does not support native Windows installation. It refuses to overwrite
existing environments. See [installation and recovery](docs/REFERENCE.md#install-and-activate)
for custom paths, activation problems and dependency updates.

For a **package-only update**, with the main environment active, enter the
checkout containing the revision you want and run:

```bash
python -m pip install .
ldsc-gpca --version
```

This installs the checkout's Python package and resolves its Python requirements;
it does not update the separate LDSC runtime or installed GenomicSEM/R software.
Keep the environments and code revision together when reproducing a run.

The [Docker guide](DOCKER.md) describes building a local image without host Conda.
Some examples there retain older option names and image tags: use this README's
current `--ldsc_results` and `--gwama_output_info` names, and rebuild from your
chosen checkout. The guide is not evidence of a published container image.

## Choose a workflow

Choose based on the files you already have and the result you need. If you
still need to run LDSC and want PC1 SNP association results, start with
**`pipeline`**. If LDSC is already finished, go straight to
[PCA/GWAMA](#continue-from-completed-ldsc-results).

The commands below start with `ldsc-gpca`. Each linked guide includes complete
commands, required manifest columns, options and defaults. Use at least two
traits for PCA/GWAMA; preparation alone can process one.

### Run the complete analysis

[`ldsc-gpca pipeline`](docs/PIPELINE.md) prepares inputs where needed, runs
LDSC, checks the estimates, calculates PCA, runs PC1 GWAMA and exports SNP
results. Supply a trait manifest and LD references, then choose the LDSC input
route:

| Your LDSC starting files | Pipeline setting | Guide |
| --- | --- | --- |
| GWAS summary-statistics VCFs with the required fields | Python LDSC is the default. List VCF paths in `vcf_files`. | [Python from VCFs](docs/PIPELINE.md#python-from-vcfs) |
| Existing Python LDSC munged files | Add **both** `--ldsc_only` and `--munged_dir DIR`. | [Python from munged files](docs/PIPELINE.md#python-from-munged-files) |
| Raw GWAS summary-statistics tables for GenomicSEM | Add `--ldsc_backend genomicsem`; list table paths in `sumstats_file`. | [GenomicSEM from raw tables](docs/PIPELINE.md#native-genomicsem-from-raw-tables) |
| GWAS summary-statistics VCFs for GenomicSEM | Add `--ldsc_backend genomicsem --vcf_input`; list VCF paths in `vcf_files`. Binary traits need an appropriate manifest N and both prevalences. | [GenomicSEM from VCFs](docs/PIPELINE.md#native-genomicsem-with-explicit-vcf-input) |
| Existing GenomicSEM-compatible munged files | Add `--ldsc_backend genomicsem`, then either `--munged_dir DIR` or `--munged_input` with manifest `munged_file` paths. | [GenomicSEM from munged files](docs/PIPELINE.md#native-genomicsem-from-munged-files) |

For GWAMA, also supply per-trait SNP tables through `--gpca_input_folder`, or
include VCF paths in `vcf_files` so the pipeline can prepare them. The same VCFs
can serve both stages if they contain [all required fields](#vcf-fields-and-sample-size).
**Munged LDSC files alone are not enough for SNP-level GWAMA.**

The explicit native VCF route retains FORMAT/SI as INFO for native filtering.
Its [sample-size rules](docs/GENOMICSEM_LDSC.md#start-from-gwas-vcfs) differ from
Python's binary extraction. The older [quantitative VCF fallback](docs/PIPELINE.md#native-genomicsem-from-quantitative-vcfs)
remains available without `--vcf_input`; it uses NEF and exports no INFO.

A full run saves LDSC estimates in `ldsc/`, and PCA reports, PC1 weights and
combined SNP results in `gpca/`. Choose a new output directory; pipeline requires
one that does not already exist. See [output files](#find-and-interpret-results).
To stop after LDSC and PCA, add `--validate_only`: LDSC still runs, while GWAMA
and final SNP export are skipped. Preparation runs only when needed for LDSC
inputs or explicitly requested.

### Run LDSC only

Use these commands when you want heritabilities, genetic correlations and
intercepts first. Each needs a trait manifest and LD references; routes that
perform munging also need the HapMap allele reference through `--hm3`.

| Your starting files | Command and input mode | What it does and produces |
| --- | --- | --- |
| GWAS summary-statistics VCFs | [`ldsc-gpca ldsc`](docs/PYTHON_LDSC.md#start-from-vcf-files) | Extracts and filters VCF records, munges them, runs Python LDSC and writes `ldsc_results.csv` plus QC reports. |
| Existing Python LDSC `{traitname}.sumstats.gz` files | [`ldsc-gpca ldsc --ldsc_only`](docs/PYTHON_LDSC.md#start-from-munged-files); use `--munged_dir DIR` to select their directory. | Skips extraction/munging, runs Python LDSC and writes the same final result format. |
| GWAS summary-statistics VCFs for native LDSC | [`ldsc-gpca genomicsem ldsc --vcf_input`](docs/GENOMICSEM_LDSC.md#start-from-gwas-vcfs); list paths in `vcf_files`. | Converts VCFs to raw tables with INFO, runs native munging/LDSC and writes `genomicPCA_LDSC.RData`, trait reports and VCF QC. |
| Raw GWAS text tables | [`ldsc-gpca genomicsem ldsc`](docs/GENOMICSEM_LDSC.md#basic-command-from-raw-tables); list paths in `sumstats_file`. | Munges the tables, runs native LDSC with trait QC and writes the final `genomicPCA_LDSC.RData` and `Selected_Traits.csv`. |
| Existing GenomicSEM-compatible munged tables | [`ldsc-gpca genomicsem ldsc`](docs/GENOMICSEM_LDSC.md#reuse-munged-files) with `--munged_dir DIR`, or `--munged_input` and manifest `munged_file` paths. | Skips munging, runs native LDSC with trait QC and writes the same final RData and selected-trait manifest. |

Select the input mode explicitly; the program does not choose raw versus munged
mode from the filename. In standalone Python `ldsc`, `--munged_dir` without
`--ldsc_only` changes where **new** munged files are written. None of these LDSC
commands runs PCA or GWAMA.

The Python LDSC guide explains the [final columns with a real-data example](docs/PYTHON_LDSC.md#results-and-next-step)
and links to the [folder, intermediate-file and audit reference](docs/PYTHON_LDSC_OUTPUTS.md).

### Continue from completed LDSC results

Choose the command that matches the LDSC results you have. Both read those
estimates, check the selected traits and calculate PCA without rerunning LDSC.

| Completed result | Command | What else you need |
| --- | --- | --- |
| Python `ldsc_results.csv`, containing all required self-pairs and trait pairs | [`ldsc-gpca gpca`](docs/PYTHON_GPCA.md) | A trait-selection manifest; for GWAMA, per-trait SNP tables or VCFs for automatic preparation. |
| Final native `genomicPCA_LDSC.RData`, containing `LDSCoutput` | [`ldsc-gpca genomicsem gpca`](docs/GENOMICSEM_GPCA.md) | A trait-selection manifest; for GWAMA, per-trait SNP tables or VCFs for automatic preparation. |

Pass the completed result with `--ldsc_results`. For a full GWAMA run, use
`--gpca_input_folder` for prepared tables, or include `vcf_files` in the manifest
and omit that option to prepare them automatically. The command saves matrices,
PC1 weights, QC reports and SNP results. Add `--validate_only` if you want only
QC and PCA; this skips preparation, GWAMA and SNP export. Validation-only runs
do not check per-SNP GWAMA input tables.

### Prepare GWAMA inputs separately

Use [`ldsc-gpca prepare`](docs/PREPARE.md) when you have GWAS VCFs and want to
make the per-trait SNP tables for a later analysis. Supply a manifest with
`traitname,vcf_files`. It writes nine-column tables under `gpca_inputs/` and QC
reports. Preparation uses its own [VCF field requirements](docs/PREPARE.md#files-you-need)
and does not run LDSC, PCA or GWAMA. With `--write_munge_inputs` and `--hm3`, it
can also write raw tables for later munging; these are not yet munged LDSC files.

For direct access to the upstream Python scripts, see the
[raw LDSC commands](docs/REFERENCE.md#raw-ldsc-commands). These bypass the
package's managed workflow and result collection.

## Run a complete analysis

### Python LDSC from GWAS VCFs

This example uses four quantitative traits. **The data and references are not
bundled.** Replace `/data`, `/references` and `/results` paths with your own.
Before running, supply:

1. One GWAS summary-statistics VCF per trait, satisfying **both** VCF field sets
   [below](#vcf-fields-and-sample-size). These are not genotype VCFs.
2. An ancestry-appropriate LD-score directory containing
   `1.l2.ldscore.gz` through `22.l2.ldscore.gz` and their `*.l2.M_5_50` files.
   By default it also supplies regression weights; `--ld_weights` can select a
   separate chromosome LD-score directory for weights.
3. A whitespace-delimited HapMap allele table with `SNP,A1,A2` headers for
   `--hm3`. Its IDs must match the summary statistics. Genome build and allele
   conventions must also agree with your chosen references.

Save `/data/traits_vcf.csv`:

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

The empty prevalence cells indicate quantitative traits. The Python pipeline
adds `ref=yes` for every trait when that column is absent, ensuring complete
self/pair coverage.

The bundled GWAMA function does not output INFO. Full export therefore requires
an explicit constant for its selected-column summary. Set the shell variable
`GWAMA_INFO` to a value between 0 and 1 justified by your downstream policy.
**It is export metadata, not measured imputation quality or an INFO filter.**
The guard below stops if the variable is unset.

```bash
ldsc-gpca pipeline \
  --input /data/traits_vcf.csv \
  --outdir /results/four_traits_python \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --dataset_id four_traits \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

`nosplit` uses one autosomal GWAMA table per trait. The default `split` requires
retained variants on every chromosome 1–22 for every trait. The pipeline output
directory must not already exist.

First inspect `/results/four_traits_python/Pipeline_Run_Status.json`, then
`ldsc/LDSC_Trait_Status.csv`, `gpca/Python_LDSC_Input_Validation_Summary.csv`
and `gpca/GenomicPCA_PC1_QC.csv`. Check the retained traits and warnings, then
review `gpca/GenomicPCA_PC1_Weights_Used.csv`. A full successful run also creates
`gpca/four_traits_GWAMA_combined_results.txt.gz` and
`gpca/harmonisation_input/four_traits_GPCA_inputs.txt.gz`.

To run LDSC and PCA only, replace the INFO line with `--validate_only`.
For optional INFO/MAF filters, chi-square cutoffs, correlation methods, worker
counts and export settings, see [analysis settings](#choose-analysis-settings)
and the [pipeline recipes](docs/PIPELINE.md#change-filters-pca-and-export-settings).

### Native GenomicSEM from raw tables

Native LDSC has a different input contract. For each trait, supply a
whitespace-delimited raw table with SNP IDs, A1/A2, valid P values, a signed
effect such as BETA or Z, and suitable N. Recognized INFO/frequency columns
allow native munging filters. The main environment must have the pinned
GenomicSEM installation from setup.

Save `/data/traits_native_raw.csv`:

```csv
traitname,sumstats_file,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A_raw.tsv,,
Trait_B,/data/Trait_B_raw.tsv,,
Trait_C,/data/Trait_C_raw.tsv,,
Trait_D,/data/Trait_D_raw.tsv,,
```

Supply GWAMA tables separately. If the same GWAS data are available as
preparation-compatible VCFs in the earlier manifest, create those tables with:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared \
  --splitby_chr nosplit
```

Then run the native pipeline using those tables:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/traits_native_raw.csv \
  --outdir /results/four_traits_native \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --info_filter 0.9 \
  --maf_filter 0.01 \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --splitby_chr nosplit \
  --dataset_id four_traits_native \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Here the INFO/MAF values are native defaults and apply only when the raw tables
contain recognized fields. Omitted `--chisq_max` selects native automatic
filtering. An optional positive `N` column in the manifest **replaces** raw-file
N during native munging. For binary traits, supply both sample and population
prevalence; this route does not infer the missing value.

For reuse, explicitly select `--munged_dir` or `--munged_input` and omit `--hm3`;
raw munging filters and manifest N are not reapplied. To start from VCFs, select
`--vcf_input` on `genomicsem ldsc`, or combine it with `--ldsc_backend genomicsem`
on `pipeline`. This conversion retains INFO and uses quantitative NEF unless
manifest N overrides it; binary VCFs require an appropriate explicit N and both
prevalences. See the [native VCF guide](docs/GENOMICSEM_LDSC.md#start-from-gwas-vcfs).

## Use completed LDSC results

Save `/data/selected_traits.csv` in the order you want to analyze:

```csv
traitname
Trait_A
Trait_B
Trait_C
Trait_D
```

For Python estimates, check the results and calculate PCA:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/four_traits_python/ldsc/ldsc_results.csv \
  --outdir /results/python_pca_check \
  --validate_only
```

For native estimates, use the final RData object:

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/four_traits_native/ldsc/genomicPCA_LDSC.RData \
  --outdir /results/native_pca_check \
  --validate_only
```

Review retained traits, matrices, warnings and loadings. If LDSC removed traits,
use its retained-trait manifest instead of requesting traits that are absent.
These commands do **not** inspect SNP tables. To continue from Python estimates
to GWAMA, supply prepared inputs and use a new output directory:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/four_traits_python/ldsc/ldsc_results.csv \
  --gpca_input_folder /results/four_traits_python/prepare/gpca_inputs \
  --outdir /results/python_pc1_gwama \
  --splitby_chr nosplit \
  --dataset_id four_traits \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

This recalculates QC/PCA from the saved estimates and runs GWAMA; it does not
rerun LDSC. Repeat any nondefault PCA/QC options used in the check. For native
results, use `genomicsem gpca` with the final RData file and your prepared folder;
see the [native GWAMA example](docs/GENOMICSEM_GPCA.md#run-pc1-gwama).
If a previous run used `--validate_only`, GWAMA inputs may still need preparation.

## Prepare the right inputs

Manifests are comma-separated files with unique, nonempty `traitname` values.
Use simple names such as `Trait_A`, without spaces or path characters. **Manifest
row order defines trait order** in matrices, loadings and GWAMA inputs; check
retained-trait reports after any removal. Absolute file paths work consistently.
Relative paths inside manifests resolve beside the manifest, except standalone
Python `ldsc` VCF paths, which resolve from the working directory.

| Input | Required content or mode |
| --- | --- |
| Pipeline manifest | `traitname,population_prevalence,sample_prevalence`, plus route-specific paths: `vcf_files`, native `sumstats_file`, or explicit munged reuse. For full GWAMA, supply prepared inputs or preparation-compatible VCFs. Python `ref`, if present, must be `yes` for every trait. |
| Standalone Python LDSC manifest | `traitname,ref,population_prevalence,sample_prevalence`; add `vcf_files` unless using `--ldsc_only`. Set all `ref=yes` for complete GPCA coverage. Reuse expects `{traitname}.sumstats.gz`. |
| Native LDSC manifest | `traitname,population_prevalence,sample_prevalence`; default raw mode adds `sumstats_file`. `--vcf_input` uses `vcf_files`. `--munged_dir` expects one `{traitname}.sumstats` or `.sumstats.gz`; `--munged_input` uses `munged_file` paths. Choose one mode. |
| GPCA selection | `traitname` plus completed LDSC results supplied through `--ldsc_results`. Add `vcf_files` only for automatic GWAMA preparation. |
| Python LDSC results | Headered CSV/TSV/whitespace table, optionally gzip-compressed: `p1,p2,rg,se,z,p,h2_int,h2_int_se,gcov_int,gcov_int_se`, plus at least one h2/SE pair: `h2_obs,h2_obs_se` or `h2_liab,h2_liab_se`. See [scale and pair requirements](docs/REFERENCE.md#python-ldsc-results-table). |
| Native LDSC results | Final RData containing `LDSCoutput` with `S,I,V,S_Stand,V_Stand`. The intermediate `genomicPCA_LDSC_raw.RData` is not the final GPCA input. |
| Munged LDSC files | Per-trait `SNP,A1,A2,N,Z` columns; native reuse requires tab separation. They are inputs to regression, not its outputs. Preserve Python munging provenance sidecars. |

### VCF fields and sample size

Each supported VCF represents **one GWAS sample**, not individual genotypes.
A `.vcf.gz` extension does not establish compatibility.

| Stage | Required VCF fields beyond coordinates, REF and ALT | Values used |
| --- | --- | --- |
| GWAMA preparation | FORMAT `AF,ES,SE,LP,NEF` | EA=ALT, OA=REF, EAF=AF, N=NEF, Z=ES/SE, P from LP. |
| Python LDSC extraction | IDs; INFO `AF,EUR`; FORMAT `SI,AF,EZ,LP,NEF` | SI for INFO filtering, FORMAT/AF for MAF, INFO/AF versus INFO/EUR for frequency difference, EZ/LP for association statistics. |
| Python extraction with population prevalence | Python fields above plus FORMAT `NC,NCO` | N=NC+NCO; missing sample prevalence can be inferred from median case fraction. |
| Native LDSC with `--vcf_input` | IDs; FORMAT `AF,ES,SE,LP,SI`; NEF unless manifest N is supplied | INFO=SI; MAF from AF; P from LP and effect direction from ES. Quantitative N=NEF by default; binary traits require an appropriate explicit manifest N and both prevalences. |

A Python VCF pipeline needs the union of the first two field sets, plus the
binary fields when applicable. Its LDSC filtering does not also filter GWAMA
inputs. Preparation retains valid autosomal rows; its default P floor is
`1e-300` and does not change the Z used by GWAMA.

There is no universal sample-size override. Python extraction without population
prevalence uses NEF. Native raw munging uses raw-file N unless manifest N replaces
it; munged reuse uses the N already in the files. Prepared GWAMA tables use NEF
independently of Python's binary LDSC rule. Quantitative traits leave both
prevalences empty. For binary traits, sample prevalence means the fraction of
cases in the GWAS, and population prevalence controls liability scaling; supplied
values must be strictly between 0 and 1. Native LDSC requires both or neither.
Read the [sample-size rules](docs/INPUTS.md#sample-size-and-prevalence) before
combining binary traits or reusing effective-N inputs.

### Per-trait GWAMA tables and allele handling

Preparation writes tab-separated tables with these nine columns:

```text
SNPID CHR BP EA OA EAF N Z P
```

The reader also accepts extra/reordered columns and aliases `A1`→`EA`, `A2`→`OA`
and `p`→`P`. Use the exact filenames:

- `split`: `{traitname}_chr{CHR}_GenomicPCA_inputs.tsv`, chromosomes 1–22.
- `nosplit`: `{traitname}_GenomicPCA_inputs.tsv`.

Default prepared SNPIDs are `CHR_POS_REF_ALT`. All traits must use compatible
IDs, genome builds and allele conventions. The bundled GWAMA uses the union of
SNPIDs and available traits at each SNP. For a shared ID, it can align an allele
swap by negating Z and replacing EAF with `1 − EAF`; unresolved allele mismatches
are excluded. It does not resolve strand flips or perform liftover. Different
IDs do not match merely
because their coordinates agree. Output reformatting adds no further allele
harmonization. See [full input schemas](docs/INPUTS.md) and
[preparation options](docs/PREPARE.md#choose-ids-layout-and-p-value-floor).

## Choose analysis settings

### Filters and backend differences

Python LDSC and native GenomicSEM are distinct estimation routes. Matching
inputs and normalization does not make their regressions or results identical.
Use the same appropriate ancestry/build references and document backend,
munging, N/prevalence conventions and filtering when comparing them.

| Setting | Python LDSC | Native GenomicSEM LDSC |
| --- | --- | --- |
| Raw input INFO/MAF | `--info_min 0.7`, `--maf_min 0.01`; later munging uses `--munge_maf_min 0.005`. Extraction also defaults to `--max_af_difference 0.2`. | `--info_filter 0.9`, `--maf_filter 0.01`, when recognized raw columns exist. |
| Large chi-square SNPs | Default: no additional managed cutoff. `--chisq_max auto` or a positive integer enables a per-trait filter on munged Z² before regression. | Omission selects native automatic filtering; `--chisq_max` accepts a positive number, **not** the string `auto`. Applied inside native LDSC after merging with LD scores/weights. |
| Automatic cutoff | `max(80, 0.001 × max(N))`, using complete rows matched to the LD-score/weight references. Filtered copies leave source munged files intact. | The same cutoff formula is computed for each trait's merged data by the pinned native implementation. |
| Invalid estimates | Managed collection defaults to `--result_failure_action error`; `report` and `drop_traits` have different downstream consequences. | A first native LDSC pass identifies invalid nonpositive/nonfinite h2; `--invalid_h2_action drop` is the default, with `error` available. A subsequent pass computes standardized results. |

The Python cutoff is a managed per-trait filter, not a passthrough of the raw
upstream bivariate `--chisq-max` setting. Reusing munged data skips raw munging
filters; regression-stage cutoff choices still apply. See
[Python filtering](docs/PYTHON_LDSC.md#limit-large-chi-square-values),
[native filtering](docs/GENOMICSEM_LDSC.md#filters-sample-size-and-prevalence) and
[the pinned native implementation](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R).

### Genetic matrices, PC1 and GWAMA weights

The **genetic matrix** determines the PCA. The **CTI intercept matrix** determines
error correlation in GWAMA; it is not the genetic matrix.

| Choice | Behavior |
| --- | --- |
| `--pca_matrix correlation` (default) | Uses `S_Stand`: diagonal 1 and off-diagonal genetic correlations. Python builds it from the selected rg column; native GPCA uses native `S_Stand`. |
| Python `--rg_normalization pair` (default) | Uses the original Python `rg`, normalized using each pair's fitted heritabilities. |
| Python `--rg_normalization trait_wide` | Uses the extra `rg_trait_wide` column: native pair covariance divided by `sqrt(h2_A_self × h2_B_self)` on the unconverted scale. Self-pairs set the diagonal to 1. Compilation adds this without changing original estimates or rerunning regression. |
| `--pca_matrix covariance` | Python reconstructs covariance from selected rg and self-pair h2; `--heritability_scale` controls eligible h2 columns. Native GPCA uses `S`. Units and trait scales affect this PCA. |
| CTI | Python: self-pair `h2_int` on the diagonal, pairwise `gcov_int` off diagonal. Native: `I`. Correlation normalization does not replace CTI. |

Trait-wide selection requires valid calculated normalization values for all
selected pairs. Its original SE/Z/P columns remain diagnostics of the original
rg, not uncertainty estimates for the new ratio. “Unconverted” does not certify
that effective-N binary estimates have a conventional population observed-scale
interpretation. Older exports have a reconstruction fallback with additional
requirements; see [normalization details](docs/REFERENCE.md#optional-trait-wide-correlation-normalization).
Choosing h2 columns does not itself convert their scale. Native GPCA has no
Python-style normalization or heritability-scale selector.

To inspect optional trait-wide normalization and covariance PCA together:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/four_traits_python/ldsc/ldsc_results.csv \
  --outdir /results/python_trait_wide_covariance_check \
  --rg_normalization trait_wide \
  --pca_matrix covariance \
  --heritability_scale auto \
  --validate_only
```

PCA uses symmetric eigen decomposition. The PC1 loading vector is
`v1 × sqrt(lambda1)`; GWAMA's per-SNP weight for trait t is
`sqrt(N_t) × loading_t`, with CTI in the Z-score variance calculation. The
bundled function's argument named `h2` contains these **loadings**, not SNP
heritabilities. `--pc1_orientation tutorial` is the default and flips the
whole vector if its median is negative; `as_computed` preserves the eigenvector
sign. A global sign flip reverses directional interpretation, so retain the
reported orientation when comparing analyses.

The reports describe a PCA of the selected genetic matrix, not confirmatory
factor analysis or causal trait contributions. Covariance loadings depend on
units. Only PC1 is used for GWAMA; other components are reported for inspection.

### Checks that affect whether analysis continues

- Python GPCA requires every self-pair and unordered pair: `k × (k + 1) / 2`
  unique pairs, or **10 for four traits**. Triangles and full symmetric tables
  are accepted; consistent duplicates collapse and conflicting duplicates fail.
  A raw `ldsc.py --rg A,B,C` invocation alone does not produce B–C.
- Required selected values must be usable, selected self-pair h2 positive,
  ordinary required SEs positive and P values in `[0,1]`. The narrow self-pair
  exception accepts `rg` near 1 with **`se=0, z=+Inf, p=0` together**, preserving
  and auditing those values. Other checks still apply; zero between-trait SE
  is invalid. Default self-rg tolerance is `0.01`.
- Pairwise P, Z and SE are diagnostics. They do not significance-filter pairs
  or weight the PCA. Low h2/SE (default warning threshold 2) prompts review.
  Missing/failed selected traits stop GPCA by default; explicit removal options
  record what was dropped and require at least two retained traits.
- Out-of-range rg and negative genetic-matrix eigenvalues warn by default;
  corresponding action options can make them fatal. Values are not silently
  clamped or replaced by a nearest-positive-definite matrix. The loading
  calculation guards square roots of nonpositive eigenvalues and reports it;
  PC1 itself must have a finite positive eigenvalue.
- Matrices and trait ordering must agree. CTI must be positive definite for
  GWAMA. A PC1 QC pass does not prove the full GWAMA/export run succeeded.

Python long-format results supply marginal SEs, but not GenomicSEM's full
cross-estimate sampling covariance matrices `V` and `V_Stand`. Do not use that
input for `paLDSC` or other procedures requiring those matrices, or fabricate
diagonal substitutes. This limitation does not block the PC1/GWAMA calculation.
See [QC and scale choices](docs/REFERENCE.md#understand-qc-and-scale-choices).

### Export settings are not analysis corrections

`--gwama_output_info` supplies a constant in `[0,1]` for the selected-column
summary. `--gwama_output_n_eff` optionally replaces its N_eff with a positive
constant. **Neither changes analysis N, loadings, Z, P, BETA or SE**, and the
N_eff override does not change the full combined GWAMA result. Omit it to
retain GWAMA's reported N_eff. These options do not recalibrate the bundled
BETA/SE/N_eff formulas or validate their interpretation.

See [export options](docs/INPUTS.md#export-metadata) for examples and
[worker/thread settings](docs/REFERENCE.md#cpu-workers-and-numerical-threads)
for parallelism. Worker counts and numerical-library threads are separate;
more workers are not always faster. The guides list all command-specific
options, including INFO filters, N inputs, retries, compression and archiving.

## Find and interpret results

Pipeline paths below are relative to `--outdir`. For standalone stages,
outputs are directly under that stage's `--outdir`. Check status and trait
removals before interpreting loadings or SNP associations.

| Location / file | What to inspect |
| --- | --- |
| `Pipeline_Run_Status.json`; `manifests/` | Overall completion, stage commands, failures and resolved trait selections. |
| `prepare/Preparation_Status.csv`; `prepare/GPCA_Input_QC_Summary.csv` | Preparation outcomes and retained/removed variants. Prepared inputs are in `prepare/gpca_inputs/`. |
| Python `ldsc/ldsc_results.csv`; `LDSC_Compilation_Status.csv`, `LDSC_Pair_Status.csv`, `LDSC_Trait_Status.csv` in `ldsc/` | Final pair/self estimates and collection/QC outcomes. [Columns and real example](docs/PYTHON_LDSC.md#results-and-next-step). `drop_traits` also writes `LDSC_Retained_Traits.csv` and `LDSC_Dropped_Traits.csv`. |
| Native `ldsc/genomicPCA_LDSC.RData`; `Selected_Traits.csv`, `GenomicSEM_LDSC_Trait_QC.csv` in `ldsc/` | Final native estimates and traits retained after LDSC h2 QC. |
| Python `gpca/Python_LDSC_Input_Validation_Summary.csv`; `Python_LDSC_Dropped_Missing_Traits.csv`, `Python_LDSC_Dropped_Failed_Traits.csv` in `gpca/` | Input validation and any explicit GPCA trait removals. |
| Native `gpca/GenomicSEM_QC_Events.csv`; `GenomicSEM_QC_Removed_Traits.csv`, `GenomicSEM_Retained_Manifest.csv` in `gpca/` | Native GPCA QC and retained trait order. |
| `gpca/GenomicPCA_Correlation_Matrix_Used.csv`, `GenomicPCA_PCA_Matrix_Used.csv`, `GenomicPCA_CTI_Used.csv` | Exact correlation, selected PCA and intercept matrices used. Covariance mode also reports its covariance matrix. |
| `gpca/GenomicPCA_PC1_Weights_Used.csv`; `GenomicPCA_PC1_QC.csv` | Trait loadings, their order and numerical PC1 checks. `PASS` here is not whole-run success. |
| `gpca/GenomicPCA_Selected_Traits_Eigenvalues.csv`; `GenomicPCA_All_PCs_Variance_Explained.csv` | Eigenvalues and variance shares of the selected genetic matrix, not phenotypic variance explained or SNP heritability. Review negative eigenvalues before interpretation. |
| `gpca/GWAMA_Run_Status.csv` | Success/failure for chromosome or whole-genome GWAMA jobs; validation-only runs do not run these jobs. |
| `gpca/{dataset_id}_GWAMA_combined_results.txt.gz` | Full GWAMA columns, sorted by chromosome/position, plus `count_question,count_plus,count_minus` from `Direction`. |
| `gpca/harmonisation_input/{dataset_id}_GPCA_inputs.txt.gz`; `gpca/{dataset_id}_postprocess.json` | Selected-column summary and export audit: sources, row counts, overrides and timings. |

The final selected-column summary has **12 tab-separated columns**:

```text
SNPID CHR BP EA OA EAF N_eff BETA SE Z PVAL INFO
```

It is a PC1 result, not the nine-column per-trait GWAMA input. The directory
name `harmonisation_input` does not mean it has already been harmonized with
another dataset. Export checks structural requirements such as IDs, positions
and Direction; it does not comprehensively validate every BETA, SE, PVAL or
preserved N_eff/INFO value. See the [full output catalog](docs/REFERENCE.md#find-and-interpret-the-outputs).

## Validation and reproducibility

The committed [pipeline validation report](tests/PIPELINE_VALIDATION.md)
identifies its v0.7.0 work and baseline commit. Its evidence has distinct scopes:

| Evidence | What it establishes | Boundary |
| --- | --- | --- |
| Unit/interface tests | Parsing, input contracts, failure handling and selected numerical invariants | Many pipeline tests mock stage execution. Passing them is not a scientific comparison of backends. |
| Injected native-estimator fixture | Two-pass orchestration, cutoff forwarding and retained-trait order | Does not run native LDSC regressions. |
| Saved-estimate replay of two real-data subsets | Pipeline handoffs reproduce separate-stage matrices, PC1 and GWAMA/export results using the same saved estimates | LDSC outputs were supplied, not newly fitted. SNP subsets were reconstructed fixtures; this is not a full raw-VCF or whole-genome analysis. |
| Fresh regression/integration checks | Require the actual pinned runtimes and appropriate inputs | Optional Python LDSC integration tests were skipped in the reported environment; native GenomicSEM was unavailable there. |

The replay does not establish cross-backend equivalence, BETA/SE or N_eff
calibration, or whole-genome performance. Its private-data artifacts are not
bundled as a public reproducible dataset. No such claims should be inferred
from pipeline completion or PC1 QC.

Retain the code commit, `ldsc-gpca --version`, commands, manifests, input/reference
checksums, munging provenance, logs and all QC/status files. Setup writes Conda
explicit specifications, a pip freeze and R session information in its provenance
location; native runs also write session information. Keep these with your
outputs because package-only updates do not synchronize external runtimes.

## Troubleshooting

| Symptom | What to check or do |
| --- | --- |
| Command, R package or LDSC runtime not found | Activate the main environment. Check setup completion and the configured isolated LDSC prefix; follow [environment recovery](docs/REFERENCE.md#activation-custom-locations-and-recovery). |
| Missing VCF/table fields | Match the fields to the actual stage above. Genotype VCFs, Python extraction VCFs, native raw tables and GWAMA inputs are not interchangeable. |
| Trait absent, order mismatch or incomplete pairs | Match exact trait names, use the retained manifest after removal, and obtain all self/unordered pairs. Reordering or dropping rows cannot supply missing estimates. |
| Split input files missing | Use the matching layout in preparation and GPCA. Default split requires every chromosome 1–22 for each trait; use `nosplit` for one autosomal file per trait. |
| Invalid h2, genetic matrix or CTI | Inspect pair/trait QC and scale/N choices. Do not fill missing estimates, clamp rg or substitute a repaired matrix to bypass the error. Use documented removal/error policies only with justification. |
| Pipeline output directory already exists | Use a new directory. Pipeline has no automatic resume mode; preserve the failed run and reuse completed stages explicitly. |
| Interrupted Python LDSC | Standalone `ldsc --restart` can reuse batches whose input/reference hashes, settings and runtime/output provenance still match. Failed or changed batches rerun. See [restart rules](docs/PYTHON_LDSC.md#reuse-munged-files-or-restart). |
| Worker failed | Python LDSC command failures and native munging worker failures normally get one retry (two attempts). Invalid numerical estimates follow the QC policy instead. After exhausted retries, inspect the reported failure and worker/status logs. |
| GWAMA completed but export failed | Check `GWAMA_Run_Status.csv`, INFO requirements, SNPIDs, Direction and existing export paths. Preserve successful outputs and read the [export/recovery details](docs/REFERENCE.md#gwama-results-and-final-exports) before rerunning; export refuses existing target files. |

Native standalone LDSC expects a fresh or empty output directory; preparation
also protects existing named outputs. A successful completed LDSC stage can be
used with a standalone GPCA command, avoiding another regression run. There is
no promise that rerunning a whole pipeline will skip completed work.

## Help, methods and attribution

```bash
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca genomicsem ldsc --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
```

Managed flags use underscores. Enable switches with the flag alone, such as
`--validate_only`. Pipeline and both GPCA commands also provide `--prepare_help`
and `--postprocess_help`. For complete option/default tables, use the workflow
guides linked in [Choose a workflow](#choose-a-workflow).

- [Full reference](docs/REFERENCE.md), [input schemas](docs/INPUTS.md) and
  [naming conventions](docs/NAMING.md).
- [Docker build/run guide](DOCKER.md) and [Nextflow image/QC examples](examples/nextflow/README.md).
  The Nextflow examples demonstrate tool checks and GPCA QC, not the full
  analysis pipeline; reconcile their local image tags with the image you build.
- [Validation report](tests/PIPELINE_VALIDATION.md) and
  [upstream source attribution](docs/REFERENCE.md#sources-and-license).
- [GitHub Issues](https://github.com/JIBINJOHNV/ldsc-gpca/issues) for reproducible
  problems. Include the command, error and versions without private data.

The scientific procedure follows [Anna Fürtjes' genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
The bundled modified GWAMA source credits Hill Ip and Bart Baselmans; see its
[provenance](src/ldsc_gpca/r/vendor/README.md). Python LDSC and GenomicSEM retain
their own methods and attribution. Cite the methods actually used in your
analysis. This repository does not declare a project-wide license; upstream
components retain their own terms. See [sources and license](docs/REFERENCE.md#sources-and-license)
before redistribution.
