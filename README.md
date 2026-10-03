# ldsc-gpca

Combine GWAS summary statistics from several traits into a SNP-level association
analysis of their first genetic principal component (PC1). Run the whole analysis
with one command, or continue from LDSC results you already have.

The package supports Python LDSC and native GenomicSEM LDSC and follows the
[Fürtjes genomicPCA procedure](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
The modified GWAMA function is included. PCA reports describe all components;
**only PC1 is used for GWAMA**.

## Contents

- [Understand the workflow](#understand-the-workflow)
- [Install and activate](#install-and-activate)
- [Choose where to start](#choose-where-to-start)
- [Run a complete analysis](#run-a-complete-analysis)
- [Use LDSC results you already have](#use-ldsc-results-you-already-have)
- [Choose your analysis settings](#choose-your-analysis-settings)
- [Find your results](#find-your-results)
- [Get help](#get-help)

## Understand the workflow

A full analysis uses two kinds of information from each trait:

1. **LDSC inputs** are used to estimate heritability, genetic relationships and
   intercepts. Raw summary statistics must first be filtered and converted into
   the format required by LDSC; this is called *munging*.
2. **GWAMA inputs** contain each SNP's alleles, frequency, sample size and
   association statistics. These are separate per-trait tables. The `prepare`
   command can create them from suitable GWAS VCFs.

LDSC produces estimates for the selected traits. PCA uses their genetic
correlation matrix by default to calculate PC1 loadings. GWAMA then uses those
loadings, the LDSC intercept matrix and the per-SNP tables to calculate PC1
association results.

**A munged file is an input to LDSC, not an LDSC result.** It also does not replace
the per-SNP table needed for GWAMA. Supplying completed LDSC results is enough
for QC and PCA; running GWAMA requires the per-SNP tables as well.

The [input guide](docs/INPUTS.md) shows the file formats, VCF fields, sample-size
rules and reference files. Check it before choosing a command: the two LDSC
backends accept different starting files.

## Install and activate

From a terminal with Git and internet access:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
conda activate ldsc-gpca
```

Activate `ldsc-gpca` in each new terminal. The installer sets up Python, R,
GenomicSEM, bcftools, pigz and a separate Python LDSC environment; the package
starts the LDSC environment automatically. If Conda is missing, follow the
installer's printed initialization instructions before activating.

The installer refuses to overwrite existing environments. See
[installation, updates and recovery](docs/REFERENCE.md#install-and-activate)
for existing installations, custom locations and platform limitations. Use the
[Docker guide](DOCKER.md) to run without host Conda.

## Choose where to start

### I want to run LDSC, PCA and GWAMA together

Use **`ldsc-gpca pipeline`**. With the default Python backend, start from suitable
GWAS VCFs, or reuse Python LDSC munged files. With
`--ldsc_backend genomicsem`, start from raw GWAS tables or compatible munged
files. A quantitative-trait VCF route is also available for GenomicSEM.

For either backend, also provide prepared GWAMA tables or suitable VCFs from
which to make them. The [complete pipeline guide](docs/PIPELINE.md) gives a
separate recipe for each starting situation and explains how to stop after PCA.

### I want to prepare files or run LDSC as a separate stage

- **`ldsc-gpca prepare`** reads GWAS VCFs and creates per-trait GWAMA tables.
  It does not run LDSC or PCA. See [preparation](docs/PREPARE.md).
- **`ldsc-gpca ldsc`** runs Python LDSC from supported VCFs, or from existing
  munged files when `--ldsc_only` is set. Its analysis output is
  `ldsc_results.csv`. See [Python LDSC](docs/PYTHON_LDSC.md).
- **`ldsc-gpca genomicsem ldsc`** munges raw GWAS tables and runs native LDSC.
  Select `--munged_dir` or `--munged_input` to reuse munged files instead.
  Its final analysis output is `genomicPCA_LDSC.RData`. See
  [GenomicSEM LDSC](docs/GENOMICSEM_LDSC.md).

Input modes are explicit. These commands do not decide that a file is already
munged from its filename or contents. Standalone `genomicsem ldsc` does not
read VCFs directly.

### I already have completed LDSC results

Use **`ldsc-gpca gpca`** for a complete Python LDSC results table, or
**`ldsc-gpca genomicsem gpca`** for the final native GenomicSEM RData object.
Neither command reruns LDSC.

Add `--validate_only` for QC and PCA reports alone. For SNP-level GWAMA, also
supply prepared per-trait tables with `--gpca_input_folder`, or provide VCFs in
the manifest for automatic preparation. Follow the
[Python GPCA/GWAMA guide](docs/PYTHON_GPCA.md) or
[GenomicSEM GPCA/GWAMA guide](docs/GENOMICSEM_GPCA.md).

## Run a complete analysis

This example uses **four quantitative traits and the Python backend**. Replace
all example paths with your files. Save `/data/traits_vcf.csv` as:

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

Each VCF must contain the fields needed for **both** Python LDSC extraction and
GWAMA preparation. You also need an ancestry-appropriate LD-score directory and
a HapMap allele table with `SNP,A1,A2` headers. See the
[VCF and reference requirements](docs/INPUTS.md#vcf-fields).

The bundled GWAMA output has no INFO column. Its selected-column export therefore
needs a constant you have justified for your downstream use. Set the shell
variable `GWAMA_INFO` to that value before this command; the guard below stops
if it is unset. This value is export metadata, **not an INFO filter or an
estimated quality score**.

```bash
ldsc-gpca pipeline \
  --input /data/traits_vcf.csv \
  --outdir /results/four_traits_python \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --splitby_chr nosplit \
  --gwama_output_info "${GWAMA_INFO:?Set a justified INFO export constant}"
```

Here `nosplit` uses one autosomal GWAMA table per trait. The default `split`
requires retained variants on every chromosome 1–22 for every trait. The output
directory must not already exist.

The pipeline prepares GWAMA tables, runs Python LDSC, checks the estimates,
calculates PC1, runs GWAMA and exports the results. To run LDSC and PCA without
GWAMA, replace the `--gwama_output_info` line with `--validate_only`; the pipeline
**still runs LDSC**. To start with raw tables, munged files, native GenomicSEM or
nondefault settings, use the [complete pipeline recipes](docs/PIPELINE.md).

## Use LDSC results you already have

Save `/data/selected_traits.csv` in the order you want to analyse:

```csv
traitname
Trait_A
Trait_B
Trait_C
Trait_D
```

Check a complete Python LDSC results table and calculate PCA:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_pca_check \
  --validate_only
```

For native GenomicSEM results:

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/native_ldsc/genomicPCA_LDSC.RData \
  --outdir /results/native_pca_check \
  --validate_only
```

Review retained traits, warnings, matrices and PC1 loadings. These commands do
not inspect the per-SNP GWAMA files. The workflow guides show the subsequent
[Python GWAMA command](docs/PYTHON_GPCA.md#run-pc1-gwama) and
[native GWAMA command](docs/GENOMICSEM_GPCA.md#run-pc1-gwama).

## Choose your analysis settings

Every workflow guide includes complete examples and grouped option tables with
defaults. Start with the setting you want to change:

- **INFO, MAF and variant filters:** [Python LDSC](docs/PYTHON_LDSC.md#change-variant-filters)
  or [GenomicSEM munging](docs/GENOMICSEM_LDSC.md#filters-sample-size-and-prevalence).
  LDSC filters do not also filter the GWAMA tables.
- **Sample size and binary traits:** [where N and prevalence come from](docs/INPUTS.md#sample-size-and-prevalence).
  `--gwama_output_n_eff` changes only the exported summary, not the analysis N.
- **High chi-square SNPs:** [Python cutoff and automatic mode](docs/PYTHON_LDSC.md#limit-large-chi-square-values)
  or [native automatic cutoff](docs/GENOMICSEM_LDSC.md#filters-sample-size-and-prevalence).
- **Original or trait-wide Python correlations:** [normalization choices](docs/PYTHON_GPCA.md#choose-the-correlation-and-pca-method).
  The default is original Python `rg`; `trait_wide` is optional.
- **Correlation or covariance PCA:** [Python](docs/PYTHON_GPCA.md#choose-the-correlation-and-pca-method)
  or [GenomicSEM](docs/GENOMICSEM_GPCA.md#choose-the-pca-method).
- **Variant IDs, P-value floor and chromosome layout:** [preparation](docs/PREPARE.md#choose-ids-layout-and-p-value-floor).
- **Failed estimates, missing traits and workers:** the option tables in your
  workflow guide explain the relevant stage and defaults. See also
  [worker and thread settings](docs/REFERENCE.md#cpu-workers-and-numerical-threads).
- **Output INFO, N_eff, filenames and compression:**
  [export settings](docs/INPUTS.md#export-metadata).

## Find your results

A pipeline run writes `prepare/` when preparation is needed, `ldsc/`, `gpca/`,
resolved manifests and `Pipeline_Run_Status.json`. Check the status file and QC
reports before interpreting results. Trait removal changes the trait set and
can change PC1.

In `gpca/`, look for the matrices and `GenomicPCA_PC1_Weights_Used.csv`. A full
successful run also writes:

- `<dataset_id>_GWAMA_combined_results.txt.gz`: the combined GWAMA result table.
- `harmonisation_input/<dataset_id>_GPCA_inputs.txt.gz`: the selected-column
  summary for downstream use, including INFO and N_eff export settings.
- `<dataset_id>_postprocess.json`: export audit information.

The selected-column summary has **12 columns**. It is a PC1 result, not one of
the nine-column per-trait inputs used to run GWAMA. The folder name does not
mean the file has already been harmonised with another dataset. See the
[output catalog and interpretation](docs/REFERENCE.md#find-and-interpret-the-outputs).

## Get help

```bash
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca genomicsem ldsc --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
```

Flags use underscores. Enable an on/off setting with the flag alone, such as
`--validate_only`. Omit it to keep its default. Both GPCA commands and pipeline
also support `--prepare_help` and `--postprocess_help`.

For errors, start with [troubleshooting](docs/REFERENCE.md#troubleshooting) and
the failed stage's log. Report reproducible problems through
[GitHub Issues](https://github.com/JIBINJOHNV/ldsc-gpca/issues), including the command,
error and software versions, without sharing private data.

Further details: [full reference](docs/REFERENCE.md), [Docker](DOCKER.md),
[Nextflow image and QC examples](examples/nextflow/README.md),
[raw upstream LDSC commands](docs/REFERENCE.md#raw-ldsc-commands), and
[methods and source attribution](docs/REFERENCE.md#sources-and-license).
