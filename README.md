# ldsc-gpca

Run LD score regression (LDSC), genomic principal component analysis (genomicPCA),
and SNP-level PC1 GWAMA from GWAS summary statistics. Use one pipeline command
for the complete analysis, or start from existing LDSC results.

The package supports **Python LDSC** and **GenomicSEM LDSC** and follows the
[Fürtjes genomicPCA procedure](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
It reports all principal components and uses **PC1 only** for GWAMA. The modified
GWAMA function is included.

This guide covers **v0.7.0**. Use `ldsc-gpca --version` to check your installation.

## Contents

- [Install and activate](#install-and-activate)
- [Choose a workflow](#choose-a-workflow)
- [Prepare your inputs](#prepare-your-inputs)
- [Run the complete pipeline](#run-the-entire-analysis-in-one-command)
- [Run individual stages](#run-individual-stages)
- [Run GPCA and GWAMA from existing LDSC results](#run-gpca-and-gwama)
- [Options and defaults](#options-and-defaults)
- [Find and interpret the outputs](#find-and-interpret-the-outputs)
- [Troubleshooting](#troubleshooting)
- [Further documentation](#further-documentation)

## Install and activate

From a terminal with Git and internet access:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
conda activate ldsc-gpca
ldsc-gpca --version
```

The installer sets up Python, R, GenomicSEM, bcftools, pigz and an isolated
Python LDSC environment. **Activate only `ldsc-gpca`**, including in each new
terminal; the package launches the LDSC environment automatically. If Conda is
missing, follow the installer's printed initialization instructions.

The installer refuses to overwrite existing environments. For a package-only
update, run the following from your updated clone:

```bash
conda activate ldsc-gpca
python -m pip install .
ldsc-gpca --version
```

This does not update R packages, bcftools or the child LDSC environment. Check
dependency changes when upgrading. The package requires Python 3.10–3.12;
the supplied main environment uses Python 3.11.

See [installation and recovery details](docs/REFERENCE.md#install-and-activate)
for custom locations and platform limitations, or [Docker installation](DOCKER.md)
to run without host Conda.

## Choose a workflow

| Your starting point | Use | Result |
| --- | --- | --- |
| GWAS VCFs, or munged files plus prepared GWAMA inputs | [`ldsc-gpca pipeline`](#run-the-entire-analysis-in-one-command) | LDSC → PCA → PC1 GWAMA → final exports |
| VCFs that need GWAMA input preparation | [`ldsc-gpca prepare`](#prepare-per-trait-gpca-inputs) | Per-trait GWAMA tables |
| VCFs or existing Python LDSC munged files | [`ldsc-gpca ldsc`](#run-python-ldsc) | Pairwise `ldsc_results.csv` |
| Raw GWAS tables or GenomicSEM-compatible munged files | [`ldsc-gpca genomicsem ldsc`](#run-genomicsem-ldsc) | `genomicPCA_LDSC.RData` |
| Complete Python LDSC results | [`ldsc-gpca gpca`](#run-gpca-and-gwama) | QC, PCA and optional PC1 GWAMA |
| Native GenomicSEM LDSC results | [`ldsc-gpca genomicsem gpca`](#run-gpca-and-gwama) | QC, PCA and optional PC1 GWAMA |

If LDSC is already complete, use a GPCA command to avoid rerunning it. Start with
`--validate_only` to review the matrix and PC1. Python LDSC and GenomicSEM can
produce different estimates; selecting a backend is an analysis choice.

## Prepare your inputs

You supply the GWAS data, ancestry-appropriate LD-score references and HapMap
reference. Check genome build, allele direction, sample-size conventions and,
for binary traits, prevalence. The package does not perform liftover or establish
allele alignment across your datasets.

- **Manifest:** comma-separated CSV with a header and one row per trait. Use
  unique, non-empty `traitname` values without whitespace or path separators.
  Names must match input filenames and LDSC results exactly. Manifest row order
  determines the analysis order; at least two traits must remain after QC.
- **LDSC input:** raw or munged summary statistics, depending on the selected
  command. A munged LDSC file is not a GWAMA input file.
- **GWAMA input:** one set of per-trait TSVs, prepared from VCFs or supplied by
  you. Required columns are `SNPID,CHR,BP,EA,OA,EAF,N,Z,P`. Extra columns and
  a different column order are accepted.
- **Paths:** replace the example `/data`, `/references` and `/results` paths.
  Prefer absolute paths and fresh output directories. Quote shell paths with
  spaces and CSV fields containing commas.

Full file schemas and sample-size rules are in the
[input reference](docs/REFERENCE.md#file-and-command-conventions).

## Run the entire analysis in one command

`ldsc-gpca pipeline` prepares inputs when needed, runs the selected LDSC backend,
then runs PCA, PC1 GWAMA and export. Its `--outdir` **must not already exist**.

### 1. Create a trait manifest

For quantitative GWAS VCFs, save `/data/pipeline_traits.csv`:

```csv
traitname,vcf_files,sample_prevalence,population_prevalence
protein1,/data/protein1.vcf.gz,NA,NA
protein2,/data/protein2.vcf.gz,NA,NA
protein3,/data/protein3.vcf.gz,NA,NA
```

Both prevalence headers are required. Relative file paths in a pipeline manifest
resolve beside that manifest. Each VCF must contain one GWAS sample and satisfy
the chosen route's requirements:

| Route | Required VCF fields and restrictions |
| --- | --- |
| GWAMA preparation, either backend | FORMAT `AF,ES,SE,LP,NEF`; preparation uses `EA=ALT`, `OA=REF`, `EAF=AF`, `Z=ES/SE`, `N=NEF`. |
| Python LDSC from VCFs | Also needs FORMAT `SI,EZ` and INFO `AF,EUR`. Binary traits using population prevalence additionally need FORMAT `NC,NCO`; follow the [Python sample-size rules](docs/REFERENCE.md#python-ldsc-input). |
| GenomicSEM LDSC from VCFs alone | Quantitative traits only, with both prevalences blank/NA. Generated LDSC tables use `N=NEF` and lack INFO: use previously QCed VCFs and verify what NEF represents. |

For GenomicSEM binary traits, supply raw tables through a `sumstats_file` manifest
column, or reuse properly munged files. See [GenomicSEM input modes](docs/REFERENCE.md#genomicsem-ldsc-inputs-choose-one-mode).
Do not substitute biologically different VCF fields just to satisfy the schema.

### 2. Select the backend and run

Use a whitespace-separated HapMap allele reference with `SNP,A1,A2` headers.
LD-score and weight directories must contain the required chromosome files;
see [reference-file requirements](docs/REFERENCE.md#all-python-ldsc-options).

The bundled GWAMA function does not output INFO. Full runs therefore need a
scientifically justified constant for the summary's INFO column. Set the shell
variable `GWAMA_INFO` to that value before using these commands. It is export
metadata, **not an estimate of imputation quality or a filtering threshold**;
there is no universal default. The shell guard prevents a run with an unset value.

Python LDSC:

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
  --input /data/pipeline_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/protein_python \
  --n_cores 4 \
  --dataset_id protein_pc1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

GenomicSEM LDSC, using the quantitative VCF manifest above:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/pipeline_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/protein_genomicsem \
  --n_cores 4 \
  --dataset_id protein_pc1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

The examples retain each backend's default chi-square filter; those defaults
differ. See [options and defaults](#options-and-defaults) before comparing backends.
Use a worker count appropriate to your CPU and memory allocation.

To run LDSC and inspect PCA first, add `--validate_only` and omit
`--gwama_output_info`. **Pipeline validation still runs LDSC**, but skips GWAMA
and export. Review the outputs, then use the matching [GPCA command](#run-gpca-and-gwama)
with the saved LDSC results to continue.

Already have munged and prepared files? Follow the
[pipeline reuse examples](docs/REFERENCE.md#starting-with-munged-and-prepared-files).
Python reuse requires both `--ldsc_only` and `--munged_dir`, including required
prevalence sidecars. Supply `--gpca_input_folder` for existing GWAMA tables.

### 3. Check completion and results

The pipeline writes `Pipeline_Run_Status.json`, resolved manifests in `manifests/`,
and stage outputs in `prepare/` when needed, `ldsc/` and `gpca/`.
Review the status, warnings, retained traits and PC1 reports before interpreting
the [GWAMA results](#find-and-interpret-the-outputs).

A failed stage stops the pipeline and leaves completed outputs available.
There is no automatic pipeline resume. Use individual commands to reuse those
outputs, or choose a new directory for another full run. Standalone Python LDSC
has a separate [verified-checkpoint restart option](docs/REFERENCE.md#restarting-interrupted-ldsc-runs).

## Run individual stages

Skip stages whose outputs you already have. Each stage has its own manifest
requirements, described below and in the linked reference.

### Prepare per-trait GPCA inputs

Use a CSV containing `traitname,vcf_files` (the quantitative pipeline manifest
above also works):

```bash
ldsc-gpca prepare \
  --input /data/pipeline_traits.csv \
  --outdir /results/prepared \
  --n_cores 4
```

Use `/results/prepared/gpca_inputs` as the later `--gpca_input_folder`.
Default `--splitby_chr split` produces files for chromosomes 1–22 and requires
retained variants on every chromosome for every trait. For whole-genome files,
use `--splitby_chr nosplit` during preparation **and** GPCA.
See [preparation options and QC](docs/REFERENCE.md#prepare-per-trait-gpca-inputs).

### Run Python LDSC

Save `/data/python_ldsc_traits.csv`. Unlike pipeline mode, this command requires
the `ref` column. Use `yes` for **every trait** to obtain complete GPCA coverage:

```csv
traitname,vcf_files,ref,population_prevalence,sample_prevalence
protein1,/data/protein1.vcf.gz,yes,NA,NA
protein2,/data/protein2.vcf.gz,yes,NA,NA
protein3,/data/protein3.vcf.gz,yes,NA,NA
```

```bash
ldsc-gpca ldsc \
  --input /data/python_ldsc_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/python_ldsc \
  --n_cores 4
```

The main result is `/results/python_ldsc/ldsc_results.csv`. GPCA needs every
unordered pair and self-pair: `k × (k + 1) / 2` unique rows for `k` traits.
Consistent symmetric duplicates are accepted. See
[munged-file reuse, filters and restart](docs/REFERENCE.md#run-python-ldsc).

### Run GenomicSEM LDSC

For raw GWAS tables, save `/data/genomicsem_traits.csv`:

```csv
traitname,sumstats_file,sample_prevalence,population_prevalence
protein1,/data/protein1.tsv,NA,NA
protein2,/data/protein2.tsv,NA,NA
protein3,/data/protein3.tsv,NA,NA
```

Each table needs `SNP,A1,A2,P`, a signed statistic (`BETA` or `Z`), and sample
size. See the [full schema and binary-trait rules](docs/REFERENCE.md#genomicsem-ldsc-inputs-choose-one-mode).

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/genomicsem_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/genomicsem_ldsc
```

Use the final `genomicPCA_LDSC.RData` and its retained-trait manifest
`Selected_Traits.csv` for GPCA. Do not use `genomicPCA_LDSC_raw.RData`.
By default, native LDSC drops traits with invalid heritability; review its audit
or choose `--invalid_h2_action error` to stop instead.

## Run GPCA and GWAMA

These commands read existing LDSC results and do not rerun LDSC.
For Python results, create `/data/selected_traits.csv`:

```csv
traitname
protein1
protein2
protein3
```

### Check the matrix and PC1 first

Python LDSC results:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_gpca_qc \
  --validate_only
```

GenomicSEM results:

```bash
ldsc-gpca genomicsem gpca \
  --input /results/genomicsem_ldsc/Selected_Traits.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicPCA_LDSC.RData \
  --outdir /results/genomicsem_gpca_qc \
  --validate_only
```

Review PC1 weights, eigenvalues, CTI and any warning/removal reports.
`--validate_only` does **not** inspect the per-variant GWAMA files.

### Run SNP-level PC1 GWAMA

Supply prepared per-trait files and a justified INFO export value as described
in the [pipeline instructions](#2-select-the-backend-and-run):

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/python_gpca_gwama \
  --n_cores 4 \
  --dataset_id protein_pc1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

For native results, use `ldsc-gpca genomicsem gpca`, the native RData and the
matching retained-trait manifest. Keep other options appropriate to your files.
Both commands can also [prepare inputs automatically from VCFs](docs/REFERENCE.md#automatic-vcf-preparation).

| Input layout | Required filename for each retained trait |
| --- | --- |
| `split` (default) | `{traitname}_chr1_GenomicPCA_inputs.tsv` through `{traitname}_chr22_GenomicPCA_inputs.tsv` |
| `nosplit` | `{traitname}_GenomicPCA_inputs.tsv` |

Existing tables must use compatible alleles, positions and sample-size conventions.
Schema validation alone does not establish their scientific suitability.

## Options and defaults

Managed options use underscores, for example `--n_cores`. Enable a flag with
the flag alone: `--validate_only`, without `true` or `false`. Command help lists
all choices and defaults; use the backend-specific help for pipeline mode:

```bash
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca genomicsem ldsc --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
```

| Setting | Default | When to change it |
| --- | --- | --- |
| Pipeline `--ldsc_backend` | `python` | Use `genomicsem` for native GenomicSEM LDSC and its GPCA reader. |
| `--splitby_chr` | `split` | Use `nosplit` for one whole-genome GWAMA file per trait; preparation and GPCA must agree. |
| Python LDSC `--chisq_max` | Extra filter disabled | Set a positive integer, or `auto` for `max(80, 0.001 × maximum eligible N)` per trait. |
| GenomicSEM LDSC `--chisq_max` | Native automatic cutoff | Set a positive finite number to override it; omit for automatic mode. The string `auto` is Python-only. |
| Python GPCA `--rg_normalization` | `pair` | Use `trait_wide` to select the exported `rg_trait_wide` values instead of original `rg`. |
| GPCA `--pca_matrix` | `correlation` | Choose `covariance` only with an appropriate heritability scale and interpretation. |
| GPCA `--pc1_orientation` | `tutorial` | Use `as_computed` to keep the eigendecomposition's arbitrary sign. Tutorial mode flips the whole vector when its median is negative. |
| GPCA `--failed_ldsc_action` | `error` | `drop_traits` permits audited removal; removing traits changes the component. |
| Export `--gzip_level` | `1` | Choose 1–9 to trade compression time for file size. Numeric content is unchanged. |

**Chi-square filtering** removes variants with `Z²` above the cutoff from LDSC
inputs only. It does not filter GWAMA variants or apply a PCA significance cutoff.
Matching cutoffs does not guarantee identical backend estimates.

**Trait-wide normalization** divides pairwise genetic covariance by the square
root of the product of the two self-pair heritabilities. Original `rg` and CTI
remain unchanged in the results; selecting `trait_wide` changes the PCA input
and may change PC1/GWAMA results.
Original SE/Z/P still describe original `rg`. This option does not make Python
regression weighting identical to GenomicSEM. See the
[formula and requirements](docs/REFERENCE.md#optional-trait-wide-correlation-normalization).

**Workers:** standalone defaults are preparation `4`, Python LDSC `5`,
GenomicSEM LDSC `1`, and GPCA `0` (automatic). Pipeline uses the selected LDSC
default for LDSC and GPCA/export, with separate `--prepare_workers 4`.
Pipeline requires positive `--n_cores`. Native LDSC regression is sequential;
split GWAMA has at most 22 chromosome jobs. More workers need more memory and
can increase storage contention. See [worker details](docs/REFERENCE.md#cpu-workers-and-numerical-threads).

All remaining options, including failure policies, input reuse and heritability
scale selection, are in the [complete option reference](docs/REFERENCE.md#options-defaults-and-help).

## Find and interpret the outputs

GPCA reports are under the GPCA command's `--outdir`, or `<pipeline_outdir>/gpca/`.
Start with these files:

| File | What to check |
| --- | --- |
| `GenomicPCA_PC1_Weights_Used.csv` | Signed PC1 loadings used for GWAMA, in trait order. |
| `GenomicPCA_PC1_QC.csv` | PC1 validity and direction; this is not overall run status. |
| `GenomicPCA_All_PCs_Variance_Explained.csv` | Eigenvalues, explained percentages and negative-eigenvalue flags. |
| `GenomicPCA_PC1_Protein_Contributions.csv` | Trait contributions; these are not causal contributions. |
| `GenomicPCA_PCA_Matrix_Used.csv` | Exact genetic matrix decomposed. |
| `GenomicPCA_CTI_Used.csv` | LDSC intercept/error-covariance matrix used by GWAMA. |
| `GWAMA_Run_Status.csv` | GWAMA success/failure and current-run output paths. |
| `{dataset_id}_GWAMA_combined_results.txt.gz` | Combined GWAMA results, sorted by chromosome and position. |
| `harmonisation_input/{dataset_id}_GPCA_inputs.txt.gz` | Selected-column downstream summary. |
| `{dataset_id}_postprocess.json` | Export sources, row count and metadata overrides. |

The downstream summary has **12 tab-separated columns**:

```text
SNPID  CHR  BP  EA  OA  EAF  N_eff  BETA  SE  Z  PVAL  INFO
```

It is not interchangeable with the nine-column per-trait GWAMA inputs. Export
preserves GWAMA-reported statistics; it does not harmonise alleles or establish
effect-size/sample-size calibration. `--gwama_output_info` and the optional
`--gwama_output_n_eff` replace only those summary columns, without changing GWAMA
weights, Z, P or the combined table.

For interpretation, keep these distinctions in mind:

- Correlation PCA uses genetic correlations; GWAMA's CTI uses LDSC intercepts.
  Pairwise P values and SEs are diagnostics, not PCA filters or weights.
- Negative eigenvalues and correlations outside `[-1,1]` are reported. The
  package does not silently repair the matrix; CTI must be positive definite.
  Negative eigenvalues also complicate interpretation of explained percentages.
- The Python table lacks GenomicSEM's full sampling covariance matrices `V`
  and `V_Stand`. Procedures requiring them cannot be run from that table;
  neither GPCA command runs `paLDSC`.
- Review retained/dropped traits, warnings and command exit status. A PC1 PASS
  or an existing output file does not establish that the whole run succeeded.

Keep your commands, version, manifests, input/reference provenance and audit
files. See [all outputs](docs/REFERENCE.md#find-and-interpret-the-outputs) and
[QC and scale interpretation](docs/REFERENCE.md#understand-qc-and-scale-choices).

## Troubleshooting

| Problem | What to check |
| --- | --- |
| Command missing or new options unavailable | Activate `ldsc-gpca`, check `--version`, then reinstall the package from the updated clone. |
| Manifest rejected | Check comma separation, exact headers, trait names and the selected input mode. |
| Missing VCF fields | Preparation and Python LDSC require different fields; check the route's schema. |
| Missing chromosome inputs | Use matching `split`/`nosplit` settings. Split mode needs chromosomes 1–22 for every retained trait. |
| Missing LDSC pairs | Use `ref=yes` for every trait in standalone Python LDSC. One raw multi-trait `--rg` command does not generate all combinations. |
| Invalid heritability, ambiguous scale or non-positive-definite CTI | Inspect the LDSC/QC reports and input conventions before changing QC policy. |
| GWAMA finished but export failed | Check INFO/N_eff availability, required output fields and existing export paths. Preserve completed results; do not invent constants to bypass the error. |
| Interrupted pipeline | Inspect `Pipeline_Run_Status.json`. Reuse completed stage outputs with individual commands, or start in a new output directory. |

See [detailed troubleshooting](docs/REFERENCE.md#troubleshooting) for installation,
runtime and restart issues.

## Further documentation

- [Workflow and option reference](docs/REFERENCE.md): complete schemas, defaults,
  QC rules, output definitions and advanced commands.
- [Docker](DOCKER.md) and [Nextflow / Google Batch](examples/nextflow/README.md).
- [Pipeline validation](tests/PIPELINE_VALIDATION.md): tests performed and their limits.
- [Method and upstream sources](docs/REFERENCE.md#sources-and-license), including
  [bundled GWAMA provenance](src/ldsc_gpca/r/vendor/README.md).

The repository has no project-level LICENSE file. Upstream components retain
their own licenses and attribution.
