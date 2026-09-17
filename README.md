# ldsc-gpca

Run genomic principal component analysis (genomicPCA) and SNP-level PC1 GWAMA
from GWAS summary statistics, through one command-line interface.

Use either **Python LDSC results** or **GenomicSEM LDSC RData**. The package can
also run LDSC and prepare the per-trait summary-statistic files needed by GWAMA.
It calculates all PCs and reports their variability; **only PC1 is used for GWAMA**.
The method follows the [Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).

## Start here

1. [Choose your starting point](#choose-a-workflow).
2. Install using [Conda](#install-and-activate) or [Docker](#docker-install-and-run-without-host-conda).
3. [Check LDSC and inspect PCA](#run-gpca-and-gwama), then run GWAMA if needed.
4. [Review the outputs](#find-and-interpret-the-outputs) and retain the QC reports.

For detailed options and scientific caveats, use the [workflow reference](docs/REFERENCE.md).
For pipelines, see [Nextflow and Google Batch](examples/nextflow/README.md).

## Choose a workflow

| What you already have | Where to start |
| --- | --- |
| Complete Python LDSC results | [Python LDSC → GPCA](#python-ldsc-qcpca-example) |
| GenomicSEM LDSC RData | [GenomicSEM → GPCA](#genomicsem-qcpca-example) |
| GWAS files, but no LDSC results | [Choose an LDSC backend](#generate-ldsc-results) |
| LDSC results, but no per-trait GWAMA inputs | [Prepare GPCA/GWAMA inputs](#prepare-per-trait-gpca-inputs) |

Both GPCA backends need a **trait manifest and LDSC results**. To run SNP-level
GWAMA, they also need **per-trait GWAMA input files** in `--gpca_input_folder`.
An LDSC matrix does not contain those SNP-level statistics.

With `--validate_only`, you can inspect LDSC QC and PCA without GWAMA files.
That mode does **not** read or validate per-variant GWAMA inputs.

Before analysis, check ancestry, genome build, allele alignment and sample-size
conventions across your traits and references. File conversion and SNP selection
are not allele harmonisation. Use fresh output directories and consistent,
non-empty trait names; avoid whitespace and numeric-only identifiers.

## Install and activate

Choose this route for native execution. The installer supplies Python, R,
GenomicSEM, bcftools and an isolated CBIIT LDSC environment; native commands do
not use Docker. Git and internet access are required.

**Install once**, with Conda initialized in your terminal:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
```

**Each time you open a terminal:**

```bash
conda activate ldsc-gpca
ldsc-gpca --help
```

Activate only `ldsc-gpca`. The package launches the separate `ldsc-gpca-ldsc`
environment internally; you do not switch between environments.

The installer will not overwrite existing environments. If Conda is absent, it
can install Miniforge and prints the shell-initialization steps to follow.
It detects compatible Mamba arguments automatically. To use Conda explicitly,
add `--manager conda` to the installation command.
See [installation troubleshooting and updates](docs/REFERENCE.md#install-and-activate).

Linux amd64 Docker installation has been tested. Native Linux/macOS platforms
are detected, but the pinned dependencies have not been verified on every platform.
Installing with `pip` alone does not install the R, bcftools or LDSC dependencies.

## Docker: install and run without host Conda

Install and start [Docker Desktop](https://docs.docker.com/desktop/) or
[Docker Engine](https://docs.docker.com/engine/install/), then build from the repository:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
docker build --platform linux/amd64 -t ldsc-gpca:0.5.0 .
docker run --rm --platform linux/amd64 ldsc-gpca:0.5.0 ldsc-gpca --help
```

The image is **not published to a registry**; build it locally. Apple Silicon uses
amd64 emulation. No host Conda/R installation or activation is needed.

To run the Python-LDSC QC/PCA example with files from your computer:

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /absolute/path/to/data:/work -w /work \
  ldsc-gpca:0.5.0 ldsc-gpca gpca \
  --input /work/selected_traits.csv \
  --python_ldsc /work/ldsc_results.csv \
  --outdir /work/qc_results --validate_only
```

Replace the host directory with yours. Results appear in its `qc_results/` folder.
Manifest paths must also be accessible inside the container.

The image uses `USER root` and `ENTRYPOINT []`. Always include the executable
`ldsc-gpca` after the image name; the local UID/GID override above avoids root-owned
outputs on Linux. Nextflow can supply its shell command directly.

[Docker guide](DOCKER.md) · [Nextflow and Google Batch setup](examples/nextflow/README.md)

Google Batch additionally needs a registry-hosted image, GCS storage and appropriate
IAM permissions. Local Nextflow tests passed; live GCP execution remains untested.

## Run GPCA and GWAMA

The commands below assume native activation. Docker users can run the same
commands after the image name, with the necessary input/output folders mounted.
All example paths are placeholders; replace them with your files.

### Common inputs

Create `selected_traits.csv` (**comma-separated**, with a header):

```csv
traitname
protein1
protein2
```

Select at least two unique traits. Names must match your LDSC results and GWAMA
filenames exactly. Manifest row order controls the analysis order; extra traits
in the LDSC file are allowed.

| Backend | LDSC input |
| --- | --- |
| Python | Complete pairwise table, including self-pairs; CSV, TSV or whitespace-delimited, `.gz` accepted. [Required columns](docs/REFERENCE.md#python-ldsc-qcpca-example) |
| GenomicSEM | Binary RData containing `LDSCoutput` with `S,V,I,S_Stand,V_Stand`. Use the final `genomicPCA_LDSC.RData`, not the preliminary raw object. |

Managed Python LDSC runs export each batch directly from its native in-memory
result table to `<batch>.results.csv` using `%.17g` numerical formatting.
The compiler reads each CSV once; readable logs are retained for troubleshooting
but never parsed in the normal workflow. No extra regression jobs or estimate
recomputation are performed. All correlation, heritability and intercept fields
retain their available floating-point precision. Missing or invalid numerical
exports stop compilation instead of falling back to rounded log tables.

For old runs only, explicitly passing `.log` files to `compile_results()`
retains detailed-log recovery. This preserves printed precision, not every digit
of the original estimates. No SE is imputed. Existing PCA/GWAMA outputs are not
automatically regenerated. Raw `ldsc-gpca ldsc.py` remains unchanged and does
not use the managed numerical exporter.

For GWAMA, both backends additionally require **tab-separated** per-trait files
with these nine columns, in this order:

```text
SNPID  CHR  BP  EA  OA  EAF  N  Z  P
```

The spaces above separate column names for readability; the files must use tabs.
Aliases `A1,A2,p` are accepted for `EA,OA,P`.

- Default split layout: `{traitname}_chr{CHR}_GenomicPCA_inputs.tsv`, for every chromosome 1–22.
- Whole-genome layout: `{traitname}_GenomicPCA_inputs.tsv`, with `--splitby_chr nosplit`.

The GWAMA reader checks the schema, not complete row-level numerical QC. Use
already validated, consistently aligned files. If needed, [prepare them from supported VCFs](#prepare-per-trait-gpca-inputs).

### Start with QC/PCA only

Run **one** of the following, according to your LDSC backend.

#### Python LDSC: QC/PCA example

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --python_ldsc /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_gpca_qc \
  --validate_only
```

#### GenomicSEM: QC/PCA example

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_path /results/genomicsem_ldsc/genomicPCA_LDSC.RData \
  --outdir /results/genomicsem_gpca_qc \
  --validate_only
```

Inspect the [PC1, variance-explained and trait-QC reports](#find-and-interpret-the-outputs)
before proceeding. A PC1 PASS is not a guarantee that every QC check passed.

### Run SNP-level PC1 GWAMA after reviewing QC

The modified GWAMA v1.2.6 R function is bundled; no `--source_path` is needed.
Use a new output directory and omit `--validate_only`.

**Resolve export metadata first:** the bundled GWAMA function does not output INFO.
Automatic export therefore needs a scientifically justified `--gwama-output-info`
value; otherwise the command fails at export after GWAMA completes. This value is
user-supplied metadata, not an estimated INFO score. Do not invent one to make a run pass.

The examples require you to set `GWAMA_INFO` to your justified value in `[0,1]`.
The shell guard deliberately stops if it is unset. See the [export policy](docs/REFERENCE.md#automatic-gwama-export).

**From Python LDSC:**

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --python_ldsc /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/python_gpca_gwama \
  --dataset-id cluster1 \
  --gwama-output-info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

**From GenomicSEM LDSC:**

```bash
ldsc-gpca genomicsem gpca \
  --input /data/selected_traits.csv \
  --ldsc_path /results/genomicsem_ldsc/genomicPCA_LDSC.RData \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/genomicsem_gpca_gwama \
  --dataset-id cluster1 \
  --gwama-output-info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

Both default to correlation PCA, tutorial PC1 orientation and chromosome-split
GWAMA. Use `--splitby_chr nosplit` for whole-genome input files.
Automatic export follows GWAMA; it is not a separate command.

### Prepare per-trait GPCA inputs

Skip this step if you already have the validated nine-column files above.
Otherwise, create `prepare_traits.csv` (**CSV**):

```csv
traitname,vcf_files
protein1,/data/protein1.vcf.gz
protein2,/data/protein2.vcf.gz
```

Each VCF must have one sample and FORMAT fields `AF,ES,SE,LP,NEF`.
`ES` and `AF` must refer to ALT; `LP` is `-log10(P)`.

```bash
ldsc-gpca prepare \
  --input /data/prepare_traits.csv \
  --outdir /results/prepared \
  --splitby_chr split
```

This creates `gpca_inputs/` plus per-trait QC summaries and original-record issue
reports. It uses `EA=ALT`, `OA=REF`, `N=NEF` and `Z=ES/SE`; verify your N convention.
It does not harmonise alleles, lift over coordinates or apply MHC/INFO-score filters.
Split mode requires retained variants on all 22 autosomes for every trait.

Optional `--write-munge-inputs --hapmap-file FILE` also writes **unmunged**
HapMap-selected tables. This is SNP selection, not allele alignment or LDSC.
See [preparation formats, QC and options](docs/REFERENCE.md#prepare-per-trait-gpca-inputs).

GPCA can also prepare VCFs automatically when `--gpca_input_folder` is omitted
and its manifest includes `vcf_files`. Use `--prepare-help` for those options.

## Generate LDSC results

Skip this section if you already have valid LDSC results. Choose **one backend**;
neither LDSC command automatically runs GPCA or creates the nine-column GWAMA inputs.

### Run Python LDSC

Use `ldsc-gpca ldsc` for the package's VCF extraction → CBIIT LDSC workflow.
Its fixed schema requires INFO `AF,EUR` and FORMAT `SI,AF,EZ,LP,NEF`;
binary traits with population prevalence also need `NC,NCO`.
This is not a generic ancestry-independent GWAS-VCF reader.

The CSV manifest requires `gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence`.
Set `ref=yes` for every selected trait to generate complete GPCA pair/self-pair coverage.

Use `--chisq-max 80` to reproduce GenomicSEM-style extreme-statistic filtering:
each munged trait is filtered independently to `Z^2 <= 80` immediately before
pairwise LDSC. The original VCF and munged files are preserved. This option is
also supported with `--ldsc_only`; it is intentionally not forwarded to native
Python LDSC's different cross-product implementation.

Advanced users can invoke the pinned upstream console scripts directly:

```bash
ldsc-gpca ldsc.py --help
ldsc-gpca munge_sumstats.py --help
```

Everything following the script name is passed through unchanged. These raw
commands bypass manifest checks, managed per-trait filtering, retries, provenance
and result compilation; raw `ldsc.py --rg --chisq-max` retains native LDSC semantics.

[Python LDSC command, manifest example, references and prevalence rules →](docs/REFERENCE.md#run-python-ldsc)

### Run GenomicSEM LDSC

Use `ldsc-gpca genomicsem ldsc` for GWAS tables → native munging → LDSC.
The default CSV manifest requires `traitname,munge_inputs,sampleprevalence,populationprevalence`.
For quantitative traits use blank/NA prevalence entries; binary liability conversion
requires both valid sample and population prevalence.

Already munged? Use `--munge-output DIRECTORY` or `--munged-input` to skip munging.
The final outputs include `genomicPCA_LDSC.RData` and the retained `Selected_Traits.csv`.

[GenomicSEM commands, table columns, references and existing-munged-file options →](docs/REFERENCE.md#run-genomicsem-ldsc)

## Find and interpret the outputs

Start with these CSV files in your GPCA `--outdir` (written once the run reaches PCA):

| File | Question it answers |
| --- | --- |
| `GenomicPCA_All_PCs_Variance_Explained.csv` | How much estimated genetic variability does each PC describe? |
| `GenomicPCA_PC1_QC.csv` | Is PC1 valid, and which direction convention was used? |
| `GenomicPCA_PC1_Protein_Contributions.csv` | What are each retained trait's signed loading and PC1 contribution? |
| `GenomicPCA_PC1_Weights_Used.csv` | Which loadings were supplied to GWAMA? |
| `GWAMA_Run_Status.csv` | Which GWAMA jobs succeeded, or were intentionally skipped in QC/PCA-only mode? |

PC1 contribution is the squared eigenvector coefficient, expressed as a percentage;
it is not a causal contribution. Negative eigenvalues make raw variance percentages
non-standard. See [definitions and audit files](docs/REFERENCE.md#find-and-interpret-the-outputs).

After successful GWAMA **and export**, `--dataset-id cluster1` produces:

```text
outdir/
  cluster1_GWAMA_combined_results.txt.gz
  cluster1_postprocess.json
  harmonisation_input/
    cluster1_GPCA_inputs.txt.gz
```

The compressed tables are tab-separated. The file in `harmonisation_input/` is a
12-column downstream summary, **not** the nine-column input to this package's GWAMA.
Export does not harmonise alleles or comprehensively validate every numerical field.
[Export columns, INFO/N_eff overrides and limitations →](docs/REFERENCE.md#automatic-gwama-export)

## QC and scientific interpretation

- **Defaults stop on missing/invalid selected LDSC estimates.** Optional trait removal changes the analysis; inspect its audit and assess sensitivity.
- **Correlation PCA is the default.** Covariance PCA is available but depends on phenotype/heritability scales; do not silently mix or convert scales.
- **PC1's sign is arbitrary.** The tutorial default reverses the entire loading vector when its median is negative. Interpret effect signs relative to that recorded orientation.
- **No automatic matrix repair.** Out-of-range rg and negative PCA eigenvalues warn by default; a non-positive-definite CTI stops the run.
- **Python pairwise tables do not provide native `V`/`V_Stand`.** The package does not fabricate them or run `paLDSC`.

[QC defaults and scale policies](docs/REFERENCE.md#understand-qc-and-scale-choices) ·
[Troubleshooting](docs/REFERENCE.md#troubleshooting-and-limitations)

## Help and reproducibility

Use the help for your selected command, for example:

```bash
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca genomicsem ldsc --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
```

Help lists required columns, separators, defaults and choices. For either GPCA
backend, `--prepare-help` and `--postprocess-help` show the additional option groups.
Copy flag names exactly: the CLI uses both underscores and hyphens.

Keep the repository commit, image digest/environment records, commands, manifests,
reference checksums and QC/removal reports. Environment YAMLs are specifications,
not complete lockfiles. [Reproducibility and package layout](docs/REFERENCE.md#reproducibility-and-implementation)

Automated and local real-tool tests have passed, including Docker/Nextflow checks.
A full real-data LDSC-to-GWAMA run, live GCP execution and all native platforms have
not been verified. Software tests do not establish scientific suitability for publication.
[Build validation details](DOCKER.md#empty-entrypoint-build-verification)

## Sources and license

[genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html) ·
[CBIIT LDSC](https://github.com/CBIIT/ldsc) ·
[GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) ·
[bcftools](https://samtools.github.io/bcftools/bcftools.html)

Pinned source revisions and attribution are listed in the [reference guide](docs/REFERENCE.md#sources-and-attribution).
No redistribution license has been selected for this repository; third-party
software retains its own licenses.
