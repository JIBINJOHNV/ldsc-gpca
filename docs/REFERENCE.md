# Workflow and option reference

[Back to the README](../README.md) · [Docker](../DOCKER.md) · [Nextflow / Google Batch](../examples/nextflow/README.md)

Detailed input contracts, complete option tables, defaults, restart behavior,
QC policies and output definitions. Start with the [README](../README.md) to
choose a workflow, then use its command guide. The [input guide](INPUTS.md)
explains the file types and sample-size rules used across commands.

For complete commands demonstrating optional settings, use the workflow guides:
[pipeline](PIPELINE.md), [preparation](PREPARE.md), [Python LDSC](PYTHON_LDSC.md),
[GenomicSEM LDSC](GENOMICSEM_LDSC.md), [Python GPCA/GWAMA](PYTHON_GPCA.md), or
[GenomicSEM GPCA/GWAMA](GENOMICSEM_GPCA.md).

## Contents

- [Options, defaults and help](#options-defaults-and-help)
- [Install and activate](#install-and-activate)
- [File and command conventions](#file-and-command-conventions)
- [Run the entire analysis in one command](#run-the-entire-analysis-in-one-command)
- [Prepare per-trait GPCA inputs](#prepare-per-trait-gpca-inputs)
- [Run Python LDSC](#run-python-ldsc)
- [Run GenomicSEM LDSC](#run-genomicsem-ldsc)
- [Run GPCA and GWAMA](#run-gpca-and-gwama)
- [INFO filtering, sample size and export metadata](#info-filtering-sample-size-and-export-metadata)
- [Understand QC and scale choices](#understand-qc-and-scale-choices)
- [Find and interpret the outputs](#find-and-interpret-the-outputs)
- [CPU workers and numerical threads](#cpu-workers-and-numerical-threads)
- [Raw LDSC commands](#raw-ldsc-commands)
- [Docker: install and run without host Conda](#docker-install-and-run-without-host-conda)
- [Troubleshooting](#troubleshooting)
- [Help and reproducibility](#help-and-reproducibility)
- [Sources and license](#sources-and-license)

## Options, defaults and help

Every managed command groups its help by purpose: required inputs, input modes,
filters or PCA settings, diagnostics, and execution options as applicable.
Options are separated by a blank line. Each setting shows its default or explains
what happens when omitted; fixed choices are listed beside that option.
Required inputs have no default. Help/version actions simply display information.

```bash
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca ldsc --ldsc_backend genomicsem --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
```

Both GPCA commands show **analysis, automatic VCF preparation and GWAMA export**
options in `--help`. Use `--prepare_help` or `--postprocess_help` to display only
that section. These help commands work without starting an analysis or R.
The raw `ldsc.py` and `munge_sumstats.py` commands retain upstream help formatting.

`CHOICE`, `INTEGER`, `FLOAT`, `DIRECTORY` and example filenames in help are
placeholders; replace them with a listed choice, number or path. For example,
write `--rg_normalization trait_wide`, not `--rg_normalization CHOICE`.
An **off** flag is omitted; enable it with the flag alone, such as `--validate_only`.
**Unset** means omitted, and the accompanying explanation gives its behavior.
Environment-dependent defaults identify the environment variable or fallback.

Frequently used choices and defaults:

| Command | Option | Default when omitted | Accepted choices or values |
| --- | --- | --- | --- |
| `ldsc` and `pipeline` | `--ldsc_backend` | `python` | `python`, `genomicsem`; selects matching munging and regression |
| `prepare --mode ldsc`/`both` | `--munge_backend` | `python` | `python`, `genomicsem`; standalone munging only, omit with `--raw_only` |
| `prepare` and both GPCA commands | `--splitby_chr` | `split` | `split`, `nosplit` |
| Python `ldsc` | `--chisq_max` | Disabled | Positive integer, or `auto` for the per-trait automatic cutoff |
| `genomicsem ldsc` | `--chisq_max` | GenomicSEM automatic cutoff | Positive finite number; omit to use automatic mode |
| Python `ldsc` | `--result_failure_action` | `error` | `error`, `report`, `drop_traits` |
| `genomicsem ldsc` | `--invalid_h2_action` | `drop` | `drop`, `error` |
| Both GPCA commands | `--failed_ldsc_action` | `error` | `error`, `drop_traits` |
| Python-table `gpca` | `--rg_normalization` | `pair` (original Python `rg`) | `pair`, `trait_wide` (precomputed `rg_trait_wide`) |
| Python-table `gpca` | `--heritability_scale` | `auto` | `auto`, `observed`, `liability`, `mixed` |
| Both GPCA commands | `--pca_matrix` | `correlation` | `correlation`, `covariance` |
| Both GPCA commands | `--pc1_orientation` | `tutorial` | `tutorial`, `as_computed` |
| Both GPCA commands | `--rg_out_of_range_action`, `--negative_eigen_action` | `warn` | `warn`, `error` |
| `prepare` and GPCA automatic preparation | `--gpca_id_source` | `chr_pos_ref_alt` | `chr_pos_ref_alt`, `vcf_id` |
| `prepare` and GPCA automatic preparation | `--munge_id_source` | `vcf_id` | `chr_pos_ref_alt`, `vcf_id` |
| Both GPCA commands, final export | `--gzip_level` | `1` | Integers 1–9; affects compression, not numeric results |

Worker defaults depend on the command: `prepare --n_cores 4`, Python
`ldsc --n_cores 5`, `genomicsem ldsc --n_cores 1`, and GPCA `--n_cores 0`
(automatic). GPCA's optional VCF preparation separately uses `--prepare_workers 4`.
Pipeline uses the selected LDSC backend's worker default for LDSC and GPCA/export;
its `--prepare_workers` is separate and its `--n_cores` must be positive.
See [worker behavior](#cpu-workers-and-numerical-threads) before increasing them.

For all options, constraints and file requirements, use the tables for
[preparation](#all-preparation-options), [Python LDSC](#all-python-ldsc-options),
[GenomicSEM LDSC](#all-genomicsem-ldsc-options),
[shared GPCA](#all-shared-gpca-analysis-options),
[Python-only GPCA](#additional-options-for-python-table-gpca-only),
[automatic preparation](#automatic-vcf-preparation), and
[automatic export](#all-automatic-export-options).

## Install and activate

From a terminal with Git and internet access:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
conda activate ldsc-gpca
ldsc-gpca --version
ldsc-gpca --help
```

The installer creates the main `ldsc-gpca` environment and an isolated
`ldsc-gpca-ldsc` environment. It installs/checks Python, R, GenomicSEM, bcftools,
pigz and CBIIT LDSC. **Activate only `ldsc-gpca`**; the package launches the child
LDSC environment automatically. Run `conda activate ldsc-gpca` in each new terminal.

If Conda is absent, the installer can bootstrap checksum-verified Miniforge.
Follow its printed shell initialization command before activation. It does not
use sudo or edit your shell startup files. Linux x86_64/aarch64 and macOS
Intel/Apple Silicon are detected; pinned dependency availability can limit native
installation. A fresh install has not been tested on every supported architecture.

| Installer argument | Default | Meaning |
| --- | --- | --- |
| `--no-activate` | Absent | Install/check only. Recommended for the commands above; activate explicitly afterward. Without it, terminal execution opens an activated Bash shell; sourcing can activate the current Bash/zsh shell. Noninteractive execution skips activation. |
| `--manager auto` | `auto` | Environment creator: `auto`, `conda`, or `mamba`. Auto uses compatible Mamba if available, otherwise Conda. |
| Positional `ENV_ROOT` | Unset | Advanced custom-prefix installation under `ENV_ROOT/main` and `ENV_ROOT/ldsc`; activate with `source ENV_ROOT/activate.sh`. Default installation uses the named environments above. |
| `--help` | Off | Show installer usage without installing. |

The installer refuses to overwrite existing environments. For a package-only
update with unchanged dependencies, run these commands from your updated clone:

```bash
conda activate ldsc-gpca
python -m pip install .
ldsc-gpca --version
```

For an existing environment that lacks the faster gzip compressor:

```bash
conda activate ldsc-gpca
conda install -c conda-forge pigz=2.8
```

Review changes to `environment.yml` and `environment.ldsc.yml` when upgrading;
`pip install .` does not install/update R packages, bcftools or the child LDSC
environment. The Python package requires Python 3.10–3.12; the provided main
environment uses 3.11. Activation and recovery details follow below.

You supply the GWAS data, ancestry-appropriate LD references and HapMap lists.
Check genome build, allele direction, sample-size and prevalence conventions
across inputs; installation cannot establish their scientific compatibility.

### Activation, custom locations and recovery

<details>
<summary>No Conda yet, or “conda activate” does not work?</summary>

The installer can download checksum-verified
[Miniforge 26.7.2-0](https://github.com/conda-forge/miniforge/releases/tag/26.7.2-0)
when Conda is absent. This requires `curl` and `sha256sum` or macOS `shasum`.
It does not use sudo or edit your shell startup files.

After installation, follow the exact Conda initialization command printed by
the installer. For persistent setup, run its printed `conda init bash` command
(or `conda init zsh` for zsh), then reopen your terminal. After that, use:

```bash
conda activate ldsc-gpca
```

If you already have Conda but activation is unavailable, initialize that Conda
installation for your shell first. Do not mix multiple Conda installations.

</details>

<details>
<summary>Existing installations, custom locations and installation records</summary>

The default `--manager auto` checks Mamba's `env create` interface and falls back
to Conda if that interface is unavailable. Older interfaces without `--yes`
receive noninteractive confirmation through command-local environment variables;
the installer does not change your Conda configuration or upgrade your base environment.
You can explicitly select `--manager conda` or `--manager mamba`.

If an older installer stopped with `unrecognized arguments: --yes`, update the
repository and retry with Conda:

```bash
git pull --ff-only
bash scripts/setup_environments.sh --no-activate --manager conda
```

Check `conda env list` first: retry directly only if neither target environment
was created. A dependency-solving or installation failure is not automatically
retried with a different manager; partial environments are preserved for inspection.

The named installer refuses to overwrite either existing environment. Inspect
`conda env list` before reinstalling. If only the package code changed and
dependencies are unchanged, update from the repository while the main environment
is active:

```bash
conda activate ldsc-gpca
python -m pip install .
```

This updates the Python package and bundled R source, not the separately installed
R/Conda dependencies. Follow environment changes separately when upgrading.

Named environments are placed in Conda's configured environment directories.
Version records and an optional activation helper are saved under
`<repository>/.environments/named/`. No `activate.sh` is needed for normal use.

Older installations in `<repository>/.environments/main` are not renamed or moved.
Their existing activation helper still works. To obtain the new named environments,
run the installer without a custom root; this is a separate installation and uses
additional disk space.

For Docker or advanced custom-path installations, an explicit root preserves
the prefix-based mode:

```bash
bash scripts/setup_environments.sh --no-activate /absolute/path/ldsc-gpca-envs
```

That mode uses `<root>/main` and `<root>/ldsc`, and activates through
`source <root>/activate.sh` (not `conda activate ldsc-gpca`). Do not move
Conda environments after creation.

</details>

## File and command conventions

- Replace `/data`, `/references` and `/results` in examples with your paths.
  Examples illustrate layouts; the package does not supply those datasets.
- A **manifest** lists traits/files; a **variant table** lists SNPs. They have
  different headers and delimiters. All managed manifests are comma-separated
  CSVs with a header, one row per trait and unique, non-empty `traitname` values.
- Use simple identifiers such as `Trait_A`, with no whitespace, slashes or
  leading/trailing spaces. Keep them identical across manifests, LDSC estimates
  and per-trait filenames. GPCA requires at least two retained traits and preserves
  manifest row order in matrices, loadings and GWAMA inputs.
  Identifiers such as `NA`, `NaN`, `nan` and `001` are read literally, quoted
  or unquoted, in manifests and Python LDSC `p1`/`p2` columns. Leading zeros
  are preserved. Missing-value tokens apply to numeric fields; empty or
  whitespace-only identifiers are invalid.
  See [cleanup validation](../tests/CLEANUP_VALIDATION.md) for before/after checks
  of package performance changes using saved real LDSC and GenomicSEM results.
- Header spelling and case matter. Use one canonical spelling for each managed
  argument/header. Multiword managed flags use underscores; abbreviations and
  older argument spellings are rejected. Extra input columns are accepted;
  old manifest names may be annotations but cannot replace required headers.
- Quote shell paths containing spaces. Within CSVs, quote fields containing
  commas. Prefer absolute paths. Relative manifest file paths resolve beside the
  manifest for `prepare` and `genomicsem ldsc`; Python `ldsc` uses supplied VCF
  paths relative to the working directory. Command-line paths are relative to the
  working directory.
- In tables below, **Required** means no default; **Unset** means omitted.
  **Off** means a flag is absent. Enable flags such as `--validate_only` by writing
  the flag alone, without `true` or `false`. Options with values accept
  `--option value`.
- Use fresh output directories. Some commands refuse existing outputs; others
  can leave older files in place after a failure. An old CSV is not evidence that
  the latest command succeeded.

### Input formats at a glance

| Input | Delimiter / file type | Required columns or contents |
| --- | --- | --- |
| Preparation manifest | CSV, comma | `traitname,vcf_files` |
| Python LDSC VCF manifest | CSV, comma | `traitname,vcf_files,ref,population_prevalence,sample_prevalence` |
| Python LDSC reuse manifest | CSV, comma | `traitname,ref,population_prevalence,sample_prevalence` |
| GenomicSEM manifest | CSV, comma | `traitname,sample_prevalence,population_prevalence`, plus `sumstats_file` or `munged_file` when that mode needs it |
| GPCA manifest, either backend | CSV, comma | `traitname`; also `vcf_files` for automatic preparation |
| GWAMA per-trait input | TSV, tabs | Required: `SNPID,CHR,BP,EA,OA,EAF,N,Z,P`; extra columns/order accepted |
| Python pairwise LDSC results | CSV, TSV or whitespace; `.gz` accepted | Pair identifiers, estimates, SEs and h2 columns; [full schema below](#python-ldsc-results-table) |
| GenomicSEM LDSC results | Binary `.RData` | Object `LDSCoutput` containing `S,V,I,S_Stand,V_Stand` |
| Munged LDSC summary statistics | TSV, normally `.sumstats.gz` | `SNP,A1,A2,N,Z`; these are not GWAMA inputs |
| HapMap list for `prepare` | TSV, tabs | `SNP`; other columns ignored |
| HapMap reference for either LDSC workflow | Whitespace, tabs/spaces | `SNP,A1,A2` |

Commas in column lists above separate **names**, not necessarily file fields.
Use each row's stated delimiter. Each module below explains its conditional
columns, numerical conventions and outputs.

## Run the entire analysis in one command

`ldsc-gpca pipeline` joins the existing analysis stages:

```text
VCF preparation, if needed
          ↓
Python LDSC or GenomicSEM munging/LDSC
          ↓
LDSC QC → PCA → chromosome/whole-genome PC1 GWAMA
          ↓
combined results + selected-column summary + audit files
```

Choose `--ldsc_backend python` or `--ldsc_backend genomicsem`. All scientific
calculations run through the existing modules; the original GWAMA function and
variant headers are unchanged. Stages run sequentially, with each stage using
its supported internal parallelism. This simplifies execution; it does not make
GenomicSEM LDSC parallel or make the two implementations statistically
identical.

After updating the clone, install the new command into your main environment:

```bash
conda activate ldsc-gpca
python -m pip install .
ldsc-gpca --version
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
```

The following examples assume existing input/reference files. Replace the paths,
trait names and export metadata with values appropriate for your analysis.

### One manifest for VCF-based analysis

Create `/data/pipeline_traits.csv`, **comma-separated**, with a header:

```csv
traitname,vcf_files,sample_prevalence,population_prevalence
Trait_A,/data/Trait_A.vcf.gz,NA,NA
Trait_B,/data/Trait_B.vcf.gz,NA,NA
Trait_C,/data/Trait_C.vcf.gz,NA,NA
Trait_D,/data/Trait_D.vcf.gz,NA,NA
```

At least two unique traits are required. Relative file paths resolve beside this
manifest for both pipeline backends. For Python LDSC, pipeline sets `ref=yes`
for every trait to obtain all pairs/self-pairs. You may include `ref`, but every
value must be `yes`; target-only selection cannot provide complete GPCA input.
Both prevalence headers are required; the existing backend's prevalence/N rules
still apply.

| VCF route | Fields and limitations |
| --- | --- |
| Python LDSC plus automatic GWAMA preparation | Each VCF must satisfy **both** existing schemas: FORMAT `AF,ES,SE,LP,NEF` for preparation, plus FORMAT `SI,EZ` and INFO `AF,EUR` for Python extraction. Binary population-prevalence traits additionally need FORMAT `NC,NCO`. |
| GenomicSEM with `--vcf_input` | FORMAT `AF,ES,SE,LP,SI`, plus NEF unless manifest N overrides it. Converts VCFs to raw tables with INFO for GenomicSEM filtering. Binary traits require an appropriate explicit manifest N and both prevalences. Separate GWAMA preparation still needs NEF. See the [VCF guide](GENOMICSEM_LDSC.md#start-from-gwas-vcfs). |
| Older GenomicSEM VCF fallback without `--vcf_input` or LDSC tables | Quantitative traits only, with both prevalences blank/NA. Preparation needs FORMAT `AF,ES,SE,LP,NEF` and writes HapMap-selected raw tables. Those tables use `N=NEF` and have **no INFO column**; GenomicSEM INFO filtering cannot act on an absent field. |
| GenomicSEM from separately prepared LDSC tables | Add `sumstats_file` paths to the manifest, or use existing munged inputs below. These routes support quantitative and binary traits with suitable N and prevalence settings. |

VCFs must have one GWAS sample and consistent ancestry/build/alleles. The pipeline
does not harmonise alleles, infer the phenotype type or infer a defensible
sample-size convention. Preparation's default GPCA IDs use `CHR_POS_REF_ALT`;
its optional munging IDs use the VCF ID and must match the HapMap reference.
`--hm3` is a plain whitespace-separated allele reference with `SNP,A1,A2`
headers. Pipeline creates a SNP-only TSV for preparation when needed and retains
the original allele reference for munging.

### Complete Python LDSC pipeline

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
  --input /data/pipeline_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/cluster1_python_pipeline \
  --chisq_max 80 \
  --splitby_chr split \
  --n_cores 22 \
  --prepare_workers 4 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

`GWAMA_INFO` must be assigned a justified export value before running, as explained
in [the GWAMA export policy](#all-automatic-export-options). To deliberately
replace the summary's N_eff with a constant, add `--gwama_output_n_eff VALUE`;
otherwise reported N_eff is preserved. Neither override changes the analysis.

### Complete GenomicSEM pipeline

For the quantitative VCF manifest above, change the backend and output directory:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/pipeline_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/cluster1_genomicsem_pipeline \
  --chisq_max 80 \
  --splitby_chr split \
  --n_cores 22 \
  --prepare_workers 4 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

If `sumstats_file` is populated for every trait, those files supply GenomicSEM LDSC
munging instead of generating munging tables from the VCFs. Example manifest:

```csv
traitname,vcf_files,sumstats_file,sample_prevalence,population_prevalence
Trait_A,/data/Trait_A.vcf.gz,/data/Trait_A_ldsc.tsv,NA,NA
Trait_B,/data/Trait_B.vcf.gz,/data/Trait_B_ldsc.tsv,NA,NA
Trait_C,/data/Trait_C.vcf.gz,/data/Trait_C_ldsc.tsv,NA,NA
Trait_D,/data/Trait_D.vcf.gz,/data/Trait_D_ldsc.tsv,0.2,0.05
```

The example prevalence values illustrate the format, not a recommended setting.
Follow the [GenomicSEM unmunged-table requirements](#genomicsem-ldsc-inputs-choose-one-mode)
for columns, N and allele conventions. Partly supplied `sumstats_file` paths are
rejected; the pipeline does not guess a different input mode for each trait.

The GenomicSEM backend still performs its existing two LDSC passes. Pipeline reads
the retained-trait audit and automatically writes a clean `traitname`-only GPCA
manifest in original manifest order. No manual `Selected_Traits.csv` conversion
is needed in pipeline mode. The GenomicSEM default `--invalid_h2_action drop` can
remove invalid-h2 traits before GPCA; use `error` to stop instead.

### Starting with munged and prepared files

You can skip repeated preparation/munging. Use a CSV with
`traitname,sample_prevalence,population_prevalence`; VCF paths are unnecessary
when an existing GWAMA input folder is supplied. Python's munged directory must
retain its required prevalence sidecars. For GenomicSEM reuse, produce files with
`prepare --mode ldsc --munge_backend genomicsem` (or `--mode both`), using the
same `--hm3` allele reference for all traits. Do not directly reuse Python-munged
files; matching HapMap SNPs does not guarantee GenomicSEM-compatible allele orientation.

Python:

```bash
ldsc-gpca pipeline \
  --ldsc_backend python \
  --input /data/reuse_traits.csv \
  --ldsc_only \
  --munged_dir /data/python_munged \
  --gpca_input_folder /data/prepared/gpca_inputs \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/cluster1_python_reuse_pipeline \
  --chisq_max 80 \
  --n_cores 22 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

GenomicSEM:

```bash
ldsc-gpca pipeline \
  --ldsc_backend genomicsem \
  --input /data/reuse_traits.csv \
  --munged_dir /data/genomicsem_munged \
  --gpca_input_folder /data/prepared/gpca_inputs \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/cluster1_genomicsem_reuse_pipeline \
  --chisq_max 80 \
  --n_cores 22 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

GenomicSEM also accepts `--munged_input` with manifest `munged_file` paths,
mutually exclusive with `--munged_dir`. Python pipeline reuse requires
`--ldsc_only` and `--munged_dir` together. If you already have LDSC estimates,
use the [separate GPCA command](#run-gpca-and-gwama) to avoid rerunning LDSC.

### Chi-square cutoff and pipeline options

`--chisq_max 80` removes variants with **`Z² > 80` from LDSC inputs**, independently
for each trait. GenomicSEM applies its per-trait filter after merging LD
references; managed Python filters copies of munged files before regression.
The explicit threshold is forwarded to both GenomicSEM LDSC passes. It is **not** a
PCA significance filter and does **not** remove variants from the GPCA/GWAMA input
files. Equivalent thresholds do not guarantee identical SNP sets or estimates
because other backend filters/conventions can differ. See the
[pinned GenomicSEM implementation](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R).

| Pipeline option | Default | Meaning |
| --- | --- | --- |
| `--ldsc_backend` | `python` | `python` or `genomicsem`; chooses LDSC and its matching GPCA reader. |
| `--input` | Required | CSV described above; at least two traits, exact names/order. |
| `--outdir` | Required | A directory that **does not yet exist**; completed/partial runs are never overwritten. |
| `--chisq_max` | Unset | Python: positive integer or `auto`. GenomicSEM: positive finite number; omit for automatic mode. When omitted, Python's extra filter is disabled and GenomicSEM uses its automatic rule; the audit records which policy applied. |
| `--rg_normalization` | `pair` | Python backend only: `pair` uses original `rg`; `trait_wide` uses the newly computed `rg_trait_wide` column for PCA/GWAMA. |
| `--n_cores` | Python `5`; GenomicSEM `1` | Positive worker count passed to LDSC and GPCA/export. GenomicSEM LDSC itself remains sequential. `0` is not accepted by pipeline. |
| `--prepare_workers` | `4` | Positive trait-preparation worker count, separate from LDSC/GWAMA workers. |
| `--gpca_input_folder` | Unset | Reuse existing nine-column GWAMA tables. Otherwise prepare from VCFs. |
| `--validate_only` | Off | Still run upstream LDSC, then QC/PCA; skip GWAMA and export. Preparation still runs if needed to generate GenomicSEM LDSC munging tables or explicitly requested munging exports. |
| `--dataset_id` | Pipeline output-directory basename | Final filename prefix under `gpca/`. |
| `--splitby_chr` | `split` | `split` needs all chromosomes 1–22 per trait; `nosplit` uses one whole-genome file per trait. |
| `--gwama_output_info` | Unset | Required up front for a full pipeline using bundled GWAMA. A custom source supplying INFO may omit it. |
| `--help` | Off | Backend-specific complete option list; select the backend before requesting help. |
| `--prepare_help`, `--postprocess_help` | Off | Existing preparation/export help. |

The remaining options reuse their existing names, defaults and meanings:

- [Python LDSC options](#all-python-ldsc-options) or
  [GenomicSEM LDSC options](#all-genomicsem-ldsc-options), according to the backend.
  `--hm3` always needs `SNP,A1,A2` in pipeline mode. It is unnecessary when neither
  munging nor optional munging-table preparation is needed. `--rscript` configures
  GenomicSEM LDSC only; the GPCA stage uses `Rscript` on PATH. `--restart` belongs to
  standalone Python `ldsc` and is rejected by `pipeline`, which requires a fresh
  output directory.
- [Shared GPCA settings](#all-shared-gpca-analysis-options), except
  `--ldsc_results`: pipeline generates that file and passes it automatically.
  Defaults remain strict GPCA failure handling and correlation PCA.
- [Python-table GPCA options](#additional-options-for-python-table-gpca-only)
  with the Python backend, except the cosmetic `--color` option.
- [Automatic preparation options](#automatic-vcf-preparation) and
  [export options](#all-automatic-export-options). Optional `--write_munge_inputs`
  also works when a GWAMA folder is supplied, provided VCFs and HM3 are available.

Backend-inapplicable options, abbreviations and old spellings are rejected.
With Python `--result_failure_action drop_traits`, pipeline reads
`LDSC_Retained_Traits.csv` and passes only that retained subset to GPCA, in the
original manifest order. GenomicSEM uses its `Selected_Traits.csv` in the same
way. Missing or invalid retained-trait audits stop the pipeline. If Python
`--result_failure_action report` produces only diagnostic results, pipeline
stops before GPCA because a valid `ldsc_results.csv` is required.
`--failed_ldsc_action drop_traits` and `--allow_missing_traits` remain explicit
GPCA choices; removals change the analyzed component. Python munging failures
can leave fewer LDSC traits, and GPCA will stop on missing selected traits unless
its missing-trait policy permits removal.

### Pipeline outputs and failure handling

```text
<outdir>/
  Pipeline_Run_Status.json
  manifests/
    input_traits.csv
    ldsc_traits.csv
    gpca_traits.csv
    hm3_snps.tsv                 # only when preparation needs SNP selection
  prepare/                      # only when preparation runs
    gpca_inputs/
    munge_inputs/               # when requested/needed for GenomicSEM VCF input
    ... preparation QC reports
  ldsc/
    ldsc_results.csv            # Python backend
    genomicsem_LDSC.RData        # GenomicSEM backend instead
    ... backend audits/intermediates
  gpca/
    ... matrix, PC and GWAMA reports
    <dataset_id>_GWAMA_combined_results.txt.gz
    <dataset_id>_postprocess.json
    harmonisation_input/<dataset_id>_GPCA_inputs.txt.gz
```

The pipeline status JSON records version, backend, cutoff policy, original trait
order, the manifest passed to GPCA, removals before GPCA, exact stage commands,
elapsed times and completion/failure. Further GPCA removals are in its own audits.
Failure stops subsequent stages and preserves completed outputs. The pipeline
does not automatically resume; use existing individual commands to reuse those
outputs, or choose a fresh directory for another full run. Final export remains
the existing concurrent Polars/pigz implementation.

For automatic Python cutoffs, use `--chisq_max auto`. To use trait-wide
normalization downstream, add `--rg_normalization trait_wide`. These are
independent settings: the first selects the LDSC SNP filter, and the second
selects the correlation values used by PCA/GWAMA.

See [pipeline validation](../tests/PIPELINE_VALIDATION.md) for fixture coverage and
what has not been tested on real data.


## Prepare per-trait GPCA inputs

`ldsc-gpca prepare` accepts GWAS summary-statistics VCFs. Its default still writes
GPCA/GWAMA tables. New explicit modes can prepare shared raw LDSC tables and run
Python or GenomicSEM munging. No preparation mode runs LDSC regression, PCA or
GWAMA. See the [complete preparation guide](PREPARE.md).

### Preparation input

Every mode needs a CSV with `traitname,vcf_files` and one GWAS sample per VCF.
GPCA output uses FORMAT `AF,ES,SE,LP,NEF`. Shared LDSC output also needs `SI`;
a manifest N can replace NEF for LDSC only. Binary LDSC preparation requires
explicit total N and both prevalences. The [input examples and N rules](PREPARE.md#files-you-need)
explain each route. LDSC preparation requires a tab-separated `--hm3` reference;
munging also needs its `A1,A2` columns.

### Preparation commands

| Output needed | Options after `ldsc-gpca prepare --input CSV --outdir DIR` |
| --- | --- |
| GPCA only (default) | `--splitby_chr nosplit` for whole-genome tables |
| Shared raw LDSC tables | `--mode ldsc --raw_only --hm3 FILE` |
| Raw and Python-munged LDSC files | `--mode ldsc --hm3 FILE` |
| Raw and GenomicSEM-munged LDSC files | `--mode ldsc --hm3 FILE --munge_backend genomicsem` |
| GPCA plus shared raw tables | `--mode both --raw_only --hm3 FILE --splitby_chr nosplit` |
| GPCA, raw tables and munged files | `--mode both --hm3 FILE --splitby_chr nosplit`; add `--munge_backend genomicsem` for GenomicSEM R munging |

The default GPCA layout is `split` and requires valid variants on all chromosomes
1–22 per trait. `nosplit` writes one autosomal GPCA file per trait. LDSC files are
always unsplit. Use the guide's [complete commands](PREPARE.md#basic-command).

### All preparation options

See the [full options, defaults and applicability table](PREPARE.md#all-options),
or run `ldsc-gpca prepare --help`. `--info_filter` and `--maf_filter` apply only
when munging runs. Munging/runtime options are rejected with `--raw_only` and
GPCA-only mode. The existing `--write_munge_inputs` option retains its
[legacy schema and behavior](PREPARE.md#also-write-raw-tables-for-later-munging).

### Preparation outputs and filtering

| Folder/file | What it contains |
| --- | --- |
| `gpca_inputs/` | GPCA/GWAMA TSVs: `SNPID,CHR,BP,EA,OA,EAF,N,Z,P`. |
| `munge_inputs/` in new LDSC modes | Raw TSVs: `SNP,CHR,BP,A1,A2,EAF,BETA,SE,P,N,INFO`. These still need munging. |
| `munged/` | LDSC-ready `{traitname}.sumstats.gz` files with `SNP,A1,A2,N,Z`, plus provenance sidecars. |
| `Prepared_LDSC_Manifest.csv` | Reusable paths, N and prevalence metadata for new LDSC modes. |
| QC, settings and worker reports | Counts, original affected VCF rows, settings, errors and retries. |

See [file names, column meanings and reuse commands](PREPARE.md#outputs-and-next-step).
GPCA and LDSC QC are independent; HapMap selection and LDSC INFO filters do not
restrict GPCA rows. Basic preparation removes invalid rows, retains autosomes
1–22 and resolves duplicate IDs. Selected munging software then applies its own
filtering and allele rules. Source VCFs remain unchanged. Failed preparation
publishes audits but no tables; failed munging retains prepared raw/GPCA tables
and logs, without publishing `munged/`.

## Run Python LDSC

`ldsc-gpca ldsc --ldsc_backend python` extracts VCFs, munges them with CBIIT LDSC, runs requested
pairwise LDSC batches and compiles numerical estimates into `ldsc_results.csv`.
It can also reuse existing munged files. Python is the default when
`--ldsc_backend` is omitted. This option selects both munging and regression;
`--munge_backend` belongs only to standalone `prepare`.

### Python LDSC input

Save `/data/python_ldsc_traits.csv` as a comma-separated CSV:

```csv
traitname,vcf_files,ref,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,yes,NA,NA
Trait_B,/data/Trait_B.vcf.gz,yes,NA,NA
Trait_C,/data/Trait_C.vcf.gz,yes,NA,NA
Trait_D,/data/Trait_D.vcf.gz,yes,NA,NA
```

| Manifest column | Requirement / meaning |
| --- | --- |
| `traitname` | Unique, non-empty trait identifier; becomes the munged filename prefix. |
| `vcf_files` | Required when extracting VCFs; use absolute paths. May be omitted with `--ldsc_only`. |
| `ref` | Use literal `yes` or `no`. `yes` requests that trait against every listed trait, including itself; `no` is target-only. **Use `yes` for every trait for GPCA coverage.** |
| `population_prevalence` | Header always required. Blank/`NA` means no liability conversion; otherwise strictly between 0 and 1. |
| `sample_prevalence` | Header always required. For a population-prevalence trait, supply a fraction strictly between 0 and 1 or leave blank/`NA` for derivation from case counts. Ignored when population prevalence is absent. |

**This VCF schema differs from `prepare`.** The extraction commands expect
INFO fields `AF,EUR` and FORMAT fields `SI,AF,EZ,LP,NEF`. Population-prevalence
traits additionally require FORMAT `NC,NCO` (case/control counts).
`SI` supplies imputation quality, `EZ` is the signed Z statistic, and
`LP=-log10(P)`. Frequencies and effect statistics must follow the source allele
convention. This fixed EUR-field-specific workflow does not provide flags to
map a different ancestry or arbitrary VCF field names. Prepare one GWAS per VCF;
multi-sample extraction is not a supported general input layout.

The AF-difference filter excludes records with missing `INFO/AF` or `INFO/EUR`,
including absent record tags and `.` values. It never substitutes `FORMAT/AF`;
missing INFO header declarations stop extraction. Each trait writes
`munge_input/{traitname}_AF_Filter_QC.csv`, counting records after optional MHC
exclusion and missing `INFO/AF`, missing `INFO/EUR`, and missing either field
(the union, without double-counting). Counts are collected in the same VCF pass,
before the AF-difference, SI, MAF and palindromic filters; an audit is published
only after successful extraction. `--ldsc_only` reuses existing munged inputs;
rerun extraction and munging to apply this policy to previous analyses.

| Population prevalence | Sample size used | Sample prevalence |
| --- | --- | --- |
| Blank/`NA` | `NEF` | Not used |
| Supplied | `NC + NCO` | Supplied fraction, or median per-variant `NC/(NC+NCO)` |

The derived median is an implementation choice; verify its suitability for
meta-analysis or variable per-SNP sample sizes. An unconverted estimate based
on binary effective N must not automatically be interpreted as ordinary
observed-scale heritability.

Reference inputs:

- `--hm3`: whitespace-separated text with `SNP,A1,A2` headers. Unlike preparation's
  SNP-only list, this is passed to LDSC allele matching.
- `--ld_ref`: directory containing `1.l2.ldscore.gz` through `22.l2.ldscore.gz`
  and associated M files, including `<CHR>.l2.M_5_50` for the default regression.
- `--ld_weights`: optional separate directory with chromosome `.l2.ldscore.gz`
  regression weights. Default: the `--ld_ref` directory. Check that your reference
  is suitable for both roles before relying on this default.

### Python LDSC commands

From VCFs:

```bash
ldsc-gpca ldsc \
  --input /data/python_ldsc_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/python_ldsc \
  --n_cores 22
```

From existing munged files, use a manifest with
`traitname,ref,population_prevalence,sample_prevalence`:

```bash
ldsc-gpca ldsc \
  --input /data/python_ldsc_reuse_traits.csv \
  --ldsc_only \
  --munged_dir /data/munged \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/python_ldsc_reused \
  --n_cores 22
```

The reuse directory must contain `{traitname}.sumstats.gz` (tab-separated
`SNP,A1,A2,N,Z` columns) and applicable `{traitname}.prevalence.json` sidecars.
For total-N traits, matching sidecars/file hashes are required. Legacy NEF files
without sidecars are accepted with provenance warnings. Reuse does not reapply
extraction/munging filters; recorded mismatches fail and unknown settings warn.
By default it reruns all requested LDSC batches. Add `--restart` to reuse
verified completed batches in the same output directory.

### Restarting interrupted LDSC runs

Add `--restart` to the same command and keep the same `--outdir`:

```bash
ldsc-gpca ldsc \
  --input "${python_csv}" \
  --ldsc_only \
  --munged_dir "${munge_dir}/" \
  --ld_ref "${ld_ref}/" \
  --outdir "${out_folder}/" \
  --n_cores 90 \
  --ldsc_retries 1 \
  --chisq_max 80 \
  --result_failure_action report \
  --restart
```

Every new CLI run writes an atomic `.results.csv.checkpoint.json` after each
batch exits successfully and its numerical table has all requested comparisons.
`--restart` checks SHA-256 fingerprints of effective munged contents, LD scores
and M files, regression weights, the batch command/prevalences, filter settings,
the Python LDSC/exporter source and runtime versions, and the saved result file.
It reuses matching batches and reruns missing, interrupted, corrupt or changed
batches. Checkpoints are independent, so an interrupted job does not invalidate
other completed jobs. Results are recompiled and QC is reapplied in manifest order.

Completed numerical estimation failures (for example, missing rg from negative
h2) remain recorded and can be reused; restarting does not make those estimates
valid. `--result_failure_action drop_traits` can exclude unusable traits during
LDSC collection. GenomicPCA also retains its separate trait-removal option.
Changing the collection policy or
retry count does not invalidate otherwise matching regression results.

**Older outputs without these checkpoints rerun once**, even when their names
and `LDSC_Runtime.json` appear to match: those files cannot prove successful
completion with the current inputs. Keep checkpoints beside the batch CSVs.
Changing the worker count can change the batch layout and cause recomputation;
keep it unchanged for the most predictable restart. Extra manifest annotation
columns do not invalidate a batch. Gzip header timestamps are ignored when
hashing effective sumstats, so repeated chi-square filtering with identical
contents can reuse results.

Restart applies to pairwise LDSC. `--ldsc_only` skips extraction and munging;
chi-square filtering, content verification and final compilation still run.
Content verification reads each shared input once per invocation and adds I/O.
Keep inputs/runtime unchanged during a run and use one active command per output
directory. The batch status report is updated as each job finishes, with
`completed`, `reused`, `execution_failed` or `pending`, plus a restart reason.

### Parallel worker failures and retries

Every package-managed parallel stage uses **two attempts by default**: the initial
attempt and one retry if it fails. Successful jobs run once. Only failed or
unconfirmed jobs are retried. If a job still fails, the analysis stops with the
stage, trait/chromosome/job, attempt count and error; downstream stages do not run.
The current stage collects worker outcomes before reporting failure, preserving
completed outputs. This is not an immediate termination of other running jobs.

| Stage | Attempt audit |
| --- | --- |
| Python LDSC extraction | `LDSC_Extraction_Worker_Attempts.csv` |
| Python LDSC munging | `LDSC_Munging_Worker_Attempts.csv` |
| Chi-square filtering | `LDSC_ChiSquare_Worker_Attempts.csv` |
| Pairwise Python LDSC | `ldsc_results/LDSC_Batch_Status.csv`, including `Attempts` |
| VCF preparation for GPCA | `Preparation_Worker_Attempts.csv` |
| GenomicSEM munging | `GenomicSEM_Munging_Worker_Attempts.csv` |
| GWAMA, either backend | `GWAMA_Worker_Attempts.csv` and `GWAMA_Run_Status.csv` |
| GWAMA export file reads | `GWAMA_Export_Worker_Attempts.csv` |

Python LDSC retains its existing `--ldsc_retries` override: default `1` means two
total attempts; `0` explicitly disables its retries. Other managed stages use
one retry. Retries in those stages run one job at a time to reduce concurrent
memory demand. GWAMA and GenomicSEM munging use fresh child processes on Unix/macOS,
including serial retries; Windows runs these R jobs sequentially. A terminated
main/controller process cannot retry itself. No upstream estimator is modified.

`NULL`, `try-error` and malformed R worker results become explicit failures,
rather than breaking status reporting. GWAMA requires fresh non-empty output;
an old file cannot establish success. Normal runs retain their matrix, PC1 and
GWAMA calculations. Preflight validation errors stop before workers are launched.

Execution failures are fatal even with `--result_failure_action drop_traits`.
That option still handles invalid estimates from completed LDSC commands, but
does not discard traits to hide crashed/exhausted jobs. Failed munging no longer
silently removes traits. See [validation](../tests/WORKER_RETRY_VALIDATION.md).

### All Python LDSC options

| Option | Default | Meaning / accepted values |
| --- | --- | --- |
| `--ldsc_backend CHOICE` | `python` | Select `python` for these options; `genomicsem` selects the GenomicSEM parser and workflow. No separate `--munge_backend` is accepted. |
| `--input CSV` | Required | Trait manifest described above. |
| `--outdir DIRECTORY` | Required | Results, intermediate files and logs. Reuse the same directory with `--restart`. |
| `--ld_ref DIRECTORY` | Required | Chromosome LD-score reference directory. |
| `--ld_weights DIRECTORY` | Same as `--ld_ref` | Regression weights directory. |
| `--hm3 FILE` | Unset | Required unless `--ldsc_only`; HapMap `SNP,A1,A2` reference. |
| `--n_cores INTEGER` | `5` | Concurrent local extraction, munging and LDSC workers; integer ≥1. |
| `--result_failure_action {error,report,drop_traits}` | `error` | Stop on numerical failures (`error`), preserve failures for inspection (`report`), or exclude failed traits and publish a complete retained subset (`drop_traits`). All modes retain diagnostics. |
| `--ldsc_only` | Off | Skip extraction/munging and reuse munged files. |
| `--restart` | Off | Reuse completed batches with matching input/parameter/runtime checkpoints; rerun unverified or changed batches. |
| `--munged_dir DIRECTORY` | `<outdir>/ldsc_input` | Existing munged input directory; requires `--ldsc_only`. Rejected before output creation without that flag. Leave unset for fresh VCF runs. |
| `--ldsc_retries INTEGER` | `1` | Additional attempts per failed LDSC command or malformed/incomplete export; ≥0. One means two total attempts. Numerical estimation failures are not retried. |
| `--chisq_max INTEGER\|auto` | Unset; disabled | Positive integer: fixed cutoff. `auto`: per-trait `max(80, 0.001 * max(N))` after complete-row matching to reference and weight LD-score SNPs. Keep `Z² <= cutoff` before LDSC in both modes. |
| `--exclude_mhc` | Off | Exclude the configured MHC interval during VCF extraction. |
| `--mhc_chr STRING` | `6` | Chromosome label used for MHC exclusion; match the VCF's naming. |
| `--mhc_start INTEGER` | `25000000` | Inclusive MHC start; positive, relevant with `--exclude_mhc`. |
| `--mhc_end INTEGER` | `35000000` | Inclusive MHC end; must be ≥ start. Verify genome build. |
| `--info_min FLOAT` | `0.7` | Inclusive minimum FORMAT/SI during extraction; `[0,1]`. |
| `--maf_min FLOAT` | `0.01` | Extraction keeps FORMAT/AF between this value and `1-value`, inclusively; `(0,0.5)`. |
| `--munge_maf_min FLOAT` | `0.005` | Additional munging MAF threshold, strictly greater than this value; `[0,0.5)`. |
| `--max_af_difference FLOAT` | `0.2` | Maximum absolute difference between INFO/AF and INFO/EUR; `[0,1]`. Missing either INFO field is excluded and counted; no FORMAT/AF fallback. |
| `--remove_palindrome` | Off | During extraction, remove A/T and C/G variants within the configured AF interval. LDSC munging subsequently removes all palindromic SNPs regardless. |
| `--paliandromaf_lower FLOAT` | `0.45` | Inclusive lower AF bound for extraction's palindromic removal. Use this exact option spelling. |
| `--paliandromaf_upper FLOAT` | `0.55` | Inclusive upper AF bound; `0 <= lower <= upper <= 1`. |
| `--conda_executable EXECUTABLE` | `CONDA_EXE`, otherwise `conda` | Conda executable used to launch child LDSC. |
| `--ldsc_env NAME` | Saved prefix if configured; otherwise `ldsc-cbiit` | Explicit child environment name overrides the saved prefix; mutually exclusive with `--ldsc_env_prefix`. Setup configures the child automatically. |
| `--ldsc_env_prefix DIRECTORY` | `LDSC_GPCA_LDSC_PREFIX` if set | Explicit child environment prefix; mutually exclusive with `--ldsc_env`. |
| `--bcftools EXECUTABLE` | `bcftools` | Local program/path; unnecessary with `--ldsc_only`. |
| `--help` | Off | Show usage and exit. |

Managed munging preserves native computed Z and N with **17 significant digits**
(`%.17g`) when it first writes `.sumstats.gz`, before filtering or regression.
It calls native `munge_sumstats(args, p=False)` and writes the returned table once;
there is no additional file-reading pass. Native QC, allele handling, Z calculation,
column order and missing values are preserved. No installed LDSC source is edited.
The [native writer](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py#L695-L718)
otherwise rounds to three decimal places. For example, Z=8.9444 has Z²=80.00229136
and is excluded at 80, whereas rounding to 8.944 would retain it. More digits can
increase file size and writing time; this preserves computed floating-point
precision, not precision already lost in the source GWAS.

`LDSC_Runtime.json`, each munging log and its `.prevalence.json` sidecar record
the managed format. Raw `ldsc-gpca munge_sumstats.py` remains an unchanged upstream
passthrough. `--ldsc_only` reuses existing numbers and does not restore lost digits.
To replace previously rounded results, rerun from original inputs **without
`--ldsc_only`**, then rerun PCA/GWAMA. Restart checkpoints fingerprint both the
munging wrapper and native source, as well as actual regression inputs, so old
checkpoints are invalidated by this change. A restart still cannot improve the
precision of a reused rounded file. See [validation](../tests/MUNGE_PRECISION_VALIDATION.md).

Choose chi-square filtering explicitly for `ldsc-gpca ldsc`:

| Setting | Behavior |
| --- | --- |
| Omit `--chisq_max` | No additional chi-square filtering (unchanged default). |
| `--chisq_max 80` | Use cutoff 80 for every trait; any positive integer is accepted. |
| `--chisq_max auto` | Calculate a separate cutoff for each trait using GenomicSEM's automatic rule. |

For `auto`, first omit rows with missing fields and match SNPs to **both** the
reference and regression-weight LD-score files. Take the largest `N` among those
rows, multiply by 0.001, and use the larger of that result and 80. For example,
maximum matched N of 50,000 gives cutoff 80; 200,000 gives 200; 500,250 gives
500.25. Automatic cutoffs are not rounded. `N` is the value already in the munged
file, including effective N when that convention was used. The cutoff is computed
before excluding large Z² values and before intersecting trait pairs.
[GenomicSEM source](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R#L172-L195).

Both settings keep `Z² <= cutoff` and write separate filtered copies, preserving
source files. LD matching determines the automatic cutoff; Python LDSC still
performs its usual LD matching and missing-row removal on the copies. The summary
`LDSC_ChiSquare_Filter_Summary.csv` records each trait's `threshold_mode`, numeric
`chisq_max`, `maximum_matched_n`, and `complete_ld_matched_rows`; the last two are
blank for fixed cutoffs. Reference SNPs are loaded once per run, and automatic N
selection uses a memory-bounded preliminary pass over each trait file.

Filtered copies and exclusion files use `pigz` when installed, with Python gzip
as the fallback. Both use compression level 9 and preserve the uncompressed
row text. Compressor failures follow the existing two-attempt worker policy;
failed compression is never published as a completed filtered file. Compression
threads are divided between the two output streams and active trait jobs, with
at least one compressor worker per stream. The final combined exclusion file
can use all `--n_cores` workers. Final GWAMA export's `--gzip_level` setting does
not change filtering compression.

Standard tab-separated LD-score SNP columns are read with Polars; files with
mixed whitespace, quoting or empty tab fields retain the pandas whitespace
reader. Numerical LDSC result CSVs also use Polars, followed by the existing
numeric conversion and QC. This preserves numeric precision and trait names
such as `001` and `NA`. These paths support the pinned Polars **0.20.31** and
require no dependency upgrades or PyArrow. See the
[I/O validation](../tests/IO_ACCELERATION_VALIDATION.md) for measured timings and scope.

The option requires a value: `auto` or an integer greater than zero. Zero,
negative values, decimals (including `80.0`), and a bare `--chisq_max` are rejected.
This replaces the previous floating-point CLI cutoff syntax; write `80` instead
of `80.0`. Changing filtering requires rerunning LDSC; `--restart` reuses batches
only when the selected mode, inputs, references, and other settings match.
Previously filtered inputs cannot recover SNPs already removed: use the original
munged files when comparing cutoffs. Matching this filter alone does not remove
other differences between Python LDSC and GenomicSEM.

The managed option uses per-trait filtering and is not forwarded to raw Python
LDSC's cross-product `--chisq-max`; see [raw commands](#raw-ldsc-commands).
Missing Z placeholders (blank, `NA`, `nan`, `NaN`, or `.`) pass through unchanged
for Python LDSC to discard. They are counted as `variants_missing_z` in the filter
summary; `variants_after` includes these rows, so finite-Z survivors equal
`variants_after - variants_missing_z`. Malformed/non-finite non-missing values
still fail, as does a trait with no finite Z surviving the threshold.
The `prepare --p_min` setting does not apply to Python LDSC extraction. Extreme
LP can underflow to P=0 here and be removed by munging; inspect its logs.

### Python LDSC outputs

| Output under `--outdir` | Purpose |
| --- | --- |
| `ldsc_results.csv` | Compiled estimates; in `drop_traits` mode contains only the retained complete subset. Full GPCA matrix validation is still required. |
| `ldsc_results_diagnostic.csv` | All compiled estimates, including missing values; suitable for explicit downstream failure handling, not automatic approval for analysis. |
| `LDSC_Pair_Status.csv`, `LDSC_Trait_Status.csv` | Numerical failures/warnings, exact trait names, manifest order and source file/row provenance. |
| `LDSC_Compilation_Status.csv` | `compiled`, `compiled_with_trait_exclusions`, `estimation_failures`, `insufficient_traits`, `structural_failure`, or `incomplete`; action, result path and selection counts when applicable. |
| `LDSC_Retained_Traits.csv` | In `drop_traits` mode, the retained manifest in original order, with required headers and source annotations. Use this manifest for downstream genomicPCA. |
| `LDSC_Dropped_Traits.csv` | Excluded traits, reasons, exclusion order, failed-pair counts and self-h2 diagnostics. Written even when fewer than two traits remain. |
| `ldsc_results/LDSC_Batch_Status.csv` | Completed, reused, failed and pending jobs; restart reasons and execution errors. |
| `ldsc_results/*.results.csv.checkpoint.json` | Per-batch completion, request provenance and numerical-output checksum, used by `--restart`. |
| `ldsc_results/` | Per-batch `.results.csv` numerical exports and readable logs. |
| `ldsc_input/` | Newly munged `{traitname}.sumstats.gz` files and prevalence sidecars in the default VCF workflow. |
| `munge_input/` | Intermediate extracted tables. |
| `LDSC_Runtime.json` | Runtime, reference paths, filters and export contract. |
| `LDSC_Trait_Prevalence_Metadata.csv` | N/prevalence choices and derivation. |
| `execution_errors.log` | Command failures, when present. |
| `ldsc_input_chisq_filtered/` | Filtered copies when `--chisq_max` is supplied. |
| `LDSC_ChiSquare_Filter_Summary.csv` | Before/removed/retained row counts and a separate `variants_missing_z` count; retained rows include missing-Z placeholders. |
| `LDSC_ChiSquare_Excluded_Variants.tsv.gz` | Excluded `gwas_name,SNP,Z,CHISQ` records; this existing audit header is unchanged. |

Managed LDSC exports each native result table numerically with `%.17g` formatting
before readable-log rounding. Compilation reads each numerical CSV once; it
does not reconstruct estimates from logs. Missing/malformed exports fail.
Numerical estimation failures follow the explicit policy below.
Extraction and munging jobs retry once; exhausted failures stop the run.
Exhausted LDSC command retries stop before compilation in every mode, including
`drop_traits`. Completed batch files and checkpoints are retained for restart.
Malformed result exports, conflicting estimates, runtime setup failures and
input/provenance errors still stop. Check retained traits and
logs before using the results: compilation is not the full GPCA QC check.

### Handling failed LDSC estimates in future runs

To remove failed traits within the **LDSC command** and continue with the
remaining set, use:

```bash
ldsc-gpca ldsc \
  --input "${python_csv}" \
  --ldsc_only \
  --munged_dir "${munge_dir}/" \
  --ld_ref "${ld_ref}/" \
  --outdir "${out_folder}/" \
  --n_cores 50 \
  --ldsc_retries 1 \
  --chisq_max 80 \
  --restart \
  --result_failure_action drop_traits
```

Set `ref=yes` for every trait so each has a self-pair. The command removes
traits with missing/invalid self estimates first, then resolves remaining
unusable pairs by highest failed-pair count, lower self-h2 Z and later manifest
position, matching the existing R selection policy. This deterministic heuristic
does not guarantee the largest possible subset. Low positive h2/SE or finite rg
outside [-1,1] are warnings, not automatic reasons to remove a trait.

At least two traits, usable self estimates and every unordered pair must remain.
One valid orientation can cover an absent opposite orientation. Conflicting
observed orientations remain fatal, using the existing absolute duplicate
tolerances (0.001 for rg/SE/p/intercepts, 0.01 for z, with a 1e-12 comparison
epsilon). An exhausted execution failure stops the run and is recorded in
`LDSC_Batch_Status.csv`; it is not treated as a reason for biological trait removal.
The original observed estimates stay in `ldsc_results_diagnostic.csv`, with
missing-comparison reasons in `LDSC_Pair_Status.csv`. No estimate is filled in,
clipped or repaired. Failed batches can be retried with `--restart`; completed
matching batches remain reusable, including after changing the selection policy.

On success, use **both** `LDSC_Retained_Traits.csv` and `ldsc_results.csv`
downstream, so excluded traits are not requested again:

```bash
ldsc-gpca gpca \
  --input "${out_folder}/LDSC_Retained_Traits.csv" \
  --ldsc_results "${out_folder}/ldsc_results.csv" \
  --outdir "${out_folder}/gpca_validation" \
  --validate_only
```

This runs the usual strict matrix checks on the selected set. LDSC selection
does not bypass positive-definiteness or other genomicPCA/GWAMA requirements.
If fewer than two traits remain, selection reports are saved and the command
fails without publishing `ldsc_results.csv`.

A non-positive heritability or missing correlation is an estimation outcome,
not a reason to discard valid results from other pairs. Strict behavior remains
the default. Both `--result_failure_action error` and `report` preserve all
well-formed estimates and QC reports. `error` then exits unsuccessfully;
`report` finishes collection and returns the diagnostic table path when needed.
Neither policy imputes, clips or repairs estimates. Raw batch CSVs remain intact.
Malformed numeric text, missing columns/files and trait-name/coverage mismatches
still fail in report mode. Missing, malformed, incomplete or unchanged batch exports trigger the
same bounded retries as command failures; valid numerical failures are not retried.

Alternatively, to defer trait selection until genomicPCA, add
`--result_failure_action report` to your existing **`ldsc`** command, then run:

```bash
ldsc-gpca gpca \
  --input selected_traits.csv \
  --ldsc_results python_ldsc_results/ldsc_results_diagnostic.csv \
  --outdir gpca_validation \
  --failed_ldsc_action drop_traits \
  --validate_only
```

This uses the existing deterministic removal policy: failed self-pairs first,
then traits incident to the most remaining failed pairs, with ties resolved by
lower self-h2 Z and then later manifest position. It is a greedy complete-subset
heuristic, not a guarantee of the largest subset. Low positive h2/SE and finite
out-of-range rg are warnings rather than automatic significance filters.
Conflicting duplicate estimates are rejected **before** any trait can be removed.
The retained set must still have at least two traits, complete pair coverage,
positive-definite CTI and a finite positive PC1 eigenvalue. Genetic-matrix negative
eigenvalues remain reported under the existing warning/error policy; no nearest-PD
replacement is made. No correlation is replaced by zero.

`Python_LDSC_Retained_Traits.csv`, `Python_LDSC_Dropped_Failed_Traits.csv`,
`Python_LDSC_Failed_Pairs.csv` and `Python_LDSC_Self_Pair_QC.csv` audit the selection.
Selection reports are saved even if too few traits remain or later checks fail.
`GenomicPCA_Run_Status.csv` and `GWAMA_Run_Status.csv` distinguish validation,
completion and failure. Use the same selected policy across clusters; an external
multi-cluster runner can record each exit status and continue independent clusters.
This single-cluster command does not silently skip failed clusters.

Use a separate output directory per analysis. If compilation is repeated in an
existing directory, previous result/QC tables are preserved with `.previous-<id>`
suffixes so a failed rerun cannot advertise a stale `ldsc_results.csv` as current.
Result collection success never means that genomicPCA/GWAMA has passed validation.

All managed tabular inputs accept additional annotation columns and reordered
columns. Required canonical headers and scientific value checks still apply.
Old manifest aliases are allowed as extra annotations but cannot substitute for
missing canonical headers. GWAMA receives only the required nine fields in its
original order; its weighting implementation is unchanged.

## Run GenomicSEM LDSC

`ldsc-gpca ldsc --ldsc_backend genomicsem` runs GenomicSEM munging and LDSC from raw
tables or explicitly selected VCFs, or starts from existing munged files.
It produces the GenomicSEM sampling covariance matrices
needed by GenomicSEM. This command does not run Python LDSC, PCA or GWAMA.

The existing `ldsc-gpca genomicsem ldsc` entry point remains available with the
same workflow. In unified `ldsc`, select `--ldsc_backend genomicsem` explicitly;
its raw-table and VCF modes use GenomicSEM munging automatically. Reuse modes
skip munging. No separate `--munge_backend` override is accepted.

### GenomicSEM LDSC inputs: choose one mode

Input mode is selected explicitly, not detected from file contents or extensions.
Without an input-mode flag, the command expects raw paths in `sumstats_file`
and runs munging. `--vcf_input`, `--munged_dir` and `--munged_input` are mutually exclusive,
and one mode applies to every trait. `--input` is always the trait manifest CSV.
See the [input-mode explanations and examples](GENOMICSEM_LDSC.md#choose-the-input-mode-explicitly).

All four modes require at least two traits and these manifest columns:

| Column | Meaning |
| --- | --- |
| `traitname` | Unique, non-empty identifier; no whitespace or path separators. |
| `sample_prevalence` | Sample case fraction, or blank/`NA` for a quantitative trait. |
| `population_prevalence` | Population prevalence, or blank/`NA` for a quantitative trait. |

Supply **both** prevalence values strictly between 0 and 1 for binary liability
conversion, or leave both blank/`NA`. Mixed quantitative/binary rows are allowed;
a partially supplied prevalence pair fails. Unlike Python extraction, this route
does not derive case fractions or choose NC+NCO versus NEF for you.

**VCF mode — `--vcf_input`.** Include `vcf_files`. VCFs need one GWAS sample,
IDs and FORMAT `AF,ES,SE,LP,SI`; quantitative N defaults to NEF unless manifest
N overrides it. Binary VCFs require an explicit appropriate N and both
prevalences. The adapter retains INFO and allele frequency for GenomicSEM filters,
writes `vcf_input/` and `GenomicSEM_VCF_*` audits, then runs the unchanged
GenomicSEM workflow. See [commands, fields, sample-size rules and output columns](GENOMICSEM_LDSC.md#start-from-gwas-vcfs).

**Mode 1 — unmunged tables (default).** Also include `sumstats_file`:

```csv
traitname,sumstats_file,sample_prevalence,population_prevalence
Trait_A,/data/Trait_A.tsv,NA,NA
Trait_B,/data/Trait_B.tsv,NA,NA
Trait_C,/data/Trait_C.tsv,NA,NA
Trait_D,/data/Trait_D.tsv,0.2,0.05
```

The prevalence numbers illustrate the file format; use values justified for your
own trait/cohort. Each GWAS file is a whitespace-delimited table with a header:

| Variant column | Meaning / requirement |
| --- | --- |
| `SNP` | Variant ID matching the HapMap reference. |
| `A1`, `A2` | Effect and other allele. |
| `P` | Association P value. |
| `BETA` or `Z` | Signed effect relative to A1; supply a recognized signed statistic. |
| `N` | Sample size, unless a positive constant `N` is supplied in the manifest. |
| `INFO` | Optional recognized imputation-quality field for INFO filtering. |
| `MAF` or recognized allele-frequency field | Optional field used for MAF filtering. |

Optional positive manifest `N` supplies one constant per trait and replaces
raw-file N during GenomicSEM munging, even when that file already contains N. In
munged reuse mode, manifest N is ignored and existing munged N is used. GenomicSEM
munging interprets recognized column/sample-size conventions;
inspect its logs to confirm the interpretation and which filtering fields were
available. Do not assume every form of N is passed through unchanged.

**Mode 2 — munged directory.** Use `--munged_dir /data/munged`. Only the three
shared manifest columns above are required:

```csv
traitname,sample_prevalence,population_prevalence
Trait_A,NA,NA
Trait_B,NA,NA
Trait_C,NA,NA
Trait_D,NA,NA
```

The directory must contain exactly one `{traitname}.sumstats.gz` **or**
`{traitname}.sumstats` per trait. Having both for a trait is ambiguous and fails.

**Mode 3 — explicit munged paths.** Use `--munged_input` and a `munged_file` column:

```csv
traitname,munged_file,sample_prevalence,population_prevalence
Trait_A,/data/munged/Trait_A.sumstats.gz,NA,NA
Trait_B,/data/other/Trait_B.sumstats.gz,NA,NA
Trait_C,/data/munged/Trait_C.sumstats.gz,NA,NA
Trait_D,/data/munged/Trait_D.sumstats.gz,NA,NA
```

For both reuse modes, munged files must be **tab-separated** with
`SNP,A1,A2,N,Z` headers; column order may vary. Produce them with
`prepare --mode ldsc --munge_backend genomicsem` (or `--mode both`), using the
same `--hm3` allele reference for all traits. Do not directly reuse Python-munged
files: even with the same HapMap reference, strand-complement coding can remain
and be misaligned by GenomicSEM LDSC. `prepare --write_munge_inputs`
produces unmunged files and does not satisfy this contract directly.
Keep the entire preparation bundle. The tool verifies completed status, the
manifest/file/sidecar/reference checksums, pinned GenomicSEM preparation identity,
N source/convention/override, both prevalences, unique SNP IDs and exact reference
A1/A2 order. Missing provenance, swapped/complemented/incompatible pairs and
invalid N/Z fail with `GenomicSEM_Input_QC.csv`; no input is edited or silently
dropped. Older unverified files must be prepared again from raw inputs.
See the [GenomicSEM preparation and reuse example](GENOMICSEM_LDSC.md#reuse-munged-files).

All modes need `--ld_ref`, containing `<CHR>.l2.ldscore.gz` and
`<CHR>.l2.M_5_50` for chromosomes 1 through `--chromosomes`. Optional
`--ld_weights` contains the chromosome `.l2.ldscore.gz` weight files and defaults
to `--ld_ref`. Default munging additionally requires a whitespace-separated
`--hm3` reference with `SNP,A1,A2` headers. Omit `--hm3` when skipping munging.

### GenomicSEM LDSC commands

From unmunged tables:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/genomicsem_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/genomicsem_ldsc \
  --n_cores 4
```

From a munged directory:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/genomicsem_reuse_traits.csv \
  --munged_dir /data/munged \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/genomicsem_ldsc_reused
```

From explicit munged paths:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/genomicsem_munged_paths.csv \
  --munged_input \
  --ld_ref /references/eur_ld_chr \
  --ld_weights /references/eur_weights_chr \
  --outdir /results/genomicsem_ldsc_paths
```

### All GenomicSEM LDSC options

| Option | Default | Meaning / accepted values |
| --- | --- | --- |
| `--ldsc_backend CHOICE` | `python` on `ldsc` | Set `genomicsem` for these options. The existing `genomicsem ldsc` command implies GenomicSEM and does not take this selector. |
| `--input CSV` | Required | Manifest appropriate for one of the four modes above. |
| `--outdir DIRECTORY` | Required | Fresh or empty directory; existing non-empty directories are refused. |
| `--ld_ref DIRECTORY` | Required | LD-score and M reference files. |
| `--ld_weights DIRECTORY` | Same as `--ld_ref` | Separate regression weights. |
| `--hm3 FILE` | Unset | Required for raw-table or VCF munging; omit in either reuse mode. |
| `--munged_dir DIRECTORY` | Unset | Reuse per-trait files from `prepare --munge_backend genomicsem`; skip munging. Mutually exclusive with `--munged_input`. |
| `--munged_input` | Off | Reuse manifest `munged_file` paths from `prepare --munge_backend genomicsem`; skip munging. |
| `--vcf_input` | Off | Extract manifest `vcf_files`, then munge. Mutually exclusive with `--munged_dir` and `--munged_input`. |
| `--bcftools EXECUTABLE` | `bcftools` | VCF query executable; used only with `--vcf_input`. |
| `--p_min FLOAT` | `1e-300` | VCF P floor, finite and in `(0,1)`; adjusted records are audited. Can affect GenomicSEM Z reconstructed from P. |
| `--n_cores INTEGER` | `1` | VCF extraction and munging workers; ≥1. **Does not parallelize GenomicSEM LDSC regression.** |
| `--info_filter FLOAT` | `0.9` | GenomicSEM munging INFO threshold; `[0,1]`. Requires a recognized INFO field to filter on. |
| `--maf_filter FLOAT` | `0.01` | GenomicSEM munging MAF threshold; `[0,0.5]`. Depends on recognized frequency fields. |
| `--chromosomes INTEGER` | `22` | Read chromosomes 1 through this value; allowed 1–22. Does not change downstream split GWAMA's 22-file requirement. |
| `--n_blocks INTEGER` | `200` | Requested jackknife blocks; ≥2. GenomicSEM overrides the count for more than 18 traits; inspect its log. |
| `--chisq_max FLOAT` | Unset | Finite positive threshold. Unset uses GenomicSEM's automatic chi-square rule. |
| `--invalid_h2_action ACTION` | `drop` | `drop`: audit/remove non-positive or non-finite raw h2 traits before the second pass; `error`: stop. This default differs from downstream GPCA's strict default. |
| `--rscript EXECUTABLE` | `Rscript` | GenomicSEM Rscript program name/path. This option belongs to `genomicsem ldsc`; GPCA uses Rscript on PATH. |
| `--help` | Off | Show usage and exit. |

Munging filters are not reapplied in reuse modes. The wrapper runs an initial
`stand=FALSE` LDSC, audits h2 and applies `--invalid_h2_action`, then reruns
`stand=TRUE` on retained traits. Thus, it performs two GenomicSEM LDSC passes.
Fewer than two retained traits, failed munging or invalid final matrices stops
the run. No sampling covariance matrices are fabricated or repaired.

### GenomicSEM LDSC outputs

| Output under `--outdir` | Use |
| --- | --- |
| `genomicsem_LDSC.RData` | Final `LDSCoutput` with `S,V,I,S_Stand,V_Stand`; input to `genomicsem gpca`. |
| `genomicsem_LDSC_raw.RData` | Preliminary first-pass object; **do not use as GPCA input**. |
| `Selected_Traits.csv` | Retained traits plus GenomicSEM audit fields; extract `traitname` for the GPCA manifest as shown below. |
| `GenomicSEM_LDSC_Trait_QC.csv`, `GenomicSEM_LDSC_Events.csv` | h2 QC, removals and events. |
| `Resolved_Manifest.csv`, logs, `sessionInfo.txt` | Resolved inputs and runtime provenance. |
| `munge_output/` | Newly munged files from raw tables or VCFs. |
| `vcf_input/`, `GenomicSEM_VCF_QC_Summary.csv`, `GenomicSEM_VCF_QC_Issues.csv`, `GenomicSEM_VCF_Worker_Attempts.csv` | VCF mode only: extracted raw inputs, per-trait conversion counts/settings, original records affected by QC and worker attempts. [Schemas and interpretation](GENOMICSEM_LDSC.md#start-from-gwas-vcfs). |

The retained-trait audit `Selected_Traits.csv` can be used directly as a GPCA
manifest; extra internal columns are ignored. Alternatively, create a minimal
manifest with this command, preserving exact trait strings and row order:

```bash
python - /results/genomicsem_ldsc/Selected_Traits.csv /results/genomicsem_ldsc/Selected_Traits_for_GPCA.csv <<'PY'
import csv
import sys

with open(sys.argv[1], newline='', encoding='utf-8-sig') as source:
    traits = [row['traitname'] for row in csv.DictReader(source)]
with open(sys.argv[2], 'x', newline='', encoding='utf-8') as target:
    writer = csv.writer(target)
    writer.writerow(['traitname'])
    writer.writerows([trait] for trait in traits)
PY
```

The destination must not already exist. If using automatic VCF preparation
later, add matching `vcf_files` paths to this clean manifest. Otherwise, supply
`--gpca_input_folder`.

## Run GPCA and GWAMA

Use `ldsc-gpca gpca` for a Python pairwise table or
`ldsc-gpca genomicsem gpca` for GenomicSEM RData. Each command reads
existing LDSC estimates, computes PCA and optionally runs PC1 GWAMA; neither
reruns LDSC. The common options and outputs below apply to both.

### Common inputs

Save `/data/selected_traits.csv` as CSV:

```csv
traitname
Trait_A
Trait_B
Trait_C
Trait_D
```

Each row selects one trait; order is authoritative. Extra traits in the LDSC
results are allowed. For GWAMA, supply a folder of per-trait TSVs or use
[automatic preparation](#automatic-vcf-preparation).

GWAMA inputs must contain the **nine required columns below**, separated by
actual tabs. Extra columns and any input column order are allowed; the reader
selects and reorders the required fields before calling GWAMA. This example contains tabs:

```tsv
SNPID	CHR	BP	EA	OA	EAF	N	Z	P
1_1000_A_G	1	1000	G	A	0.25	33000	2	0.0455002638963584
```

| Column | Meaning |
| --- | --- |
| `SNPID` | Consistent variant identifier across traits. |
| `CHR` | Chromosome. Split mode uses autosomes 1–22. |
| `BP` | Base-pair position in the common genome build. |
| `EA`, `OA` | Effect and other allele; Z must refer to EA. |
| `EAF` | Frequency of EA, between 0 and 1. |
| `N` | Positive per-variant sample-size value used by GWAMA. Verify its convention. |
| `Z` | Signed association Z statistic for EA. |
| `P` | Association P value. |

Existing variant-header handling still accepts `A1`, `A2`, `p` for `EA`, `OA`,
`P` in memory when their canonical counterparts are absent. When both are
present, the canonical field is used. The canonical output headers above are
unchanged. Missing or duplicate required fields fail. The reader checks file existence and schema;
**it does not provide complete per-variant numerical QC**. Validate externally
prepared files, sample sizes, alleles and genome builds before use.

| Layout | Required filename for every retained trait |
| --- | --- |
| `--splitby_chr split` (default) | `{traitname}_chr1_GenomicPCA_inputs.tsv` through `{traitname}_chr22_GenomicPCA_inputs.tsv` |
| `--splitby_chr nosplit` | `{traitname}_GenomicPCA_inputs.tsv` |

Choose the same layout during preparation and GPCA. A single munged
`.sumstats.gz` cannot substitute for these files because it lacks required
position, frequency and P columns.

### Python LDSC results table

For `ldsc-gpca gpca --ldsc_results`, use a headered CSV, TSV or whitespace-delimited
text table; `.gz` compression is accepted and the delimiter is detected from the
header. Supply the following columns (order need not match this table):

| Column(s) | Meaning / requirement |
| --- | --- |
| `p1`, `p2` | Trait identifiers matching the selected manifest; managed compiled results already use trait names. |
| `rg` | Genetic correlation estimate. Self-pair rg must be approximately 1. |
| `se` | Positive SE of rg. |
| `z` | Reported rg Z statistic; checked against rg/SE. |
| `p` | Reported rg P value in `[0,1]`; diagnostic, not a selection threshold. |
| `h2_int`, `h2_int_se` | Univariate LDSC intercept and positive SE for target trait p2. Self-pairs provide the GWAMA CTI diagonal. |
| `gcov_int`, `gcov_int_se` | Cross-trait LDSC intercept and positive SE; off-diagonal CTI entries. |
| `h2_obs`, `h2_obs_se` | Observed/unconverted h2 estimate and SE, when this scale is supplied. |
| `h2_liab`, `h2_liab_se` | Liability h2 estimate and SE, when this scale is supplied. |

The first seven rows of the table are always required. In addition, supply at
least one complete h2/SE column pair. A trait's selected scale must have usable
values; self-pair h2 must be positive. h2 fields describe **p2**, and self-pairs
are authoritative for each trait's h2 QC. Do not manually relabel effective-N
binary estimates as liability or ordinary observed-scale estimates.

For `k` selected traits, include every unordered pair and every self-pair:
**`k × (k + 1) / 2` unique pairs**. For three traits this is A–A, A–B, A–C,
B–B, B–C and C–C; for 180 traits it is **16,290**. Upper/lower triangles and
full symmetric tables are accepted. Consistent duplicate orientations collapse;
conflicting duplicates fail. A raw `ldsc.py --rg A,B,C` command does not supply
B–C, so it is insufficient by itself.

Every unordered pair is checked for agreement in `rg,se,z,p,gcov_int,gcov_int_se`.
Repeated self-pairs additionally require agreement in each supplied
`h2_obs,h2_obs_se,h2_liab,h2_liab_se,h2_int,h2_int_se` field. The scales are
compared separately, including a populated alternative scale; an entirely empty
alternative scale is allowed. Off-diagonal reverse-row heritability fields refer
to different target traits and are not required to agree. Conflicts are fatal
under every result/trait-removal policy, before averaging or excluding traits.
Each source self estimate must be valid: a positive average cannot rescue a
non-positive source heritability. Identical duplicates remain accepted.

Required selected numeric values must be finite, required SEs positive and P
in range, with one self-correlation exception: a self-pair with finite `rg`
within the existing tolerance of 1 may have `se=0`, `z=+Inf`, and `p=0` together.
All other checks still apply, including positive heritability/intercept SEs.
Zero SE for a between-trait correlation, negative SE, missing values, or a
different zero-SE/Z/P combination remains invalid. Python result collection
uses tolerance `0.01`; R uses `--self_rg_tolerance` (default `0.01`).

The exception preserves the native values and is recorded in
`LDSC_Pair_Status.csv` (`Warning`) and `Python_LDSC_Self_Pair_QC.csv`
(`Self_RG_Zero_SE_Exception`, `Self_QC_Note`). The existing audit fields
`Required_Numeric_Finite` and `Required_SE_Positive` retain their literal meaning;
`Required_Numeric_Valid` and `Required_SE_Valid` include the exception. The PCA
diagonal remains 1 and its calculations do not use self-correlation SE/Z/P.
No CLI opt-in or original Python LDSC modification is needed.

Missing/failed estimates stop by default. No values are filled in or
selected using a significance threshold. Managed LDSC outputs may include both
h2 scale pairs with unused cells empty and additional provenance columns;
preserve these rather than deleting/relabeling them to force a scale.

### Optional trait-wide correlation normalization

Managed Python LDSC exports three additional values directly from each fitted
pair, before formatting the readable log. Compilation retains them and adds the
normalization columns in memory using Polars, without rereading the final CSV:

| New column | Meaning |
| --- | --- |
| `gcov_native_obs` | Native genetic covariance: `RG.gencov.tot`, on the unconverted observed scale. |
| `h2_p1_pair_obs` | First trait's pair-specific heritability: `RG.hsq1.tot`, on the unconverted observed scale. |
| `h2_p2_pair_obs` | Second trait's pair-specific heritability: `RG.hsq2.tot`, on the unconverted observed scale. |
| `gcov_pair` | Covariance reconstructed from original rg and the pair-specific h2 values, on their reported scales. Retained for compatibility. |
| `rg_trait_wide` | Native covariance divided by the square root of the two **self-pair** observed h2 estimates; self-correlations are set to 1. Older exports use reconstructed covariance. |
| `normalization_status` | `calculated` or the reason the new values are unavailable. |

For new native exports, each pair already contains both heritabilities:

```text
rg(A,B)           = gcov_native_obs(A,B) / sqrt(h2_p1_pair_obs × h2_p2_pair_obs)
rg_trait_wide(A,B) = gcov_native_obs(A,B) / sqrt(h2_A_self_obs × h2_B_self_obs)
```

The reconstruction uses the [ratio definition in the pinned Python LDSC source](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py#L702-L713).

The native fields are exported without liability conversion, even when the
original `h2_liab` column is liability-scale. If the analysis uses effective
sample size (NEF), “observed” here means LDSC's unconverted output; it does not
assert that NEF-based h2 has the usual population observed-scale interpretation.
Self denominators come from `h2_p2_pair_obs` on A–A and B–B rows. A triangular
table with these native fields and all self-pairs supports trait-wide normalization.
Native covariance and pair h2 must reproduce original rg, and reported/native
h2 scale factors must agree with the corresponding self-pair.

Older files without native fields retain the existing fallback: the A–B row
supplies B's pair h2 and the B–A row supplies A's; `gcov_pair = rg ×
sqrt(h2_A_pair × h2_B_pair)`, then divide by `sqrt(h2_A_self × h2_B_self)`.
This fallback requires full precision results, both orientations and self-pairs
from the same analysis/settings. It does not manufacture native export columns.

Neither export nor normalization changes original `rg`, `se`, `z`, `p`, h2,
intercepts or readable logs, and no extra regressions are added to a fresh run.
Existing restart checkpoints are invalidated by the changed exporter fingerprint;
rerunning LDSC obtains native fields for old batches. The standalone annotation
command below can instead reuse older CSVs through the reconstruction fallback.
Missing, non-positive, non-finite or conflicting h2, and ambiguous/inconsistent
scales produce empty derived values
and an explanatory status. Each trait's pair and self h2 must use the same scale;
different traits may use different scales without conversion. The reconstructed
covariance inherits those units. Values outside [-1,1] are not clamped.

**PCA/GWAMA continues to use original Python `rg` by default**
(`--rg_normalization pair`). Select the additional values explicitly:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /data/genomicpca_inputs \
  --splitby_chr nosplit \
  --outdir /results/python_gpca_trait_wide \
  --rg_normalization trait_wide \
  --validate_only
```

This example inspects PCA without GWAMA. To run GWAMA, omit `--validate_only`
and supply the [required export metadata](#all-automatic-export-options). Trait-wide mode requires
finite `rg_trait_wide` values and `normalization_status=calculated` for every
selected row. Existing coverage, duplicate and original-estimate QC still apply.
It changes the PCA matrix and potentially PC1/GWAMA results; CTI stays unchanged.
The matrix output contains the selected values; PC1/GWAMA and validation audits
record `RG_Normalization`.
`Global_Genetic_Correlations.csv` retains original rg/SE/Z/P and adds `rg_used`.
**Original SE/Z/P describe original rg only**; no new SE, P value or `V_Stand`
is inferred for the normalized estimate. This option addresses normalization,
not LDSC regression-weight differences, and does not guarantee GenomicSEM equality
or greater statistical accuracy. With covariance PCA, the selected rg is multiplied
by selected-scale self-pair h2. With observed self h2 and native-export
normalization, this recovers the native observed covariance numerator; with
older reconstructed normalization, it recovers that covariance on the matching
scales. Other selected h2 scales rescale the covariance.

To append the same columns to an older saved CSV without rerunning LDSC:

```bash
python -m ldsc_gpca.normalization \
  --input /results/python_ldsc/ldsc_results.csv \
  --out /results/python_ldsc/ldsc_results_with_normalization.csv
```

Input can be gzip-compressed; output is plain CSV. Original cell text and row order
are preserved. Reusing the input path performs an atomic replacement; a separate
output is useful for comparisons. Repeated annotation recalculates the three columns.
See [normalization validation](../tests/NORMALIZATION_VALIDATION.md) and
[native-export validation](../tests/NATIVE_EXPORT_VALIDATION.md) for real-data results.

### GenomicSEM results object

For `ldsc-gpca genomicsem gpca --ldsc_results`, use a binary `.RData` containing
an object named `LDSCoutput` with `S,V,I,S_Stand,V_Stand`, appropriate trait names
and matching dimensions. The managed `genomicsem ldsc` command produces this as
`genomicsem_LDSC.RData`. For external GenomicSEM results, obtain a standardized
`stand=TRUE` output and save the correctly named object. A renamed CSV or the
first-pass `genomicsem_LDSC_raw.RData` is not suitable.

### Start with QC/PCA only

Python results:

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
  --input /results/genomicsem_ldsc/Selected_Traits_for_GPCA.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --outdir /results/genomicsem_gpca_qc \
  --validate_only
```

Review PC1 weights, eigenvalues, CTI and trait-removal/warning audits.
`--validate_only` skips VCF preparation, GWAMA and final export. It **does not
open or validate the per-variant GWAMA files**.

### Run SNP-level PC1 GWAMA after reviewing QC

The bundled GWAMA function does not output INFO. Its automatic summary export
therefore needs a scientifically justified `--gwama_output_info` value. There
is no universally correct constant: this is supplied metadata, not estimated
imputation quality. The examples use a shell variable `GWAMA_INFO`; assign it
your justified value first. The guard stops before analysis if it is unset.
If no constant is scientifically defensible, resolve the downstream INFO policy
before using this exporter; the current command otherwise fails at export after
GWAMA. A custom source that provides INFO can omit the override.

Python LDSC backend:

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/python_gpca_gwama \
  --splitby_chr split \
  --n_cores 22 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

GenomicSEM backend:

```bash
ldsc-gpca genomicsem gpca \
  --input /results/genomicsem_ldsc/Selected_Traits_for_GPCA.csv \
  --ldsc_results /results/genomicsem_ldsc/genomicsem_LDSC.RData \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/genomicsem_gpca_gwama \
  --splitby_chr split \
  --n_cores 22 \
  --dataset_id cluster1 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

To use whole-genome files, change `--splitby_chr split` to `nosplit`. GWAMA then
runs as one job; increasing `--n_cores` does not split that job into chromosomes.
To explicitly permit removal of failed LDSC traits, add
`--failed_ldsc_action drop_traits` and review its audit. Dropping traits changes
the component being analyzed; it is not merely a speed setting.

### All shared GPCA analysis options

These options apply to **both** GPCA commands. Preparation and export options
follow in separate tables; Python-table-specific options are listed afterward.

| Option | Default | Meaning / accepted values |
| --- | --- | --- |
| `--input CSV` | Required | Selected-trait manifest; `traitname` required. |
| `--ldsc_results FILE` | Required | Python table for `gpca`; GenomicSEM RData for `genomicsem gpca`. |
| `--outdir DIRECTORY` | Required | QC, PCA, GWAMA and export destination. Use a fresh directory. |
| `--gpca_input_folder DIRECTORY` | Unset | Existing nine-column per-trait tables. Omit for automatic VCF preparation, or for `--validate_only`. |
| `--source_path FILE` | Bundled `N_weighted_GWAMA.function.1_2_6.R` | Optional trusted R source defining the already-modified `my_GWAMA` or `multivariate_GWAMA` function. Normally unnecessary. |
| `--splitby_chr MODE` | `split` | `split` for chromosomes 1–22; `nosplit` for one whole-genome file per trait. |
| `--n_cores INTEGER` | `0` | Chromosome GWAMA workers; `0` selects automatically, positive values request workers. Also controls final export concurrency; [details below](#cpu-workers-and-numerical-threads). |
| `--validate_only` | Off | Run LDSC QC/PCA and audits; skip preparation, GWAMA and export. |
| `--allow_missing_traits` | Off | Remove manifest traits absent from LDSC input with an audit; otherwise stop. Distinct from missing/failed estimates among present traits. |
| `--failed_ldsc_action ACTION` | `error` | `error` stops on invalid selected traits/pairs; `drop_traits` applies audited deterministic trait removal to seek a usable subset. |
| `--h2_z_warn_threshold FLOAT` | `2` | Warn for retained h2/SE below this value; finite ≥0. `0` disables. Does not remove traits. |
| `--rg_out_of_range_action ACTION` | `warn` | `warn` or `error` for finite off-diagonal rg outside `[-1,1]`; never clamp. |
| `--pca_matrix MODE` | `correlation` | `correlation` uses standardized genetic correlations; `covariance` uses GenomicSEM S or reconstructed covariance. Scale interpretation differs. |
| `--pc1_orientation MODE` | `tutorial` | `tutorial`: reverse the whole PC1 loading vector if its median is negative; `as_computed`: retain the eigenvector sign. |
| `--negative_eigen_action ACTION` | `warn` | `warn` or `error` for substantive negative PCA eigenvalues. Raw eigenvalues remain reported. |
| `--matrix_eigen_tolerance FLOAT` | `1e-8` | Finite positive relative cutoff for substantive negative eigenvalues. |
| `--help` | Off | Show grouped analysis, automatic preparation and export options without running an analysis. |
| `--prepare_help` | Off | Show automatic preparation options and file requirements. |
| `--postprocess_help` | Off | Show automatic export options. |

### Additional options for Python-table GPCA only

These are accepted by `ldsc-gpca gpca`, **not** `ldsc-gpca genomicsem gpca`.

| Option | Default | Meaning / accepted values |
| --- | --- | --- |
| `--rg_normalization MODE` | `pair` | `pair` uses original `rg`; `trait_wide` uses precomputed `rg_trait_wide`. Changes PCA and downstream PC1/GWAMA, not CTI. |
| `--heritability_scale MODE` | `auto` | `auto`, `observed`, `liability`, or `mixed`; policy described below. Selects h2 columns, not a conversion procedure. |
| `--duplicate_tolerance FLOAT` | `0.001` | Finite positive maximum absolute difference between duplicate pair estimates, except z. |
| `--duplicate_z_tolerance FLOAT` | `0.01` | Finite positive maximum absolute duplicate-z difference. |
| `--self_rg_tolerance FLOAT` | `0.01` | Finite positive allowed absolute difference between self-rg and 1. |
| `--comparison_epsilon FLOAT` | `1e-12` | Finite non-negative allowance in boundary comparisons. |
| `--z_consistency_tolerance FLOAT` | `0.01` | Finite positive relative tolerance comparing reported z with rg/SE. |
| `--z_consistency_action ACTION` | `warn` | `warn` or `error` when the z consistency check fails. |
| `--ldsc_chunk_size INTEGER` | `250000` | Positive number of LDSC rows per read chunk; tune reading/memory, not statistical results. |
| `--color MODE` | `auto` | Help colour: `auto`, `always`, `never`. Auto respects `NO_COLOR`. |

Scale policy:

| Choice | Behavior |
| --- | --- |
| Correlation PCA + `auto` | Select each trait's only populated h2/SE scale; ignore entirely empty column pairs. Different traits can use different scales for h2 QC. Ambiguous populated scales fail. This scale option does not rescale rg or CTI; correlation normalization is controlled separately. |
| Covariance PCA + `auto` | Require one complete common self-h2/SE scale. Two complete scales or no complete common scale fails. |
| `observed` or `liability` | Use that scale for all selected traits; no fallback or conversion. |
| Covariance PCA + `mixed` | Explicitly permit trait-specific scales when each is unambiguous; warns about scale dependence. This does not establish common scientific units. |

### Automatic VCF preparation

Omit `--gpca_input_folder` and put `traitname,vcf_files` in the selected manifest
to run preparation before GPCA. The VCF schema and filtering are exactly those
of [the preparation module](#prepare-per-trait-gpca-inputs). Example:

```bash
ldsc-gpca gpca \
  --input /data/prepare_traits.csv \
  --ldsc_results /results/python_ldsc/ldsc_results.csv \
  --outdir /results/gpca_from_vcf \
  --prepare_workers 4 \
  --n_cores 22 \
  --gwama_output_info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

Prepared files go to `<outdir>/gpca_inputs/`. All original manifest traits must
prepare successfully before R starts, including traits that might later fail
LDSC QC. `--validate_only` skips preparation. With an existing input folder,
non-default preparation settings are rejected; prepare new files separately.
After an R failure, reuse successfully prepared files by explicitly supplying
their folder on the next run.

| Automatic preparation option | Default | Meaning |
| --- | --- | --- |
| `--prepare_workers INTEGER` | `4` | Trait-preparation workers, ≥1. Standalone `prepare` calls this `--n_cores`; here it is separate from GWAMA/export workers. |
| `--p_min FLOAT` | `1e-300` | Finite P floor strictly between 0 and 1. |
| `--gpca_id_source MODE` | `chr_pos_ref_alt` | `chr_pos_ref_alt` or `vcf_id`. |
| `--write_munge_inputs` | Off | Also export unmunged HapMap-selected tables. |
| `--hm3 FILE` | Unset | Tab-separated `SNP` list; required with `--write_munge_inputs`. |
| `--munge_id_source MODE` | `vcf_id` | `vcf_id` or `chr_pos_ref_alt`; must match HapMap IDs. |
| `--bcftools EXECUTABLE` | `bcftools` | Local executable name/path. |

### All automatic export options

Both GPCA commands run export after successful GWAMA. No separate command is
needed. `--validate_only` skips it.

| Option | Default | Meaning / accepted values |
| --- | --- | --- |
| `--dataset_id STRING` | Basename of `--outdir` | Safe filename prefix for final combined/summary files and export JSON. |
| `--gwama_output_n_eff FLOAT` | Preserve source N_eff | Optional finite positive constant in the selected-column summary only. Can supply N_eff if absent from the source. |
| `--gwama_output_info FLOAT` | Preserve source INFO | Optional finite constant in `[0,1]` in the selected-column summary only. Needed for the bundled GWAMA source because it omits INFO. |
| `--gzip_level INTEGER` | `1` | Final-export compression level, 1–9. Level 1 favors speed; 9 generally gives smaller files at greater CPU cost. Contents are unchanged. |
| `--archive_chromosomes` | Off | After saving exports, move current-run source results/logs into `chromosome_wise/`. Otherwise keep them in place. |

These constants do not alter input N, SNP filtering, LDSC, PCA, GWAMA weights or
the combined result table. They are **user-supplied export metadata**, not new
estimates. For example, `--gwama_output_n_eff 33000` would write 33000 into the
summary's N_eff column for every SNP; use it only when justified. Without an
override, the corresponding source column must exist. Retain the JSON audit.

## INFO filtering, sample size and export metadata

These settings act at different stages. An input threshold can change the SNPs
used for LDSC, and an analysis sample size can change estimates or GWAMA weights.
An export override only writes a constant into the final selected-column summary.

| Setting or input | Default / source | Where it acts |
| --- | --- | --- |
| Python `--info_min` | `0.7`; VCF FORMAT/SI | Inclusive extraction threshold, before Python LDSC munging. Accepted values: `[0,1]`. |
| GenomicSEM `--info_filter` | `0.9`; recognized INFO field | GenomicSEM munging threshold. Accepted values: `[0,1]`; an absent field cannot be filtered. |
| Python LDSC sample size | NEF without population prevalence; NC+NCO with it | N supplied during managed Python munging. Binary conversion also needs sample/population prevalence. |
| GenomicSEM sample size | Recognized source N; optional manifest `N` for a per-trait constant | GenomicSEM raw-table munging. Inspect logs for recognized fields and N conventions. |
| GWAMA sample size | Each input TSV's `N`; VCF preparation copies NEF | SNP-specific GWAMA weights. The LDSC manifest N does not replace GWAMA N. |
| `--gwama_output_info` | Preserve source INFO; bundled GWAMA omits it | Optional finite `[0,1]` constant in the summary only; required when the source has no INFO. |
| `--gwama_output_n_eff` | Preserve source N_eff | Optional finite positive constant in the summary only. It does not estimate a new effective sample size. |

There is no general managed `--sample_size` option. Choose the correct input N
and prevalence convention before analysis. In munged-file reuse mode, extraction
and munging INFO/MAF filters are not reapplied; GenomicSEM manifest N does not replace
N already in a munged file. Preparation does not itself apply an INFO threshold.

For example, **if scientifically justified for a particular downstream summary**:

- `--gwama_output_info 0.9` writes INFO=0.9 for every exported summary row. It
  does not assert that every variant passed a 0.9 imputation-quality threshold.
- `--gwama_output_n_eff 33000` writes N_eff=33000 for every summary row. It
  does not change GWAMA input N, weights, BETA, SE, Z or P.

Neither override changes the combined GWAMA table. Omitting an override preserves
the corresponding source column, which must exist. No universally appropriate
INFO or N_eff constant is supplied. The workflow examples use shell variables
`GWAMA_INFO` and `GWAMA_N_EFF` with guards so values must be chosen before running;
remove the N_eff option when no override is intended. A trusted custom GWAMA
source supplying INFO may omit the INFO override.

See the full examples for [Python GPCA/GWAMA](PYTHON_GPCA.md#exports-and-what-to-inspect),
[GenomicSEM GPCA/GWAMA](GENOMICSEM_GPCA.md#results-and-export-settings),
and [pipeline mode](PIPELINE.md#change-filters-pca-and-export-settings).

## Understand QC and scale choices

| Quantity | Python-table backend | GenomicSEM backend |
| --- | --- | --- |
| Correlation PCA matrix | Diagonal 1; off-diagonal original rg by default, or explicitly selected `rg_trait_wide` | `S_Stand` |
| Covariance PCA matrix | Self h2 on diagonal; off-diagonal `selected rg × sqrt(h2_i × h2_j)` | GenomicSEM `S`, on its supplied scales |
| GWAMA CTI/error covariance | Self `h2_int` on diagonal; pairwise `gcov_int` off-diagonal | GenomicSEM `I` |

CTI must be symmetric and positive definite; it is not the genetic correlation
matrix. PCA uses symmetric eigendecomposition. PC1 loadings equal the PC1
eigenvector multiplied by the square root of its positive eigenvalue. Required
trait order is preserved across all matrices, loadings and GWAMA input lists.
Pairwise P, z and SE are diagnostics, not PCA significance filters or weights.

Negative eigenvalues and rg outside `[-1,1]` are reported. The square-root guard
uses `max(eigenvalue,0)` for loadings but does not repair the underlying matrix.
PC1 must be finite/positive. The package never silently replaces CTI or the
genetic matrix with a nearest-positive-definite matrix. A PC1 PASS is not an
overall scientific-validity or run-success guarantee.

Opt-in trait removal is deterministic but is not guaranteed to find the largest
valid subset. Missing traits, failed estimates and structurally invalid matrices
are distinct problems; removal options do not solve every error. Review removal
reports because the retained set defines the phenotype. Covariance PCA depends
on h2 scale/units and needs additional care with mixed binary/quantitative traits.

The bundled modified GWAMA function receives PC1 loadings in its argument named
`h2`; these are **not SNP heritabilities**. A custom source must retain the
tutorial's weighting:

```r
W <- t(t(sqrt(N)) * h2)
sqrt_W <- W
```

The Python pairwise table does not contain GenomicSEM's full cross-estimate
sampling covariance matrices `V` and `V_Stand`. The package does not invent
diagonal substitutes. Neither GPCA command runs `paLDSC`; the Python input cannot
support procedures that require those missing matrices. See the
[default GPCA checks](#default-gpca-checks) for
backend-specific tolerances and removal rules.

### Default GPCA checks

| Condition / setting | Default action |
| --- | --- |
| Selected trait absent | Stop; `--allow_missing_traits` enables audited removal |
| Missing/invalid selected estimates | Stop; `--failed_ldsc_action drop_traits` enables a deterministic removal heuristic |
| Positive self h2 and valid required SEs | Required; scale selection never rescues invalid values by switching scales |
| Weak h2/SE | Warn below `--h2_z_warn_threshold 2`; no removal; `0` disables |
| Finite off-diagonal abs(rg)>1 | `--rg_out_of_range_action warn`; alternative `error`; never clamp |
| Substantive negative PCA eigenvalues | `--negative_eigen_action warn`; alternative `error`; relative tolerance `1e-8` |
| Asymmetric/structurally invalid matrices or non-positive-definite CTI | Stop; no nearest-positive-definite repair |
| PC1 eigenvalue | Must be finite and positive |

Python-specific defaults are self-rg tolerance **0.01**, duplicate absolute tolerance
**0.001**, duplicate-z tolerance **0.01**, and boundary epsilon **1e-12**.
Reported z versus rg/SE inconsistency warns by default at relative tolerance 0.01.
The Python table QC requires finite selected numerical values, positive required
SEs and result p in `[0,1]`, with one self-pair exception: finite rg
within tolerance of 1, SE=0, Z=+Inf and P=0 is accepted and audited
when all other checks pass; this **result-table** rule differs from munging's input
P requirement `0 < P <= 1`. GenomicSEM standardized diagonal checks use absolute
tolerance `1e-8`, not the Python self-pair tolerance.

Removal is not guaranteed to retain the largest possible subset. Python removes
failed self traits first, then traits with the most remaining failed pairs; ties
use lower self h2/SE and then later manifest position. GenomicSEM pair-removal ties
use manifest order. Removing traits changes the matrix and can change PC1/GWAMA;
inspect audits and compare sensitivity analyses rather than assuming negligible impact.
Structural errors and ambiguous scale selection are not solved by enabling removal.

Negative eigenvalues remain reported. The loading calculation uses
`pmax(eigenvalue, 0)` to avoid square roots of negative values; **the matrix itself
is not repaired**. A PC1 PASS is not a guarantee that the matrix or entire run passed.

## Find and interpret the outputs

### QC and PCA reports

Both GPCA backends write these CSV reports under `--outdir` once validation
reaches the corresponding stage, including in validation-only runs:

| File | What it tells you |
| --- | --- |
| `GenomicPCA_CTI_Used.csv` | Exact LDSC intercept/error-covariance matrix used. |
| `GenomicPCA_Correlation_Matrix_Used.csv` | Selected genetic-correlation matrix. |
| `GenomicPCA_PCA_Matrix_Used.csv` | Exact matrix decomposed; check this when choosing covariance PCA. |
| `GenomicPCA_Covariance_Matrix_Used.csv` | Covariance matrix when available; Python writes it for covariance mode, GenomicSEM reporting includes S. |
| `GenomicPCA_PC1_Weights_Used.csv` | Signed PC1 loadings and orientation used for GWAMA. |
| `GenomicPCA_Selected_Traits_Eigenvalues.csv` | Raw eigenvalues and loading-guard diagnostics. |
| `GenomicPCA_All_PCs_Variance_Explained.csv` | Every PC's eigenvalue, raw/cumulative percentages, positive-normalized percentages and negative-eigenvalue flags. |
| `GenomicPCA_PC1_QC.csv` | PC1 validity/direction and reasons; not overall run status. |
| `GenomicPCA_PC1_Protein_Contributions.csv` | Retained traits in manifest order, signed coefficients/loadings, contribution percentages and ranks. |
| `Global_Heritability_Estimates.csv`, `Global_Genetic_Correlations.csv` | Heritability/correlation summaries; GenomicSEM reports cover the retained selected set. |
| `GWAMA_Run_Status.csv` | Per-chromosome or whole-genome GWAMA success/failure and output paths when GWAMA runs. |

Raw explained percentage is `100 × eigenvalue / sum(all eigenvalues)`. Negative
eigenvalues make this unlike a conventional non-negative variance partition;
percentages can be negative or exceed 100%. Positive-normalized percentages
describe the positive part and do not repair the matrix. These are summaries of
the estimated genetic matrix, not measured phenotypic variance.

The PC1 trait contribution is `100 × (PC1 eigenvector coefficient)²`. It sums to
100% for valid PC1 and is invariant to whole-vector sign reversal. It is not a
causal contribution or the final SNP-specific GWAMA weight, which also depends
on sample size. Interpret effect signs relative to the recorded PC orientation.

Python-specific audits include:

- `Python_LDSC_Input_Validation_Summary.csv`
- `Python_LDSC_Heritability_Scales.csv`, `Python_LDSC_Self_Pair_QC.csv`
- `Python_LDSC_RG_SE_Matrix.csv`, `Python_LDSC_Intercept_SE_Matrix.csv`
- `Python_LDSC_Dropped_Missing_Traits.csv`, `Python_LDSC_Dropped_Failed_Traits.csv`

GenomicSEM audits include `GenomicSEM_QC_Events.csv`, `GenomicSEM_QC_Warnings.csv`,
`GenomicSEM_QC_Removed_Traits.csv`, `GenomicSEM_Trait_QC_Summary.csv`,
`GenomicSEM_Retained_Manifest.csv`, `GenomicSEM_LDSC_Used.RData` and settings/session
files. Audit availability depends on how far the run reached; early input/scale
errors can prevent later reports. Always check command exit status and errors.

### GWAMA results and final exports

Only current-run successful files listed in the updated `GWAMA_Run_Status.csv`
are exported. Bundled correlation-GWAMA chromosome files use
`Chr{CHR}_GenomicPCA_PC1.N_weighted_GWAMA.results.txt.gz`; covariance mode uses
`GenomicPCA_Covariance_PC1` in the name. Whole-genome mode omits the chromosome
prefix. A custom source must follow the expected result-file convention.

| Final output | Location and contents |
| --- | --- |
| `{dataset_id}_GWAMA_combined_results.txt.gz` | Under `--outdir`: gzip-compressed TSV with all source columns plus `count_question,count_plus,count_minus` derived from `Direction`; sorted by chromosome/position. |
| `{dataset_id}_GPCA_inputs.txt.gz` | Under `<outdir>/harmonisation_input/`: gzip-compressed selected-column TSV described below. |
| `{dataset_id}_postprocess.json` | Under `--outdir`: sources, row count, overrides, engine, actual workers, compression and timings. |
| `chromosome_wise/` | Source results/logs moved here only if `--archive_chromosomes` is enabled. |

The selected-column export has these **12 tab-separated columns**, in order:

```tsv
SNPID	CHR	BP	EA	OA	EAF	N_eff	BETA	SE	Z	PVAL	INFO
```

**This is a downstream summary, not the nine-column per-trait GWAMA input.**
Its `_GPCA_inputs` filename does not make it interchangeable with those inputs.
`N_eff,BETA,SE,Z,PVAL,INFO` are the corresponding GWAMA-reported fields, except
for explicit N_eff/INFO overrides. Export itself does not harmonise alleles,
change genome build or recalculate statistics.

Export validates required columns, unique/non-empty SNPIDs, positions, Direction,
current-run files and override ranges. It does not comprehensively QC every
BETA, SE, PVAL or preserved N_eff/INFO value. Existing export/archive paths are
refused. If archiving fails, saved exports and some moved originals can remain;
inspect the error and files before retrying.

## CPU workers and numerical threads

### What `--n_cores` controls

| Command/stage | Default | What is parallel |
| --- | --- | --- |
| `prepare` | `4` | Preparation of separate traits; also the default shared Polars pool at CLI startup. Explicit `POLARS_MAX_THREADS` is preserved. |
| `ldsc` | `5` | Separate extraction/munging tasks and pairwise LDSC batches. |
| `genomicsem ldsc` | `1` | Munging only. Its two GenomicSEM LDSC passes remain sequential. |
| Either GPCA, split GWAMA | `0` = auto | Up to 22 chromosome jobs; auto uses `min(22, max(1, detected physical cores - 1))`, falling back to 1 if detection fails. |
| Either GPCA, whole-genome GWAMA | `0` = auto | One GWAMA job. Windows GWAMA is sequential even in split mode. |
| Either GPCA, automatic preparation | Separate `--prepare_workers 4` | Preparation workers, independent of GWAMA/export concurrency. |
| Either GPCA, final export | Follows `--n_cores` | Concurrent chromosome reads, shared Polars pool and pigz compression workers. Auto uses up to 22 available logical CPUs, respecting CPU affinity where supported. |

For a machine with 100 CPUs and 900 GB RAM, `--n_cores 22` is a reasonable
**starting point** for one chromosome-split GPCA run, because there are only 22
chromosome jobs. It is not a measured optimum for your data/storage. Python LDSC
can have more runnable batches; benchmark progressively larger worker counts
and watch CPU, memory and disk throughput. Use only CPUs allocated to your job
on a scheduler, even if the machine has more.

When running several clusters at once, budget workers across them. Each GWAMA
worker reads many trait files, and multiple jobs reading the same storage can
be slower despite spare CPU/RAM. The package does not automatically tune worker
counts from RAM size, file sizes or disk speed. Reuse prepared inputs and LDSC
results when scientifically applicable; avoid repeating upstream work solely
to analyze a different selected subset.

### Why numerical libraries default to one thread

Standalone `prepare` sets the shared Polars pool from `--n_cores` before importing
Polars. This pool is shared by preparation workers; it is not a separate pool per
trait. As with GPCA export, an explicit `POLARS_MAX_THREADS` is respected. Python
API callers with Polars already initialized retain their existing pool. Polars
documents that this setting must be applied before initialization in its
[thread-pool reference](https://docs.pola.rs/api/python/stable/reference/api/polars.thread_pool_size.html).

In `prepare --mode both`, GPCA and LDSC retain independent filters and QC reports.
Their original-record reports share a streaming VCF read, and the transformed
GPCA frame is released before constructing the LDSC frame. Automatic chi-square
filtering still scans all complete LD-matched rows before filtering; repeated N
text reuses only the preceding validated value within that scan. No cutoff is
cached across files or runs.

The CLI sets the following variables to `1` **only when they are unset**, before
numerical imports and child launches:

| Variable | Controls |
| --- | --- |
| `OPENBLAS_NUM_THREADS` | OpenBLAS numerical-library threads |
| `OMP_NUM_THREADS` | OpenMP thread default |
| `MKL_NUM_THREADS` | Intel MKL numerical-library threads |
| `VECLIB_MAXIMUM_THREADS` | Apple Accelerate/vecLib numerical threads |
| `NUMEXPR_NUM_THREADS` | NumExpr expression-evaluation threads |

This reduces nested parallelism: for example, many LDSC/GWAMA workers should not
each start a large numerical thread team and compete for the same CPUs. It does
not reduce `--n_cores` to one. Existing values are respected, so shell settings
can override these defaults. Only libraries using a given variable are affected;
this is not a universal cap on every library thread.

### Faster chromosome combination and export

The exporter uses Polars and gzip compression as follows:

1. Read gzip chromosome files with up to
   `min(resolved n_cores, number of files)` concurrent readers. Each CSV reader
   has one parser thread; table operations share one Polars pool.
2. Restore run-status file order, concatenate, count Direction characters and
   sort by numeric CHR/BP. Equal-position ties keep their input order.
3. Write the combined table, then the selected summary, directly to gzip.
   N_eff/INFO constants are assigned before summary writing, without rereading
   and rewriting the entire output afterward.

`POLARS_MAX_THREADS`, if already set, is preserved. Otherwise explicit
`--n_cores` supplies its default before import; automatic mode selects up to 22
available CPUs. An already initialized Polars pool in a Python API process
cannot be resized; the actual pool size is recorded in the JSON audit.

With `pigz` on PATH, compression uses the resolved worker count. Without it,
Python gzip compresses with one worker; concurrent table reads still operate.
`--gzip_level 1` is the faster default; choose 9 when smaller output matters
more. Polars formatting and compression can overlap, so `--n_cores` is a
per-stage worker setting, **not a hard total process-thread limit**.

Source columns are read as **strings**, preserving the written decimal/scientific
notation, tiny P values, leading zeros and missing-value spellings. Only separate
CHR/BP sort keys are converted numerically. Compression is lossless. This avoids
rounding existing GWAMA numbers during export; it cannot restore precision
already lost upstream. Explicit override arguments are parsed as Python floats
and are not arbitrary-precision decimal inputs. See the installed-version
[Polars CSV documentation](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.read_csv.html)
and [thread-pool documentation](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.thread_pool_size.html).

The complete combined table still resides in memory, with additional parsing
and sorting buffers. No plain-text intermediate export is written to disk.
Check `read_workers`, `polars_threads`, `compression` and `timings_seconds` in
the export JSON to see where time was spent. Local synthetic timing/content
checks are documented in [export validation](../tests/POSTPROCESS_VALIDATION.md);
they do not predict runtime for a 10-million-variant run on shared storage.

## Raw LDSC commands

These two advanced launchers run the pinned upstream CBIIT scripts in the child
environment. All arguments after the script name are passed through unchanged.
They use the upstream spelling, including hyphenated flags; the managed tables
above do not define their interface.

For a whitespace-delimited GWAS table with `SNP,A1,A2,N,Z,P` and a HapMap
`SNP,A1,A2` allele reference:

```bash
ldsc-gpca munge_sumstats.py \
  --sumstats /data/trait.tsv.gz \
  --snp SNP --a1 A1 --a2 A2 --N-col N --p P \
  --signed-sumstats Z,0 \
  --merge-alleles /references/hm3_alleles.tsv \
  --out /results/trait
```

For one genetic-correlation comparison from already munged files:

```bash
ldsc-gpca ldsc.py \
  --rg /results/trait1.sumstats.gz,/results/trait2.sumstats.gz \
  --ref-ld-chr /references/eur_ld_chr/ \
  --w-ld-chr /references/eur_weights_chr/ \
  --out /results/trait1_trait2
```

| Option shown | Meaning |
| --- | --- |
| `--sumstats FILE` | Unmunged table for the raw munging script. |
| `--snp`, `--a1`, `--a2`, `--N-col`, `--p` | Explicit column mapping; values above match the example schema. |
| `--signed-sumstats Z,0` | Signed statistic column and its null value; here Z with null 0. |
| `--merge-alleles FILE` | HapMap allele-matching reference. |
| `--rg FILE1,FILE2,...` | First trait versus each remaining trait; not all combinations. |
| `--ref-ld-chr PREFIX` | Chromosome LD-score prefix; a directory prefix needs the trailing slash. |
| `--w-ld-chr PREFIX` | Chromosome regression-weight prefix. |
| `--out PREFIX` | Output filename prefix, unlike managed `--outdir`. |

For **all upstream Python options and defaults**, use the installed scripts' own help:

```bash
ldsc-gpca munge_sumstats.py --help
ldsc-gpca ldsc.py --help
```

The launchers add no analysis defaults to those scripts. They use `CONDA_EXE`
or `conda`, and `LDSC_GPCA_LDSC_PREFIX` when set, otherwise
`LDSC_GPCA_LDSC_ENV`, otherwise `ldsc-cbiit`. The setup script configures the
child automatically. Unlike managed help, raw help needs the child runtime.

Raw launchers bypass package manifest checks, extraction, retries, per-trait
chi-square filtering, completeness checks, provenance and numerical CSV
compilation. Raw `ldsc.py --rg ... --chisq-max VALUE` uses Python LDSC cross-product
filtering, not managed `--chisq_max`'s per-trait rule. For the complete managed
LDSC-to-GPCA workflow, use `ldsc-gpca ldsc`.

## Docker: install and run without host Conda

With Docker installed, build locally from the repository:

```bash
docker build --platform linux/amd64 -t ldsc-gpca:0.6.2 .
```

Then run a command explicitly; the image has no ENTRYPOINT:

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -e TMPDIR=/tmp \
  -v /absolute/host/data:/data:ro \
  -v /absolute/host/results:/results \
  ldsc-gpca:0.6.2 \
  ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --ldsc_results /data/ldsc_results.csv \
  --outdir /results/gpca_qc \
  --validate_only
```

Create the host result directory first and mount every required input/reference.
Manifests must contain paths visible **inside** the container, not host-only
paths. The example is validation-only; full GWAMA also needs the per-trait
inputs and the appropriate export policy. The tag above is built locally; no
registry-published image is implied. On macOS, Docker's allocated CPUs/RAM can
be lower than host capacity, and amd64 can require emulation on Apple Silicon.

See [DOCKER.md](../DOCKER.md) for mounting and environment details, and
[Nextflow / Google Batch](../examples/nextflow/README.md) for orchestration examples.
The latter includes local tests; live GCP execution has not been validated here.

## Troubleshooting

| Symptom | What to do |
| --- | --- |
| `ldsc-gpca` not found / old options shown | Activate the main environment; check `ldsc-gpca --version`; reinstall the package from the updated clone. |
| Installer refuses an environment | Inspect `conda env list` and existing installation records. Do not rerun expecting an overwrite; use a fresh target or a package-only update when appropriate. |
| Child LDSC launch fails | Inspect the configured child prefix/name and Conda executable. Test raw `ldsc.py --help`. |
| Manifest missing-column error | Confirm CSV commas, exact headers and the input mode's required fields. A tab-delimited manifest is not a CSV. |
| Missing VCF field | Check the module's fixed schema: `prepare` needs ES/SE; Python `ldsc` needs EZ/SI plus INFO/EUR. Do not simply rename biologically different fields. |
| Missing chromosome TSVs | Split mode needs chromosomes 1–22 for every retained trait. Use matching `nosplit` preparation and GPCA when whole-genome layout is intended. |
| Missing LDSC pairs | Generate every pair plus self-pairs. Use `ref=yes` for every trait in managed Python LDSC. A single raw multi-trait `--rg` is insufficient. |
| h2-scale ambiguity | Inspect populated columns and N/prevalence provenance. Select an explicit scientifically appropriate scale; column renaming is not conversion. |
| Non-positive-definite CTI | Investigate intercept estimates and inputs. Trait dropping does not guarantee a fix; the package does not repair CTI automatically. |
| LDSC/PCA validation passed, but GWAMA failed | Validation-only does not inspect GWAMA files. Check their filenames, nine-column schema, alignment and numerical inputs. |
| GWAMA finished, but export failed | Read the export error: INFO/N_eff may be absent, SNPIDs/positions invalid, or output paths already present. Do not assume the regressions failed or invent constants. Preserve completed source results. |
| Combination/adding constants is slow | Check actual Polars/pigz workers and stage timings in JSON; verify pigz is installed, use gzip level 1, and measure storage contention across concurrent runs. |
| More workers do not help | Verify the stage is parallel, available jobs and scheduler/container allocation. Shared disk throughput or per-worker memory can dominate. |

For more failure-stage details, see the [workflow reference](#contents).

## Help and reproducibility

```bash
ldsc-gpca --help
ldsc-gpca --version
ldsc-gpca pipeline --help
ldsc-gpca pipeline --ldsc_backend genomicsem --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca ldsc --ldsc_backend genomicsem --help
ldsc-gpca gpca --help
ldsc-gpca gpca --prepare_help
ldsc-gpca gpca --postprocess_help
ldsc-gpca genomicsem gpca --help
ldsc-gpca genomicsem gpca --prepare_help
ldsc-gpca genomicsem gpca --postprocess_help
```

`-h` is also available for help. Root `--version` prints the package version;
it is not a per-workflow option. Managed help does not open analysis data or
require an operational R analysis runtime. Use the installed CLI; a copied thin
R entry script without its bundled modules is insufficient.

For maintainers: Python and R help use the same formatter. R analysis help is
also shipped as text so the package CLI can display it without launching R.
After editing an R argument definition, regenerate and check these snapshots
from the repository root using an environment with the package's R dependencies:

```bash
python scripts/update_cli_help.py
python scripts/update_cli_help.py --check
```

Use `--rscript /path/to/Rscript` if Rscript is not on PATH. The check exits with
an error if either help snapshot differs from its actual R parser output.
This maintenance step does not change analysis defaults or upstream algorithms.

Keep the command, package version/repository commit, local modifications,
manifests, input/reference checksums, original and retained trait order, all
QC/removal reports, N/prevalence sources and GWAMA source used. The installer
saves environment/package records. The YAMLs are portable specifications,
not complete lockfiles; some transitive Conda/R dependencies are not pinned.
Record or lock those separately when exact restoration is required.

Verification is stage-specific: the package has synthetic/unit and selected
real-tool checks, but a full real-data LDSC-to-GWAMA analysis, every native
platform and live GCP execution have not all been validated. Documentation
examples use illustrative paths and do not imply that your data have been run.

Additional documentation: [workflow reference](#contents),
[interface naming](NAMING.md), [Docker](../DOCKER.md),
[Nextflow / Google Batch](../examples/nextflow/README.md), and
[export validation](../tests/POSTPROCESS_VALIDATION.md).

## Sources and license

- [Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html): PC1 procedure and GWAMA modification.
- [Pinned CBIIT LDSC source](https://github.com/CBIIT/ldsc/tree/6c673952cee74bd5c57aef1555a03b1c015399a0): Python LDSC implementation; derived from [Bulik-Sullivan LDSC](https://github.com/bulik/ldsc).
- [Pinned GenomicSEM source](https://github.com/GenomicSEM/GenomicSEM/tree/6b65ca5db39fdade08b0d811477be1cdd57b5039): GenomicSEM munging and LDSC; see [ldsc.R](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R) for GenomicSEM defaults and automatic rules.
- [bcftools manual](https://samtools.github.io/bcftools/bcftools.html): VCF querying/filtering.
- [Polars 0.20 CSV reader](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.read_csv.html) and [thread pool](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.thread_pool_size.html): string parsing and concurrency controls.
- [pigz manual](https://zlib.net/pigz/pigz.pdf) and [Python gzip documentation](https://docs.python.org/3.10/library/gzip.html): compression behavior.

The repository currently has no project-level LICENSE file; do not assume a
blanket redistribution license. Upstream components retain their own licenses
and author attribution. See [GWAMA provenance](../src/ldsc_gpca/r/vendor/README.md)
for the bundled source.
