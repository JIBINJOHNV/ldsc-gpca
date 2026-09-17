# Workflow reference

[Back to the starting guide](../README.md) · [Docker](../DOCKER.md) · [Nextflow / Google Batch](../examples/nextflow/README.md)

Detailed file contracts, defaults, QC policies and outputs. Start with the README
for installation and the shortest route to your first analysis.

## On this page

- [Installation details](#install-and-activate)
- [Prepare GPCA/GWAMA inputs](#prepare-per-trait-gpca-inputs)
- [Python LDSC](#run-python-ldsc)
- [GenomicSEM LDSC](#run-genomicsem-ldsc)
- [GPCA/GWAMA inputs and commands](#run-gpca-and-gwama)
- [QC, PC1 orientation and heritability scales](#understand-qc-and-scale-choices)
- [Outputs and automatic export](#find-and-interpret-the-outputs)
- [Troubleshooting and limitations](#troubleshooting-and-limitations)
- [Reproducibility and implementation](#reproducibility-and-implementation)

## Install and activate

This section is for native Conda installation. To use containers instead, jump to
[Docker installation and usage](../README.md#docker-install-and-run-without-host-conda).

### 1. Install once

With Conda already initialized in your terminal, run:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
bash scripts/setup_environments.sh --no-activate
```

The installer installs and checks Python, R, GenomicSEM, bcftools and CBIIT LDSC.
You do not need to install them separately. Keep the repository's YAML files and
`scripts/` directory together. Git and internet access are required.

### 2. Activate whenever you open a terminal

```bash
conda activate ldsc-gpca
ldsc-gpca --help
```

**Only activate this one environment.** The package automatically launches Python
LDSC in a separate environment named `ldsc-gpca-ldsc`. You never need to switch
between them. The separation keeps LDSC's pandas 1.5.0 dependencies independent
of the main package's pandas 2.x dependencies.

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

**Platform limits:** Linux x86_64/aarch64 and macOS Intel/Apple Silicon are detected,
but old dependency pins may not resolve on every platform. The installer stops
rather than changing pins or enabling emulation. Linux amd64 Docker installation
has been tested; a fresh native installation on every platform has not.

### What you still need to supply

- GWAS summary statistics, appropriate LD-score references and HapMap allele/SNP lists.
- For GWAMA: the reviewed modified v1.2.6 R function is bundled by default.
  `--source_path` optionally selects a trusted custom modified function.
- Compatible ancestry, genome build, allele orientation, phenotype definitions,
  and sample-size/prevalence conventions. Installation does not validate these.

Installing only with `python -m pip install .` installs the Python package, not R,
bcftools, GenomicSEM or the isolated LDSC environment. It is suitable only when you
manage those dependencies separately. Python 3.10–3.12 is required by this package;
the provided main environment uses Python 3.11.

## Prepare per-trait GPCA inputs

Create `prepare_traits.csv` (**comma-separated**, with a header):

```csv
traitname,vcf_files
protein1,/data/protein1.vcf.gz
protein2,/data/protein2.vcf.gz
```

Names must be unique/non-empty. Relative VCF paths resolve beside this manifest.
Each VCF must contain exactly one sample and FORMAT fields `AF,ES,SE,LP,NEF`.
`ES` and `AF` must refer to ALT; `LP` means `-log10(P)`.
VCFs must already have appropriate scientific QC and consistent build/alleles.
This command does **not** perform allele harmonisation, liftover, imputation-score
filtering or MHC filtering.

```bash
ldsc-gpca prepare \
  --input /data/prepare_traits.csv \
  --outdir /results/prepared \
  --splitby_chr split
```

The GPCA output is **tab-separated**, with exactly these columns in this order:

```text
SNPID CHR BP EA OA EAF N Z P
```

The header above lists names; spaces shown here are not the actual TSV delimiter.
The mapping is `EA=ALT`, `OA=REF`, `EAF=AF`, `N=NEF`, `Z=ES/SE`, and
`P=max(10^(-LP), p_min)` for valid LP. Verify that your supplied NEF convention is
appropriate: preparation does not infer total sample size or perform liability conversion.

| Option | Default | Other choice / meaning |
| --- | --- | --- |
| `--splitby_chr` | `split` | `nosplit`: one file per trait |
| `--gpca-id-source` | `chr_pos_ref_alt` | `vcf_id`: use existing VCF ID |
| `--prepare-workers` | `4` | Lower this when each trait uses substantial memory |
| `--bcftools` | `bcftools` on PATH | Executable path |
| `--p-min` | `1e-300` | Positive P floor below 1; adjustments are audited |
| `--write-munge-inputs` | Off | Also export unmunged HapMap-selected tables |
| `--hapmap-file` | Unset | Required with `--write-munge-inputs`; TSV with `SNP` header |
| `--munge-id-source` | `vcf_id` | `chr_pos_ref_alt`; must match the HapMap identifiers |

`chr_pos_ref_alt` constructs `CHR_POS_REF_ALT` after chromosome/position
normalization. `vcf_id` does not look up rsIDs. GPCA and munging identifier choices
are independent. A standard rsID list needs rsIDs in the VCF ID field.

To export optional munging tables, add `--write-munge-inputs --hapmap-file
/references/hm3_snps.tsv` to the command. This is **SNP selection only**, not allele
alignment. Extra allele columns in this list are ignored. GPCA outputs are not
restricted to HapMap. Zero HapMap matches stops the optional export.

| Output, under `--outdir` | Purpose |
| --- | --- |
| `gpca_inputs/{traitname}_chr{CHR}_GenomicPCA_inputs.tsv` | Default split output; chromosomes 1–22 |
| `gpca_inputs/{traitname}_GenomicPCA_inputs.tsv` | `nosplit` output |
| `munge_inputs/{traitname}_munge_inputs.txt` | Optional **space-separated, unmunged** table: `SNP,CHR,POS,A1,A2,eaf_A1,beta,se,N,p` |
| `Preparation_Status.csv`, `Preparation_Settings.json` | Execution status and settings |
| `GPCA_Input_QC_Summary.csv` | Per-trait retained/removed/adjusted counts and errors |
| `GPCA_Input_QC_Issues.csv` | Affected original VCF records, followed by `traitname,QC_action,QC_reason` |

Preparation keeps autosomes 1–22. Split mode requires at least one retained variant
on **every** chromosome for each trait, because the R loader expects all 22 files.
Invalid numerical values, non-positive SE/N, negative LP, invalid frequency or
position, invalid selected IDs, equal alleles and multiallelic/symbolic alleles are
removed at the variant level. Sequence indels may pass; this is not reference validation.
Duplicates keep the largest valid original LP, then the first original row on ties.
That rule is deterministic, not proof of superior quality.

Valid P values are not significance-filtered. Extreme LP values are floored on
conversion; Z=ES/SE is unchanged. In subsequent LDSC munging, Z magnitude is derived
from P, so a P floor can limit extreme Z magnitudes. It is a numerical policy,
not a scientifically optimal universal threshold.

The issues CSV preserves original field strings, not transformed output values;
metadata lines are excluded. Actions distinguish `removed`, `duplicate_removed`
and retained `p_adjusted` records. Different sample headers are combined using
their union. Source VCFs are not changed. If any trait fails, prepared tables are
not published; audit reports describe the failure. Unextractable traits have
unknown counts, not invented zeros. Use a fresh output directory to retry.

## Run Python LDSC

This route uses **local bcftools plus CBIIT LDSC in the isolated environment**.
Its VCF schema differs from `prepare`: INFO fields `AF,EUR` and FORMAT fields
`SI,AF,EZ,LP,NEF` are required by the extraction commands. Binary traits with
population prevalence additionally need `NC,NCO`. This is an EUR-field-specific
workflow, not a generic VCF reader. The flags do not provide ancestry/field mapping.

Create `python_ldsc_traits.csv` (**CSV**; all five headers are required for VCF input):

```csv
gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence
protein1,/data/protein1.vcf.gz,yes,NA,NA
protein2,/data/protein2.vcf.gz,yes,NA,NA
```

Use absolute VCF paths. `ref=yes` requests that trait against every listed trait,
including itself. Set it for **every trait** when you need complete GPCA coverage;
`ref=no` leaves a trait as a target only. The later GPCA manifest uses `traitname`,
not `gwas_name`; the identifier values must match.

```bash
ldsc-gpca ldsc \
  --input_file /data/python_ldsc_traits.csv \
  --output_folder /results/python_ldsc \
  --ld_ref_snp_file /references/hm3_alleles.tsv \
  --ld_ref /references/eur_ld_chr \
  --chisq-max 80 \
  --n_cores 5 \
  --ldsc-retries 1
```

`--ld_ref_snp_file` must have **whitespace-separated `SNP,A1,A2` headers**, unlike
the SNP-only preparation list. CBIIT munging checks allele compatibility with this
reference; do not equate this with whole-workflow harmonisation of your GWAMA files.
`--ld_ref` is a directory with chromosome `.l2.ldscore.gz` and associated M files.
This wrapper currently uses the **same directory for LD references and weights**;
it has no separate Python weight-directory flag. Ensure these files are appropriate
for both uses and for the study ancestry.

### Sample-size and prevalence rules

| Population prevalence | Extraction/munging N | Sample prevalence / LDSC conversion |
| --- | --- | --- |
| Blank/NA | `NEF` | Supplied sample prevalence is ignored; no conversion for that trait |
| Supplied, strictly between 0 and 1 | `NC + NCO` | Use supplied valid sample prevalence, or derive the median per-variant `NC/(NC+NCO)` |

If all population prevalences are absent, no prevalence flags are passed to LDSC.
If some are present, per-trait missing markers are used for unconverted traits.
The median case fraction is an implementation rule; verify its suitability for
cohort meta-analysis and variable per-SNP sample sizes. A population prevalence
alone does not establish that the available count fields are scientifically appropriate.

### Optional GenomicSEM-style chi-square filtering

`--chisq-max FLOAT` independently retains `Z^2 <= FLOAT` in every munged trait
after munging and before pairwise LDSC. It applies identically to VCF-generated
inputs and `--ldsc_only` inputs. The default is disabled, preserving previous
package behavior.

The package writes filtered copies under
`ldsc_input_chisq_filtered/{gwas_name}.sumstats.gz`; it never overwrites the VCF
or source munged files. `LDSC_ChiSquare_Filter_Summary.csv` records per-trait
before/removed/after counts and `LDSC_ChiSquare_Excluded_Variants.tsv.gz` records
the excluded `gwas_name,SNP,Z,CHISQ` values. `LDSC_Runtime.json` records the rule
and threshold.

This is a GenomicSEM-compatible per-trait rule, not a pass-through to CBIIT
LDSC's native `--rg --chisq-max`, which filters the cross-product using
`Z1^2 * Z2^2 < threshold^2`. The package deliberately runs pairwise LDSC on the
independently filtered copies without forwarding the native option. A fixed
threshold such as `--chisq-max 80` is explicit and reproducible; this option does
not implement GenomicSEM's automatic `max(80, 0.001*N)` threshold selection.

The compiler uses the **target trait (`p2`)** for heritability scale labels. It
writes both observed/liability column pairs, leaving unused entries empty. Its
`h2_scale=NEF_unconverted` label does not prove ordinary observed-scale binary
heritability. In particular, binary estimates based on effective N must not be
treated as compatible observed/liability values merely because of a column name.

### Important defaults and failure behavior

| Setting | Default |
| --- | --- |
| Local workers, `--n_cores` | `5` |
| MHC exclusion | Off; enable with `--exclude-mhc` |
| MHC interval if enabled | Chromosome 6, 25,000,000–35,000,000 |
| Extraction SI / MAF thresholds | `--info-min 0.7`, `--maf-min 0.01` |
| Additional munging MAF threshold | `--munge-maf-min 0.005` |
| Maximum INFO AF versus EUR difference | `--max-af-difference 0.2` |
| Extraction palindromic removal | Off; `--remove-palindrome` uses AF bounds 0.45–0.55 by default |
| Munging palindromic removal | CBIIT removes all palindromic SNPs, independently of the extraction flag |
| Per-trait chi-square filter, `--chisq-max` | Disabled; a supplied positive finite value keeps `Z^2 <= value` |
| `--ldsc-retries` | `1`: initial command plus one retry |

The P floor in **`prepare` does not apply to this extraction route**. Its existing
LP-to-P transformation can underflow for extreme values; P=0 is rejected during
LDSC munging. Review logs for resulting exclusions.

Retries rerun identical failed commands; they do not repair data or statistical
problems. Exhausted retries stop before compilation. Already submitted jobs finish
before exit. Pair-level failures detected during result validation are **not retried**.
Extraction failure stops the run; handled munging failures remove affected traits.
Compilation requires all requested comparisons and finite rg, but is not the full
GPCA QC described below.

`--ldsc_only --ldsc_input_folder /data/munged` reuses `{gwas_name}.sumstats.gz`
and reruns **all** requested batches. It is not selective failed-batch resume.
Extraction and munging filters are not reapplied: recorded mismatches fail and
unknown legacy settings warn. A requested `--chisq-max` is applied to separate
filtered copies after this validation. Total-N traits require matching
`.prevalence.json` sidecars and file hashes. In this mode, `vcf_files` may be
omitted from the manifest; `gwas_name,ref,pop_prevalence,sample_prevalence`
remain required. `--ld_ref_snp_file` is also unnecessary because no munging is
performed; `--ld_ref` remains required for the LDSC regressions.

```bash
ldsc-gpca ldsc \
  --input_file /data/python_ldsc_traits.csv \
  --output_folder /results/python_ldsc \
  --ld_ref /references/eur_ld_chr \
  --ldsc_only \
  --ldsc_input_folder /data/munged \
  --chisq-max 80 \
  --n_cores 5 \
  --ldsc-retries 1
```

Outputs include `ldsc_results.csv`, `LDSC_Runtime.json`,
`LDSC_Trait_Prevalence_Metadata.csv`, munged files/sidecars and logs.
Each managed LDSC batch also writes `ldsc_results/<batch>.results.csv` directly
from the native result DataFrame with `%.17g` formatting, before readable-log
rounding. The existing LDSC process performs this export; there is no second
regression job and no reconstruction of SEs, p-values or heritabilities.
The compiler reads each numerical CSV once and does not open the readable logs.
Missing/invalid exports stop the run rather than falling back to log parsing.
Both the original batch exports and the compiled CSV preserve numerical precision.
`LDSC_Runtime.json` records this output contract.

The installed package supplies a standalone exporter executed by the existing
isolated LDSC Python runtime. The hook observes only the native result-table
rendering, restores its temporary hooks even on failure, and leaves regression,
filtering, liability conversion and dependency files unchanged. No LDSC
environment reinstall or source patch is needed. Raw script passthrough below
intentionally remains native and does not perform this managed export.

Explicit `.log` inputs to the Python `compile_results()` API remain supported
for one-time legacy recovery; passing a missing CSV never triggers that path.
When `--chisq-max` is supplied, outputs also include
`LDSC_ChiSquare_Filter_Summary.csv`, `LDSC_ChiSquare_Excluded_Variants.tsv.gz`,
and `ldsc_input_chisq_filtered/`.
Command failures are logged in `execution_errors.log`. Old output files are not
deleted on failure; never mistake an older CSV for a successful new run.

For separately managed environments, use `--ldsc-env NAME` or
`--ldsc-env-prefix DIRECTORY` (mutually exclusive). An explicit name overrides the
saved `LDSC_GPCA_LDSC_PREFIX`; without a saved prefix the fallback is `ldsc-cbiit`.
`--conda-executable` defaults to `CONDA_EXE` or `conda` on PATH.
`--bcftools` defaults to `bcftools` on PATH and is not needed with `--ldsc_only`.

### Raw CBIIT script passthrough

The installed command can invoke the exact pinned CBIIT console scripts without
the managed pipeline:

```bash
ldsc-gpca munge_sumstats.py \
  --sumstats /data/trait.tsv.gz \
  --merge-alleles /references/hm3_alleles.tsv \
  --out /results/trait
```

```bash
ldsc-gpca ldsc.py \
  --rg /results/trait1.sumstats.gz,/results/trait2.sumstats.gz \
  --ref-ld-chr /references/eur_ld_chr/ \
  --w-ld-chr /references/eur_ld_chr/ \
  --out /results/trait1_trait2
```

All arguments after `ldsc.py` or `munge_sumstats.py` are passed through unchanged,
and the upstream exit status is returned. `--help` therefore displays the native
script help. The launcher uses `CONDA_EXE` when set, then `conda`; it uses
`LDSC_GPCA_LDSC_PREFIX` when set, otherwise `LDSC_GPCA_LDSC_ENV`, and finally the
legacy environment name `ldsc-cbiit`.

Raw passthrough intentionally performs no package manifest validation, VCF
extraction, GenomicSEM-style per-trait filtering, retries, completeness checks,
provenance recording or result compilation. In particular, raw
`ldsc-gpca ldsc.py --rg ... --chisq-max 80` uses native LDSC's cross-product rule.

## Run GenomicSEM LDSC

This is an integrated command, not an external prerequisite you must script
yourself. It runs local R/GenomicSEM and does not call Python LDSC or GWAMA.

Default input manifest (**CSV**):

```csv
traitname,munge_inputs,sampleprevalence,populationprevalence
protein1,/data/protein1.txt,NA,NA
disease1,/data/disease1.txt,0.2,0.05
```

| Input | Required format |
| --- | --- |
| Manifest | At least two unique trait names; no whitespace/path separators. Relative file paths resolve beside the manifest. |
| `munge_inputs` files | Whitespace-separated GWAS tables: `SNP,A1,A2,P`, signed effect `BETA` or `Z`, and `N`; positive manifest `N` can supply a per-trait constant instead. |
| `--hm3` | Whitespace-separated reference with `SNP,A1,A2`; native munging aligns effects to reference allele order. |
| `--ld` | Directory with `<CHR>.l2.ldscore.gz` and `<CHR>.l2.M_5_50`. |
| `--wld` | Optional separate directory of chromosome `.l2.ldscore.gz` regression weights; defaults to `--ld`. |

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/genomicsem_traits.csv \
  --hm3 /references/hm3_alleles.tsv \
  --ld /references/eur_ld_chr \
  --wld /references/eur_weights_chr \
  --outdir /results/genomicsem_ldsc \
  --cores 4
```

For quantitative traits, both prevalence entries should be blank/NA; for binary
liability conversion, both must be strictly between 0 and 1. Mixed rows are allowed,
but **partial prevalence pairs stop**. Unlike the Python route, this wrapper does
not derive case fractions or choose NC+NCO versus NEF for you. GenomicSEM performs
its own column interpretation and can transform recognized sample-size conventions;
inspect its logs rather than assuming every input N field is passed through unchanged.

To skip munging, add **one** of the following and omit `--hm3`:

- `--munge-output /data/munged`: the manifest needs
  `traitname,sampleprevalence,populationprevalence`; the directory must contain
  exactly one `{traitname}.sumstats.gz` or `{traitname}.sumstats` for each trait.
- `--munged-input`: put explicit file paths in a `traits` manifest column instead
  of `munge_inputs`.

Existing munged files must be tab-separated with `SNP,A1,A2,N,Z` headers; order may
vary. `prepare --write-munge-inputs` produces **unmunged** tables, not these files.
Native munging interprets recognized columns; check that interpretation and N/allele
conventions if using preparation exports as its inputs.

| Option | Default / choices |
| --- | --- |
| `--cores` | `1`; munging workers, not parallel LDSC |
| `--info-filter`, `--maf-filter` | `0.9`, `0.01`; filtering depends on recognized fields |
| `--chromosomes` | `22`; reads chromosomes 1 through the supplied value, allowed 1–22 |
| `--n-blocks` | `200`; GenomicSEM may override this for >18 traits |
| `--chisq-max` | Unset; uses GenomicSEM's automatic rule |
| `--invalid-h2-action` | `drop`; alternative `error` for non-positive/non-finite raw h2 |
| `--rscript` | `Rscript` on PATH |

The wrapper runs `stand=FALSE`, audits raw h2, then reruns `stand=TRUE` on retained
traits. Fewer than two retained traits or invalid final matrices stops the run.
Failed munging stops rather than silently removing traits. A fresh/empty output
directory is required. The final object contains genuine native sampling covariance
matrices; none are fabricated or repaired. Final RData validation does not replace
the downstream GPCA/CTI checks.

Important outputs:

- `genomicPCA_LDSC.RData`: final `LDSCoutput` for `genomicsem gpca`.
- `genomicPCA_LDSC_raw.RData`: preliminary object; **do not use it as GPCA input**.
- `Selected_Traits.csv`: retained trait order; usable as the GPCA manifest with
  existing GPCA files. Add `vcf_files` if automatic preparation is needed.
- `GenomicSEM_LDSC_Trait_QC.csv`, `GenomicSEM_LDSC_Events.csv`,
  `Resolved_Manifest.csv`, native logs and `sessionInfo.txt`.
- `munge_output/`: newly munged files in default mode.

## Run GPCA and GWAMA

### Common inputs

Both backends use a **CSV manifest** whose `traitname` column selects at least two
unique, non-empty traits. Row order controls matrices, PC1 weights and GWAMA input
order. Use exactly the identifiers from your LDSC results and GPCA filenames.
Extra LDSC traits are allowed; they do not have to be removed from the source file.

```csv
traitname
protein1
protein2
```

For GWAMA, supply `--gpca_input_folder` with tab-separated files containing exactly
`SNPID,CHR,BP,EA,OA,EAF,N,Z,P` in that order. The aliases `A1,A2,p` are accepted for
`EA,OA,P` in memory. Default split filenames are
`{traitname}_chr{CHR}_GenomicPCA_inputs.tsv` for all chromosomes 1–22;
`--splitby_chr nosplit` uses `{traitname}_GenomicPCA_inputs.tsv`.
The R reader checks the schema, **not complete row-level numerical QC**; externally
prepared files must already be validated and consistently aligned.

Both GPCA backends default to the bundled `N_weighted_GWAMA.function.1_2_6.R`,
under `src/ldsc_gpca/r/vendor/`. It retains the original author credits and the
reviewed genomicPCA modification unchanged. No `--source_path` is needed.
An optional `--source_path` override must define an already-modified `my_GWAMA`
or `multivariate_GWAMA`. Its weighting must use:

```r
W <- t(t(sqrt(N)) * h2)
sqrt_W <- W
```

Here the argument named `h2` receives **PC1 loadings**, not SNP heritabilities.
The bundled weighting was checked against the tutorial and synthetic calculations.
This is not validation of every inherited output formula or of custom sources.
The bundled source omits INFO: automatic export requires a scientifically justified
`--gwama-output-info` override; otherwise it stops after GWAMA. No INFO value is invented.

### Start with QC/PCA only

#### Python LDSC: QC/PCA example

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --python_ldsc /results/python_ldsc/ldsc_results.csv \
  --outdir /results/python_gpca_qc \
  --validate_only
```

The table can be CSV, TSV or whitespace-delimited text with a header; `.gz` is
accepted. Always required:

```text
p1 p2 rg se z p h2_int h2_int_se gcov_int gcov_int_se
```

Also supply `h2_obs,h2_obs_se` and/or `h2_liab,h2_liab_se` according to the scale
rules below. For k selected traits, require **k × (k + 1) / 2 unordered estimates**,
including all self-pairs. Triangular and symmetric tables are accepted. Consistent
opposite orientations are collapsed; conflicting duplicates stop. A single
`ldsc.py --rg A,B,C` command supplies A–B and A–C, not B–C.

#### GenomicSEM: QC/PCA example

```bash
ldsc-gpca genomicsem gpca \
  --input /results/genomicsem_ldsc/Selected_Traits.csv \
  --ldsc_path /results/genomicsem_ldsc/genomicPCA_LDSC.RData \
  --outdir /results/genomicsem_gpca_qc \
  --validate_only
```

This input is **binary RData**, containing `LDSCoutput` with
`S,V,I,S_Stand,V_Stand`. It is not a renamed CSV or the preliminary raw object.
This command reads existing estimates; it does not rerun GenomicSEM LDSC.

`--validate_only` checks LDSC/PCA and writes available diagnostics. It skips GWAMA
and VCF preparation, so it **does not validate your per-variant GWAMA input files**.

### Run SNP-level PC1 GWAMA after reviewing QC

Use a fresh output directory and omit `--validate_only`. With the bundled source,
first set `GWAMA_INFO` to a scientifically justified export value in `[0,1]`.
The shell guard below stops before analysis if it is unset; it is not an estimated
INFO score. See [the export policy](#automatic-gwama-export).

```bash
ldsc-gpca gpca \
  --input /data/selected_traits.csv \
  --python_ldsc /results/python_ldsc/ldsc_results.csv \
  --gpca_input_folder /results/prepared/gpca_inputs \
  --outdir /results/python_gpca_gwama \
  --dataset-id cluster1 \
  --gwama-output-info "${GWAMA_INFO:?Set a scientifically justified INFO value first}"
```

For native results, use `ldsc-gpca genomicsem gpca` and replace `--python_ldsc`
with `--ldsc_path /results/genomicsem_ldsc/genomicPCA_LDSC.RData`; the other options
are shared. Both backends default to correlation PCA, tutorial PC1 orientation,
split input and `--cores 0` (automatic chromosome-worker selection).

If GPCA files are absent, omit `--gpca_input_folder` and include `vcf_files` in
the GPCA manifest. The same preparation module runs first and writes
`<outdir>/gpca_inputs/`. Its options are shown by `--prepare-help`. All original
manifest traits must prepare successfully before R starts. With an existing folder,
non-default preparation settings are rejected. To reuse prepared files after an
R failure, explicitly supply that folder next time.

## Understand QC and scale choices

### Matrices and PC1

| Quantity | Python-table backend | Native GenomicSEM backend |
| --- | --- | --- |
| Correlation PCA matrix | Diagonal 1; off-diagonal rg | `S_Stand` |
| Covariance PCA matrix | Diagonal selected self h2; off-diagonal `rg × sqrt(h2_i × h2_j)` | Native `S`, on its supplied scales |
| GWAMA CTI | Diagonal self `h2_int`; off-diagonal `gcov_int` | Native `I` |

CTI is the LDSC intercept/error-covariance matrix used by GWAMA, **not** the genetic
correlation/covariance matrix. It must be symmetric and positive definite.
PCA uses a symmetric eigendecomposition. PC1 loadings are the first eigenvector
multiplied by the square root of its positive eigenvalue. Pairwise p, z and SE are
QC diagnostics, not PCA significance filters or precision weights.

All PCs are calculated, but only PC1 is used for GWAMA. The default
`--pc1_orientation tutorial` reverses the **whole** loading vector if its median is
negative; `as_computed` retains the arbitrary eigenvector sign. This changes the
PC's direction, not association strength. Interpret downstream effect signs and
genetic correlations relative to the recorded PC orientation, not as an intrinsic
biological “positive” direction.

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
SEs and result p in `[0,1]`; this **result-table** rule differs from munging's input
P requirement `0 < P <= 1`. Native standardized diagonal checks use absolute
tolerance `1e-8`, not the Python self-pair tolerance.

Removal is not guaranteed to retain the largest possible subset. Python removes
failed self traits first, then traits with the most remaining failed pairs; ties
use lower self h2/SE and then later manifest position. Native pair-removal ties
use manifest order. Removing traits changes the matrix and can change PC1/GWAMA;
inspect audits and compare sensitivity analyses rather than assuming negligible impact.
Structural errors and ambiguous scale selection are not solved by enabling removal.

Negative eigenvalues remain reported. The loading calculation uses
`pmax(eigenvalue, 0)` to avoid square roots of negative values; **the matrix itself
is not repaired**. A PC1 PASS is not a guarantee that the matrix or entire run passed.

### Heritability scales: Python-table backend

`--heritability_scale auto` examines populated values for selected traits, not just
column headers. Entirely empty scale pairs are ignored. Each estimate and SE must
come from the same scale. The p2 mapping is preserved; self-pairs are authoritative
for trait h2 QC and covariance construction.

| Policy | Behavior |
| --- | --- |
| Correlation + `auto` | Select the only populated scale per trait; different traits may use observed/liability scales for QC. Both populated for a trait is ambiguous and stops. rg and CTI are not rescaled. |
| Covariance + `auto` | Require exactly one complete common self-h2/SE scale across the selected traits, then finite positive values. Two complete scales or no complete common scale stops. |
| `observed` or `liability` | Explicitly select that scale for all traits; no fallback or conversion. |
| Covariance + `mixed` | Explicit opt-in to trait-specific scales; each trait still needs one unambiguous populated scale and complete valid self h2/SE. Warns about scale dependence. |

For covariance mode use `--pca_matrix covariance`. Covariance PCA is scale-dependent:
explicitly permitting mixed scales does **not** establish scientific comparability.
Verify phenotype, prevalence and N conventions, especially for unconverted binary
effective-N estimates. The reader does not infer prevalence or convert scales and
cannot verify those conventions from the h2 columns alone. Early covariance-scale
requirements can stop before trait-removal QC.

Native covariance PCA uses `S` directly, with no equivalent scale-conversion switch.
Native mixed binary/quantitative input likewise does not create common-unit covariance
simply by placing traits in one object.

## Find and interpret the outputs

### PC reports from both GPCA backends

These CSVs are written in `--outdir`, including validation-only runs, once input
validation reaches the eigendecomposition:

| File | What to inspect |
| --- | --- |
| `GenomicPCA_All_PCs_Variance_Explained.csv` | Every eigenvalue, raw/cumulative percentages, positive-normalized percentages and negative-eigenvalue flags |
| `GenomicPCA_PC1_QC.csv` | PC1 validity, direction convention and reasons; not overall run success |
| `GenomicPCA_PC1_Protein_Contributions.csv` | Retained traits in manifest order, signed coefficients/loadings, contribution percentages and ranks |
| `GenomicPCA_PC1_Weights_Used.csv` | PC1 loadings used by GWAMA and orientation metadata |
| `GenomicPCA_CTI_Used.csv`, `GenomicPCA_Correlation_Matrix_Used.csv` | Exact selected matrices |
| `GenomicPCA_Selected_Traits_Eigenvalues.csv` | Raw eigenvalues and loading-guard diagnostics |

Raw percentage is `100 × eigenvalue / sum(all eigenvalues)`. With negative eigenvalues,
this is **not a conventional non-negative variance partition**; values can be negative
or exceed 100%. Positive-normalized percentages use
`100 × max(eigenvalue,0) / sum(max(eigenvalues,0))` and are descriptive, not a repair.
The all-PC report leaves raw percentages undefined when the total is non-positive.
These percentages describe the estimated genetic matrix, not measured phenotypic variance.

Protein PC1 contribution is `100 × (PC1 eigenvector coefficient)^2`. It sums to 100%
for a valid PC1, is invariant to whole-vector sign reversal, and is not a causal
contribution or the final SNP-specific GWAMA weight, which also depends on N.

Python audits include `Python_LDSC_Input_Validation_Summary.csv`,
`Python_LDSC_Heritability_Scales.csv`, `Python_LDSC_Self_Pair_QC.csv`,
`Python_LDSC_Dropped_Missing_Traits.csv`, `Python_LDSC_Dropped_Failed_Traits.csv`,
and global heritability/correlation and SE-matrix reports.
Native audits include `GenomicSEM_QC_Events.csv`, `GenomicSEM_QC_Warnings.csv`,
`GenomicSEM_QC_Removed_Traits.csv`, `GenomicSEM_Trait_QC_Summary.csv`,
`GenomicSEM_Retained_Manifest.csv`, `GenomicSEM_LDSC_Used.RData` and settings/session files.
Native `Global_*` reports describe the retained selected set, not every unused RData trait.

Audits are stage-dependent: early input/scale errors can prevent later reports.
Native events are written on analysis errors once audit initialization succeeds;
Python reporting is not identical at every failure stage. A matrix-wide warning
does not identify a causal “bad protein.” `Problematic_Trait_Combinations.csv` is
no longer generated, and old files are not deleted automatically.

### Automatic GWAMA export

After successful GWAMA, the package reads only outputs listed as successful in
the **updated current-run** `GWAMA_Run_Status.csv`. Expected source results end in
`.N_weighted_GWAMA.results.txt.gz`; the external function must follow that convention.
Validation-only runs skip export.

| Output | Location |
| --- | --- |
| `{name}_GWAMA_combined_results.txt.gz` | `--outdir` |
| `{name}_GPCA_inputs.txt.gz` | `<outdir>/harmonisation_input/` |
| `{name}_postprocess.json` | `--outdir`; sources, row count and overrides |

`--dataset-id` supplies `{name}`; by default it is the output folder name. Both
compressed tables are tab-separated. The combined table is sorted by chromosome
and position and adds `count_question,count_plus,count_minus` from `Direction`.
The selected-column summary contains:

```text
SNPID CHR BP EA OA EAF N_eff BETA SE Z PVAL INFO
```

**Despite its filename, this is not the nine-column per-trait GWAMA input format.**
It is intended as input to a later, separate harmonisation workflow; this export
does not align alleles or convert genome builds.

| Export option | Default | Scope |
| --- | --- | --- |
| `--gwama-output-n-eff` | Preserve reported N_eff | Optional finite positive constant in the selected summary only |
| `--gwama-output-info` | Preserve reported INFO | Optional finite constant in `[0,1]` in the selected summary only |
| `--archive-chromosomes` | Off; keep originals | Move current-run source files/logs to `chromosome_wise/` after saving outputs |

Do not supply a universal N_eff or INFO constant simply to make a run pass. Without
an override the corresponding source column must exist; an explicit override can
supply an absent column but is **user-supplied metadata, not an estimate**.
Some external GWAMA scripts omit INFO, so export may fail after GWAMA succeeds.
Use an override only if scientifically justified for the downstream use, and retain
the audit. Overrides do not change the combined table, LDSC, PCA or GWAMA calculations.

The exporter checks required columns, non-empty/unique SNPIDs, positions, Direction,
current-run files and override ranges. It **does not comprehensively validate all
reported numerical fields** such as BETA, SE, PVAL, or preserved N_eff/INFO. Independently
QC the exported statistics before downstream analysis; successful export is not
scientific validation.

Existing export/archive paths are refused. Archive failure can leave saved outputs
and some moved originals; it is reported as an error. No automatic export retry
is performed. Use `--postprocess-help` for this section's CLI options.

## Troubleshooting and limitations

| Problem | What to check |
| --- | --- |
| Installer cannot solve dependencies | OS/CPU support, network and pinned package availability. Do not silently replace pins; retain the error and partial installation for diagnosis. |
| `ldsc-gpca` is not found | Run `conda activate ldsc-gpca`. For older/custom-prefix installations only, use their `activate.sh`. Installing from GitHub does not automatically update an existing environment. |
| Child LDSC cannot launch | Check the saved prefix, `--ldsc-env` / `--ldsc-env-prefix`, and Conda executable. Preflight help must run successfully. |
| Missing VCF fields | `prepare` and Python `ldsc` have different fixed schemas. A different field/ancestry needs appropriate input conversion, not relabelling without verification. |
| Missing chromosome files | Split GWAMA expects 1–22 for each trait; choose `nosplit` during both preparation and analysis if that layout is appropriate. |
| Missing/invalid LDSC pairs | Review upstream logs and complete pair coverage. Retries do not repair negative h2 or missing estimates; opt-in trait removal changes the analysis. |
| Ambiguous heritability scale | Inspect populated scales and provenance; choose a scale explicitly only when scientifically appropriate. |
| Non-positive-definite CTI | Investigate the intercept matrix/estimates. The package deliberately does not repair it automatically. |
| GWAMA succeeded but export failed | Check INFO/N_eff availability, current-run status, filename conventions and existing output paths. Do not assume GWAMA itself failed or invent replacement values. |

Neither SNP selection nor coordinate-based IDs establish allele alignment. Match
ancestry/build/alleles across traits and references before analysis. A validation-only
run, PC1 PASS or successful export cannot replace those checks.

Python pairwise output lacks GenomicSEM's full cross-estimate sampling covariance
matrices `V` and `V_Stand`. The package does not construct fake diagonal replacements.
Procedures requiring those matrices cannot be inferred from the Python table;
`paLDSC` is not run by either GPCA command. See the
[genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html)
for the method and its alternative covariance procedure, and validate its suitability
for your own traits rather than assuming tutorial reproduction establishes general validity.

## Reproducibility and implementation

The YAMLs are **portable specifications, not complete lockfiles**. Direct Python
analysis packages and upstream source revisions are pinned, but some Conda/R
dependencies and platform builds are not. The installer saves Conda explicit lists,
pip freezes and an R package/session report under its root. Also preserve:

- This repository's commit (`git rev-parse HEAD`) and any local modifications.
- Commands, manifests, all QC/removal reports and exact retained trait order.
- Input/reference checksums, ancestry/build, sample-size and prevalence sources.
- The bundled modified GWAMA source/version, or the exact custom `--source_path` used.

Conda explicit lists are platform-specific and do not capture GitHub/CRAN R package
installation completely. Version reports are provenance, not an exact R restore
mechanism; archive or lock those dependencies separately when required.

| Component | Location |
| --- | --- |
| Installer and specifications | `scripts/setup_environments.sh`, `scripts/install_genomicsem.R`, `environment.yml`, `environment.ldsc.yml` |
| Python CLI/preparation | `src/ldsc_gpca/cli.py`, `prepare.py`, `helptext.py` |
| Python LDSC workflow | `ldsc_cli.py`, `ldsc_runtime.py`, `extraction.py`, `munging.py`, `pairwise.py`, `results.py` under `src/ldsc_gpca/` |
| R launchers and export | `gpca.py`, `genomicsem.py`, `genomicsem_ldsc.py`, `postprocess.py` under `src/ldsc_gpca/` |
| Shared R method | `src/ldsc_gpca/r/shared/` |
| Backend-specific R readers/QC/workflows | `src/ldsc_gpca/r/python_ldsc/`, `src/ldsc_gpca/r/genomicsem/` |

Use the installed CLI. R entry points are bundled with their modules; copying a
single thin R script is insufficient. The Python reader and native reader remain
separate, while matrix/PCA/reporting/GWAMA helpers are shared. R help is bundled
from its parsers; parser changes need regenerated help and parity tests.

### Sources and attribution

- [Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html): method and GWAMA modification.
- [CBIIT LDSC source used by the installer](https://github.com/CBIIT/ldsc/tree/6c673952cee74bd5c57aef1555a03b1c015399a0): Python 3 LDSC fork; upstream [Bulik-Sullivan LDSC](https://github.com/bulik/ldsc).
- [GenomicSEM source used by the installer](https://github.com/GenomicSEM/GenomicSEM/tree/6b65ca5db39fdade08b0d811477be1cdd57b5039): native munging, LDSC and sampling covariance matrices.
- [bcftools manual](https://samtools.github.io/bcftools/bcftools.html): extraction/filter syntax.

### License

No redistribution license has been selected for this repository. Public availability
does not itself grant an open-source license. Third-party software retains its own licenses.
