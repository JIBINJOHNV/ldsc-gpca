# ldsc-gpca

Python LDSC and genomicPCA/GWAMA through a unified command-line interface.

* `ldsc-gpca prepare`: QC-filtered VCFs to GPCA inputs and optional LDSC munging tables.
* `ldsc-gpca ldsc`: VCF extraction, munging, and pairwise Python LDSC.
* `ldsc-gpca gpca`: genomicPCA/GWAMA from Python LDSC tables using the bundled R script.
* `ldsc-gpca genomicsem gpca`: genomicPCA/GWAMA from native GenomicSEM `LDSCoutput` RData.
* `ldsc-gpca genomicsem ldsc`: native GenomicSEM munging followed by LDSC, or LDSC from existing munged files.

This package does not replace LDSC, convert the R method to Python, or fabricate
GenomicSEM sampling covariance matrices. Both backends reuse neutral shared
PCA, reporting and GWAMA functions.
The existing workflows preserve their analysis defaults and P-value
conversion. Packaging tests are not a substitute for scientific validation on real data.

## Available commands and help

| Command | Available? | Starts from |
| --- | --- | --- |
| `ldsc-gpca prepare` | Yes | Single-sample VCF files listed in a CSV manifest |
| `ldsc-gpca ldsc` | Yes | VCF manifest, or existing Python LDSC munged files with `--ldsc_only` |
| `ldsc-gpca gpca` | Yes | Complete pairwise Python LDSC table |
| `ldsc-gpca genomicsem gpca` | Yes | Existing GenomicSEM `LDSCoutput` RData |
| `ldsc-gpca genomicsem ldsc` | Yes | Unmunged summary tables, or existing munged files |
| `ldsc-gpca genomicsem munge` | No standalone command | Munging is included in `genomicsem ldsc` by default |

Postprocessing is automatic after GWAMA; it is not a top-level `postprocess` command.
No-argument commands and `--help` display guidance without opening input files or
requiring R. Argument errors display full help and exit nonzero. Normal analyses
retain concise progress/error output; full help is not printed on every successful run.
Python help obtains defaults/choices from its parser. Full R help is generated from
the R parsers and bundled, so it remains accessible when R is missing. R parser
changes must regenerate `r/python_ldsc/help.txt` and `r/genomicsem/help.txt`;
local parity tests check that the generated text matches the actual R help.

GPCA help is scoped to the requested section (for both backends):

| Help option | Shows |
| --- | --- |
| `--help` or no arguments | Analysis/QC inputs and automatic GWAMA export options; one usage block |
| `--prepare-help` | Optional VCF preparation options, required VCF/manifest columns and separators |
| `--postprocess-help` | GWAMA export options only; no VCF preparation requirements |

Preparation options remain accepted on GPCA commands, but apply only when
`--gpca_input_folder` is omitted and the run is not `--validate_only`.
The standalone `prepare`, Python `ldsc`, and `genomicsem ldsc` commands retain
their own module-specific help. No unrelated module help is appended to them.

## Input files: columns and separators

| File / option | Separator or format | Required content |
| --- | --- | --- |
| Preparation `--input` | Comma-separated CSV, header required | `traitname,vcf_files`; names unique/non-empty; relative VCF paths resolve beside manifest |
| Python LDSC `--input_file` | Comma-separated CSV, header required | `gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence`; use absolute VCF paths; prevalence headers required even when blank |
| GPCA `--input`, either backend | Comma-separated CSV, header required | `traitname`; add `vcf_files` when using automatic preparation instead of an existing GPCA folder |
| Python GPCA `--python_ldsc` | CSV, TSV or whitespace-delimited text; `.gz` accepted; header required | `p1,p2,rg,se,z,p,h2_int,h2_int_se,gcov_int,gcov_int_se`, plus `h2_obs,h2_obs_se` and/or `h2_liab,h2_liab_se`; selected traits need every pair and self-pair |
| GenomicSEM `--ldsc_path` | Binary RData, **not** CSV/TSV | `LDSCoutput` containing `S,V,I,S_Stand,V_Stand`; generate upstream with `stand=TRUE` |
| `--gpca_input_folder` files | Tab-separated TSV, header required | Exact order: `SNPID,CHR,BP,EA,OA,EAF,N,Z,P`; `A1/A2/p` aliases accepted for `EA/OA/P` |
| Prepare `--hapmap-file` | Tab-separated text, header required | `SNP`; extra allele columns allowed; selection only, no allele alignment |
| LDSC `--ld_ref_snp_file` | Whitespace-separated (spaces or tabs), header required | `SNP,A1,A2` for standard LDSC `--merge-alleles`; not the preparation-only SNP list contract |
| LDSC `--ld_ref` | Directory, not a delimited table | Chromosome-prefixed `.l2.ldscore.gz` and associated M files; current workflow uses this directory for both LD references and weights |
| `--source_path` | R source file | Already-modified `my_GWAMA` or `multivariate_GWAMA`; required for GWAMA, not `--validate_only` |

For complete GPCA coverage, set `ref=yes` for every selected trait in the Python
LDSC manifest (`ref=no` leaves a trait as a target). Blank/NA prevalences are allowed.
Split GPCA filenames are `{traitname}_chr{CHR}_GenomicPCA_inputs.tsv`; whole-genome
filenames are `{traitname}_GenomicPCA_inputs.tsv`.

### Configurable options versus fixed contracts

Input/output paths, exposed filtering thresholds, parallelism, retries, ID choice,
P floor, PCA mode and PC1 orientation are supplied through CLI options. Their
defaults and allowed choices appear in command help; no personal paths are required.
Not every implementation setting is configurable: the Python LDSC workflow still
uses Docker image `jibinjv/ldsc:v3`, an EUR-specific VCF schema, and a shared
reference/weight directory. Preparation expects FORMAT `AF,ES,SE,LP,NEF`; LDSC
extraction expects INFO `AF,EUR` and FORMAT `SI,AF,EZ,LP,NEF` (plus `NC,NCO` for
population-prevalence traits). Autosomal/split naming and statistical invariants
are fixed contracts, not hidden filtering switches. This help update does not
introduce generic field mapping or change those scientific assumptions.

## Installation

Use Python 3.10–3.12 in an isolated environment (local validation used Python 3.10).
Clone the repository, then install the package:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
python -m venv .venv
source .venv/bin/activate
python -m pip install .
ldsc-gpca --help
ldsc-gpca prepare --help
ldsc-gpca ldsc --help
ldsc-gpca gpca --help
ldsc-gpca genomicsem gpca --help
ldsc-gpca genomicsem ldsc --help
```

For development, use `python -m pip install -e .` instead.
You can also run `python -m ldsc_gpca ...`.
Use the installed commands above for these workflows. The R script is packaged
under `src/ldsc_gpca/r/` and located automatically by `ldsc-gpca gpca`.

### External dependencies

* Preparation: local `bcftools` on PATH (or specify `--bcftools /path/to/bcftools`).
  This step uses Polars and local bcftools, not Docker. See the
  [official bcftools query documentation](https://samtools.github.io/bcftools/bcftools.html#query).
* LDSC: Docker on PATH, a running Docker daemon, and `jibinjv/ldsc:v3`.
  Docker must be able to mount all input/output/reference folders. Use absolute paths.
  The current extraction code assumes a POSIX host (Linux/macOS); native Windows
  execution is not supported. No image is pulled by package installation.
* GPCA: `Rscript` on PATH, R packages `argparse`, `data.table`, and `glue`
  (`parallel` is supplied with R). Install them in R with
  `install.packages(c("argparse", "data.table", "glue"))`.
* The modified Fürtjes GWAMA source must be supplied via `--source_path` along with
  any dependencies it requires. It is not included or altered by this package.
* Supply your own LD-score reference files, SNP list, and GWAS data.

Python dependency ranges are in `pyproject.toml`. For reproducible analyses, record
`pip freeze`, `sessionInfo()` from R, and the Docker image digest; the image tag is
not an immutable scientific environment.

## Native GenomicSEM munging and LDSC

Requires local R and an installed [GenomicSEM package](https://github.com/GenomicSEM/GenomicSEM).
This route does not use Docker, Python LDSC or an external GWAMA source.

To run munging first, supply a **comma-separated CSV** with these headers:

```csv
traitname,munge_inputs,sampleprevalence,populationprevalence
protein1,/data/protein1.txt,NA,NA
disease1,/data/disease1.txt,0.2,0.05
```

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/traits.csv \
  --hm3 /references/w_hm3.snplist \
  --ld /references/eur_ld_chr \
  --wld /references/eur_weights_chr \
  --outdir /results/native_ldsc \
  --cores 4
```

To **skip munging**, instead supply `--munge-output /data/munged` and omit `--hm3`.
The manifest then needs only `traitname,sampleprevalence,populationprevalence`;
the directory must contain exactly one `{traitname}.sumstats.gz` or `.sumstats`
for each selected trait. Alternatively use `--munged-input` with file paths in
a `traits` manifest column, matching the original LDSC script's format.
These two options are mutually exclusive. Without either, munging always runs.

| Input | Format and required columns |
| --- | --- |
| Manifest | CSV; `traitname,sampleprevalence,populationprevalence`, plus `munge_inputs` in default mode or `traits` with `--munged-input` |
| Unmunged GWAS | Whitespace-delimited; `SNP,A1,A2,P`, signed effect (`BETA` or `Z`), and `N`; optional positive manifest `N` overrides file N |
| `--hm3` | Whitespace-delimited reference with `SNP,A1,A2`; native munging performs reference-allele alignment, not just SNP selection |
| Existing munged files | Tab-separated; `SNP,A1,A2,N,Z` (column order can vary) |
| `--ld` | Directory with `1.l2.ldscore.gz` through the selected last chromosome and corresponding `.l2.M_5_50` files |
| `--wld` | Optional separate directory of chromosome `.l2.ldscore.gz` weights; omitted means use `--ld` |

Trait names must be unique, non-empty and contain no whitespace or path separators.
Relative manifest file paths resolve beside the manifest; trait order is preserved.
For quantitative traits, both prevalence entries must be blank/NA. For binary
liability conversion, both must be strictly between zero and one. Mixed rows are
allowed, but partial prevalence pairs stop: this route does not infer sample
prevalence or change supplied sample sizes. Ensure N and sample prevalence are
appropriate for the study's sampling/meta-analysis design. Native covariance S
retains the resulting observed/liability scales; it is not a common-unit covariance
matrix merely because all traits are in one file.

| Option | Default | Meaning |
| --- | --- | --- |
| `--cores` | `1` | Munging workers; does not parallelize LDSC |
| `--info-filter` | `0.9` | Native munging INFO threshold when the field is present |
| `--maf-filter` | `0.01` | Native munging MAF threshold when the field is present |
| `--chromosomes` | `22` | Read chromosomes 1 through this value; allowed 1–22 |
| `--n-blocks` | `200` | Requested jackknife blocks; GenomicSEM overrides this for >18 traits |
| `--chisq-max` | Unset | Use GenomicSEM automatic rule; otherwise a finite positive threshold |
| `--invalid-h2-action` | `drop` | Audit/drop non-positive or non-finite raw h2; alternative `error` |
| `--rscript` | `Rscript` | Executable name or path |

GenomicSEM's own logs report column interpretation, missing INFO/MAF filtering,
and any internal parameter changes. See its authoritative
[munge implementation](https://github.com/GenomicSEM/GenomicSEM/blob/master/R/munge.R)
and [LDSC implementation](https://github.com/GenomicSEM/GenomicSEM/blob/master/R/ldsc.R).

The wrapper first runs `stand=FALSE`, saves raw results and audits h2. It then
runs `stand=TRUE` on retained traits to obtain genuine `V` and `V_Stand`.
Fewer than two retained traits or invalid final matrices stop the run; no final
RData is published on failure. No covariance matrices are synthesized or repaired.
Use a fresh/empty output directory; existing results are never overwritten.

Outputs under `--outdir`:

* `munge_output/`: newly generated munged files (default mode only).
* `Resolved_Manifest.csv`, `Selected_Traits.csv`: resolved inputs and retained trait order.
* `GenomicSEM_LDSC_Trait_QC.csv`: each trait's raw h2, scale, action and reason.
* `GenomicSEM_LDSC_Events.csv`: warnings/errors/completion once R pipeline starts.
* `genomicPCA_LDSC_raw.RData`: first-pass results.
* `genomicPCA_LDSC.RData`: validated `LDSCoutput` for `genomicsem gpca --ldsc_path`.
* Native log files and `sessionInfo.txt` for provenance.

Use `Selected_Traits.csv` as the GPCA manifest with an existing GPCA input folder;
add matching `vcf_files` if automatic VCF preparation is needed. This command does
not itself run GPCA/GWAMA. Failed munging stops instead of silently removing traits.

## Prepare inputs from QC-filtered VCFs

Use this command when per-trait GPCA inputs do not already exist:

```bash
ldsc-gpca prepare \
  --input /absolute/path/prepare_traits.csv \
  --outdir /absolute/path/prepared \
  --splitby_chr split \
  --gpca-id-source chr_pos_ref_alt \
  --prepare-workers 4
```

The manifest requires these columns; row order and trait names are preserved:

```text
traitname,vcf_files
protein1,/absolute/path/protein1.vcf.gz
protein2,/absolute/path/protein2.vcf.gz
```

Relative VCF paths are resolved against the manifest directory. Each VCF must
contain exactly one GWAS sample and FORMAT fields AF, ES, SE, LP, and NEF. Inputs
must already be QC-filtered, biallelic, consistently aligned across traits, and on
the same genome build. ES and AF must refer to ALT. This command does not perform
allele harmonisation, liftover, imputation-score filtering or MHC filtering.

| Option | Default | Purpose |
| --- | --- | --- |
| `--gpca-id-source` | `chr_pos_ref_alt` | GPCA SNPID values; alternative: `vcf_id` |
| `--write-munge-inputs` | Off | Also write HapMap-filtered munging input tables |
| `--hapmap-file` | None | Required with `--write-munge-inputs`; tab-delimited with a `SNP` column |
| `--munge-id-source` | `vcf_id` | Munging SNP values; alternative: `chr_pos_ref_alt` |
| `--splitby_chr` | `split` | Per-chromosome files; `nosplit` writes one file per trait |
| `--prepare-workers` | `4` | Concurrent preparation workers |
| `--bcftools` | `bcftools` | Local executable name or path |
| `--p-min` | `1e-300` | Floor for P derived from valid LP; retain and report adjusted variants |

To also produce LDSC munging tables, add:

```bash
--write-munge-inputs \
--hapmap-file /absolute/path/w_hm3.snplist \
--munge-id-source vcf_id
```

The two identifier choices are independent. `chr_pos_ref_alt` constructs
`CHR_POS_REF_ALT` after removing a `chr` prefix and normalising chromosome/position
numbers. `vcf_id` uses VCF ID unchanged; it does not look up or convert rsIDs.
The munging identifier must match the HapMap `SNP` values. A standard rsID list
therefore needs rsIDs in VCF ID; coordinate IDs require a coordinate-ID list.
Zero HapMap matches is an error. GPCA variants are **never restricted to HapMap**.

Outputs:

* `<outdir>/gpca_inputs/{traitname}_chr{CHR}_GenomicPCA_inputs.tsv` for split mode.
* `<outdir>/gpca_inputs/{traitname}_GenomicPCA_inputs.tsv` for nosplit mode.
* `<outdir>/munge_inputs/{traitname}_munge_inputs.txt` only when requested.
* `<outdir>/Preparation_Status.csv` and `Preparation_Settings.json` for audit.
* `<outdir>/GPCA_Input_QC_Issues.csv`: original affected VCF records plus only `traitname`, `QC_action`, `QC_reason`.
* `<outdir>/GPCA_Input_QC_Summary.csv`: per-trait input, retained, removed and P-adjustment counts, plus success/error status.

GPCA columns are exactly `SNPID CHR BP EA OA EAF N Z P`, tab-delimited.
EA=ALT, OA=REF, EAF=FORMAT/AF, N=FORMAT/NEF, Z=ES/SE and P=max(10^(-LP), p_min).
Munging tables are space-delimited with `SNP CHR POS A1 A2 eaf_A1 beta se N p`.
They are **unmunged tables**, not `.sumstats.gz`; no LDSC command is executed by
`prepare`. This export preserves your supplied NEF convention. It does not derive
NC+NCO, apply prevalence conversion or use `--gwama-output-n-eff` for input N.
Confirm that NEF is appropriate for your GWAS and downstream analysis.

Only chromosomes 1–22 are retained. Split mode requires at least one variant on
each chromosome per trait because the current R loader expects all 22 files.
Invalid/missing numeric fields, SE <= 0, N <= 0, negative LP, invalid frequency,
invalid positions, equal alleles, multiallelic/symbolic alleles and non-finite Z
remove the affected rows, not the entire trait. Invalid selected GPCA identifiers
are also removed. Sequence alleles may include N or sequence indels; this is not
reference validation or allele alignment.
Duplicates keep the largest valid original LP (smallest P), then the first original
row on ties. This deterministic selection policy is not evidence of better quality.
Valid P-values are not significance-filtered. Extreme LP is capped during conversion
to prevent underflow; the resulting P floor defaults to 1e-300. This is a numerical
policy, not an LDSC-required threshold: [LDSC accepts 0 < P <= 1](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py).
Z=ES/SE is unchanged; subsequent LDSC munging derives Z magnitude from P, so the
floor limits extreme Z magnitudes in that downstream step.

The issues CSV preserves original VCF field strings (including INFO, FORMAT and
the sample column), not transformed GPCA values. It appends only `traitname`,
`QC_action`, `QC_reason`; each affected original record occurs once with combined
reasons. Actions are `removed`, `duplicate_removed`, or `p_adjusted` (retained).
VCF metadata lines are not variant records and are not copied. For different sample
column names, the combined CSV uses the union of original headers; non-applicable
cells are empty. Original VCF files are never modified. Plain and gzip VCFs are supported.
Reports survive later preparation failures; traits that cannot be extracted have
unknown (blank) QC counts and an error in the summary, not fabricated zero counts.
Issue reporting reads only the header when there are no affected records and stops
after the last affected record otherwise. Combined issue reports are streamed instead
of loading all traits' reports into memory. Inputs without duplicate identifiers skip
the LP-ranking sort; duplicate selection and the final genomic ordering are unchanged.
Summary trait names are read as text, preserving identifiers such as `001`.

Each worker uses a unique temporary TSV. If any trait fails, no prepared tables are
published and the status CSV lists failures. Unreadable/schema-invalid inputs, no
usable variants, missing required chromosomes, or invalid optional munging identifiers
still fail clearly. Existing output folders/audit files
are never overwritten; use a fresh outdir after a failure or for a different set
of settings. Each worker loads one extracted trait into memory; lower
`--prepare-workers` for large datasets.

## Python LDSC

```bash
ldsc-gpca ldsc \
  --input_file /absolute/path/traits.csv \
  --output_folder /absolute/path/ldsc_output \
  --ld_ref_snp_file /absolute/path/w_hm3.snplist \
  --ld_ref /absolute/path/eur_w_ld_chr \
  --n_cores 5 --ldsc-retries 1
```

Required manifest headers: `gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence`.
Names must be unique/non-empty. Blank prevalence values are allowed, but the headers
are required. `ref=yes` requests that trait against every listed trait, including
itself. For a complete matrix for GPCA, set `ref=yes` for **every selected trait**.

VCFs must contain the fields consumed by the original extraction filters:
INFO/AF, INFO/EUR, FORMAT/SI, FORMAT/AF, FORMAT/EZ, FORMAT/LP and FORMAT/NEF.
Traits with population prevalence additionally require FORMAT/NC and FORMAT/NCO.
Missing population prevalence keeps NEF and does not use sample prevalence.
Supplied population prevalence selects NC+NCO and fills missing sample prevalence
with the median per-variant case fraction. If all population prevalences are blank,
no prevalence flags are passed to LDSC.

MHC exclusion defaults to off. Extraction MAF defaults to 0.01 and munging MAF
to 0.005. Standard LDSC removes all palindromic SNPs in munging regardless of the
optional extraction-stage palindromic filter. See `ldsc-gpca ldsc --help`.

### Failure handling and reuse

`--ldsc-retries 1` means one initial attempt plus one retry, with identical inputs
and settings. Successful batches are not repeated during this run. Exhausted retries
stop the workflow before compilation; other submitted jobs finish before exit.
Validation failures are not retried. Retries do not repair files or statistical issues.
Execution failures are appended to `execution_errors.log` with attempt labels.

`--ldsc_only` reuses munged files and reruns **all** requested LDSC batches, not just
failed batches. Filters are not reapplied; recorded filter mismatches fail and unknown
legacy filter settings cause a warning. No selective cross-run batch resume is provided.
Extraction failures stop the run; handled munging failures remove affected traits.
Missing/malformed results or non-finite rg fail compilation. Existing result files
from earlier runs are not deleted on failure and must not be mistaken for new output.

Outputs include `ldsc_results.csv`, `LDSC_Trait_Prevalence_Metadata.csv`, per-trait
`.prevalence.json` sidecars, munged files, and LDSC logs. The compiled table reports
target-trait heritabilities and their scale. Quantitative and binary NEF-unconverted
estimates should not be assumed to be on a common liability scale.

## genomicPCA/GWAMA

Complete example, including a dataset ID and optional GWAMA-summary overrides:

```bash
ldsc-gpca gpca \
  --input /absolute/path/selected_traits.csv \
  --python_ldsc /absolute/path/ldsc_output/ldsc_results.csv \
  --gpca_input_folder /absolute/path/gpca_inputs \
  --source_path /absolute/path/my_GWAMA_26032020.R \
  --outdir /absolute/path/gpca_output \
  --splitby_chr split \
  --dataset-id cluster1 \
  --gwama-output-n-eff 330000 \
  --gwama-output-info 0.9
```

Post-processing runs automatically after successful GWAMA. With this example:

* Combined results: `/absolute/path/gpca_output/cluster1_GWAMA_combined_results.txt.gz`
* Selected-column summary: `/absolute/path/gpca_output/harmonisation_input/cluster1_GPCA_inputs.txt.gz`
* Override audit: `/absolute/path/gpca_output/cluster1_postprocess.json`

The values `330000` and `0.9` are example overrides, **not defaults**. Omit those
two options to preserve reported N_eff and INFO. They never change LDSC/GPCA
calculations or filtering.

The GPCA manifest is a separate CSV requiring `traitname`; its row order defines
matrix and GWAMA order. The LDSC file must supply every selected self-pair and pair,
with rg/diagnostic/intercept columns and supported observed or liability heritability
columns. Extra traits are allowed. The bundled R help describes QC, scale selection,
correlation/covariance modes, input file naming, and audit outputs in detail.

### Heritability scale selection for Python LDSC

`--heritability_scale auto` (default) inspects populated values **after selecting
manifest traits**, not merely column names. Entirely empty observed/liability
column pairs are ignored. Estimates and SEs are always taken from the same scale;
negative/non-finite estimates and missing/non-positive SEs still fail QC. The
reader never switches scale to rescue a failed QC value and never converts scales.
Python LDSC's pairwise heritability fields describe **p2**; self-pair values remain
authoritative for each trait's heritability QC and covariance reconstruction.

| PCA mode / policy | Behavior |
| --- | --- |
| Correlation, `auto` | Select the only populated scale per trait. Observed and liability traits may coexist for QC; `rg` and the intercept matrix are unchanged. Both scales populated for one trait is ambiguous and stops. |
| Covariance, `auto` | Select the only common scale with populated self heritability and SE for every selected trait; then require finite positive values. No complete common scale, or two complete common scales, stops. |
| `observed` or `liability` | Explicitly select that scale for all traits; no fallback to the other scale. |
| Covariance, `mixed` | Explicit opt-in to trait-specific scales. Each trait must have only one populated scale and complete finite positive self heritability/SE. Ambiguous traits still stop. A warning records that covariance PCA depends on those scales. |

For explicit mixed covariance, add `--pca_matrix covariance --heritability_scale mixed`.
This choice does not establish scientific comparability: the user must verify the
phenotype, prevalence and sample-size conventions behind each estimate. In particular,
an unconverted binary estimate calculated with effective N is not automatically an
ordinary observed-scale estimate. This table alone cannot verify those conventions.
No phenotype prevalence is inferred and no observed/liability conversion is attempted.

`Python_LDSC_Heritability_Scales.csv` records the selected scale for each manifest
trait once scale resolution succeeds, before later trait QC. Per-trait QC, removal
and global heritability reports also retain their actual scales. Derived covariance
pair reports include `Heritability_Scale_1` and `Heritability_Scale_2`. An ambiguous
scale selection stops before analysis rather than guessing a scale or dropping a trait.
This policy changes the Python-table reader only; native GenomicSEM uses its supplied S.

### Automatic preparation when GPCA inputs are absent

Omit `--gpca_input_folder` from the GPCA command to run the same preparation module
first. In that case, `--input` must contain both `traitname` and `vcf_files` as shown
above. Prepared files are saved to `<outdir>/gpca_inputs/` and passed to R automatically.
All selected manifest traits must prepare successfully before R starts. You can use
the preparation options above (including the two identifier settings and optional
munging export) directly with `ldsc-gpca gpca` in this automatic mode.

If `--gpca_input_folder` is supplied, the existing files are used and checked by R;
no conversion takes place. Non-default preparation settings are rejected rather
than silently ignored. `--validate_only` skips preparation and GWAMA. To reuse
previously generated inputs after an R failure, explicitly supply their folder on
the next invocation. No preprocessing or scientific changes are made to the R script.

The externally supplied GWAMA function must already use the Fürtjes modification:
`W <- t(t(sqrt(N)) * h2)` and `sqrt_W <- W`, where its `h2` argument receives PC
loadings. The package does not modify this function. Python pairwise LDSC output
does not supply GenomicSEM's full `V`/`V_Stand`, so it cannot support paLDSC here.

### Automatic GWAMA post-processing

After successful GWAMA, results are combined automatically. No `--postprocess`
flag is needed. The selected-column summary is always saved to `<outdir>/harmonisation_input/`.

| Setting | Default | Purpose |
| --- | --- | --- |
| Post-processing | Automatic after successful GWAMA | No `--postprocess` flag required |
| `--dataset-id` | Name of the `--outdir` folder | Prefix for exported filenames; replaces `--postprocess-name` |
| Summary folder | `<outdir>/harmonisation_input/` | Fixed location; no separate folder option |
| `--gwama-output-n-eff` | Preserve reported N_eff | Optional N_eff override in the exported summary only |
| `--gwama-output-info` | Preserve reported INFO | Optional INFO override in the exported summary only |

The optional dataset and summary-value settings are:

```bash
--dataset-id cluster1 \
--gwama-output-n-eff 330000 \
--gwama-output-info 0.9
```

Post-processing always runs after successful GWAMA. It requires successful R execution
and a successful, updated `GWAMA_Run_Status.csv`. Validation-only runs skip it.
Only current-run output prefixes in that status file are read; other matching files
in the folder are ignored. The expected source filename suffix is
`.N_weighted_GWAMA.results.txt.gz`; custom GWAMA functions must use that convention.

| Output | Location |
| --- | --- |
| `{name}_GWAMA_combined_results.txt.gz` | `--outdir` |
| `{name}_GPCA_inputs.txt.gz` | `<outdir>/harmonisation_input/` (always) |
| `{name}_postprocess.json` (sources, row count, overrides) | `--outdir` |

`--dataset-id` sets `{name}` and defaults to the name of the `--outdir` folder. Both compressed
tables are tab-delimited. Results are sorted numerically by chromosome and position;
the combined table includes `count_question`, `count_plus`, and `count_minus` from
`Direction`. Missing required columns, empty results, invalid direction strings or
positions, duplicate SNPIDs, missing files and unchanged old results cause an error.

The selected-column summary contains:

```text
SNPID CHR BP EA OA EAF N_eff BETA SE Z PVAL INFO
```

`--gwama-output-n-eff` and `--gwama-output-info` are optional, finite overrides,
**not automatic defaults**. They affect only the exported GWAMA selected-column
summary, never LDSC/GPCA calculations or variant filtering.
N_eff must be positive and INFO must lie in [0,1]. Without an override, that column
must exist and its reported values are preserved. An override can also supply an
absent column. Overrides affect only the selected-column summary, not the combined
results. They are user-supplied constants, not values estimated by this operation.
Despite its filename, this summary is not the nine-column per-trait GWAMA input
format and does not perform allele harmonisation or genome-build conversion.

Original files are kept by default. Add `--archive-chromosomes` to move the current
run's source results and matching logs into `--outdir/chromosome_wise/` only after
both tables and the audit have been saved. Existing output/archive paths are never
overwritten. An archive failure leaves saved outputs in place and is reported as an
error; some originals may already have moved. Choose a new name/output directory
for repeat exports. No automatic post-processing retry is performed.

Show these options without requiring R:

```bash
ldsc-gpca gpca --postprocess-help
```

## GPCA from native GenomicSEM LDSC output

This command accepts the `genomicPCA_LDSC.RData` produced with `stand=TRUE`,
containing an object named `LDSCoutput` with `S`, `V`, `I`, `S_Stand`, `V_Stand`.
Do not supply the preliminary `_raw.RData` file. This release adds the native
GPCA entry point only; GenomicSEM munging/LDSC estimation still runs externally.
It uses local R, not Docker. The adapter needs the same R packages as `gpca`;
GenomicSEM itself is required upstream to generate the native results.

```bash
ldsc-gpca genomicsem gpca \
  --input retained_traits.csv \
  --ldsc_path genomicPCA_LDSC.RData \
  --gpca_input_folder prepared/gpca_inputs \
  --source_path /path/to/modified_GWAMA.R \
  --outdir genomicsem_gpca_results \
  --pca_matrix correlation \
  --pc1_orientation tutorial
```

The CSV requires unique, non-empty `traitname`; its order controls the analysis.
Extra RData traits are allowed. Omit `--gpca_input_folder` to reuse automatic VCF
preparation (manifest then also requires `vcf_files`). All automatic postprocessing
options, including `--dataset-id`, `--gwama-output-n-eff` and `--gwama-output-info`,
are shared with `gpca`. Outputs use the same `harmonisation_input/` convention.

| QC / setting | Default and action |
|---|---|
| Absent manifest traits | Stop; `--allow_missing_traits` permits audited removal. |
| Invalid heritability, SE, self intercept or pair estimates | Stop; `--failed_ldsc_action drop_traits` permits audited removal. |
| Invalid pairs with removal enabled | Remove the trait involved in most remaining invalid pairs; ties use manifest order. Deterministic heuristic, **not** guaranteed maximum retention. |
| h2/SE | `--h2_z_warn_threshold 2`: warning only; `0` disables. |
| Finite abs(rg)>1 | `--rg_out_of_range_action warn`; optional `error`. Never clamp or remove solely for this warning. |
| Standardized diagonal | Require approximately 1 (absolute tolerance `1e-8`); stop on inconsistency. Not a Python self-pair rg test. |
| Structural matrix errors / inconsistent S_Stand / non-PD CTI | Stop; no matrix repair or automatic structural-error removal. |
| Negative eigenvalues | `--negative_eigen_action warn`; optional `error` for substantive negatives. `pmax` affects loadings only. |
| PC1 direction | `--pc1_orientation tutorial`: flip the entire vector when its median is negative; `as_computed` disables orientation. |
| PCA matrix | `--pca_matrix correlation` uses S_Stand; `covariance` uses native S without scale conversion. Mixed trait scales affect covariance PCA interpretation. |
| Validation only | `--validate_only` writes native QC/PCA outputs without GWAMA or VCF preparation. |

Use a fresh output directory for each attempt, including after a failed run.
QC outputs are written on normal completion **and on analysis errors** once the
output directory/audit is initialized:

* `GenomicSEM_QC_Events.csv`: severity, scope, trait, other trait, reason, value, action.
* `GenomicSEM_QC_Warnings.csv`: warning events, including retained-trait warnings.
* `GenomicSEM_QC_Removed_Traits.csv`: actual removals with reasons; multiple reasons may produce multiple rows per trait.
* `GenomicSEM_Trait_QC_Summary.csv`: one row per original manifest trait, retention/removal flags and combined reasons. Retention does not imply the overall run succeeded.
* `GenomicSEM_Retained_Manifest.csv`: selected manifest in analysis order after successful QC.
* `GenomicSEM_LDSC_Used.RData`: native selected matrices, including fully preserved/reordered cross-estimate `V` and `V_Stand` (not refitted).
* `GenomicSEM_Run_Settings.rds`, `GenomicSEM_SessionInfo.txt`: settings and R environment.

Pairwise warnings identify both traits. Matrix-wide or external-function warnings
have no trait attribution; they do not establish that a particular trait caused
the problem. Empty audit tables retain their column headers. `Global_*` diagnostics
in this adapter describe the **retained selected set**, not unused RData traits.
Existing PC1, eigenvalue, matrix and GWAMA status audit filenames are retained.
The nine-column R reader checks schema, not complete row-level QC; supply QC-filtered
GWAMA inputs. `paLDSC` remains disabled. See the
[GenomicSEM implementation](https://github.com/GenomicSEM/GenomicSEM/blob/master/R/ldsc.R)
for native matrix definitions and sampling-covariance ordering.

### Shared PC1-focused reports (both LDSC input types)

All PCs are calculated for context, but **only PC1 is used for GWAMA**. Both
`ldsc-gpca gpca` and `ldsc-gpca genomicsem gpca` write these files in `--outdir`,
also with `--validate_only`, once the selected matrix passes input validation
and its eigen decomposition is available:

| CSV | Contents |
|---|---|
| `GenomicPCA_All_PCs_Variance_Explained.csv` | Every PC's raw eigenvalue, raw percentage and cumulative percentage, separately labelled positive-eigenvalue-normalized percentages, and negative-eigenvalue flags. |
| `GenomicPCA_PC1_QC.csv` | PC1 eigenvalue/finite/length/nonzero checks, sign convention and multiplier, PASS/FAIL and reasons. Overall matrix policy is reported separately from PC1 validity. |
| `GenomicPCA_PC1_Protein_Contributions.csv` | Each retained trait/protein in manifest order, signed PC1 eigenvector coefficient and loading used for GWAMA, contribution percentage and rank (ties share rank). |

Raw percentage = `100 * eigenvalue / sum(all eigenvalues)`. If the total is
non-positive it is undefined (NA). If any eigenvalues are negative, raw percentages
are **not a conventional non-negative variance partition** and can be negative or
exceed 100%. The separate positive-normalized percentage is
`100 * max(eigenvalue, 0) / sum(max(eigenvalues, 0))`; it is a descriptive summary,
not a repair of the estimated matrix. Percentages describe the chosen correlation
or covariance matrix, not necessarily measured phenotypic variance.

Protein PC1 contribution = `100 * (PC1 eigenvector coefficient)^2`. Contributions
sum to 100% for a valid PC1 and are unchanged by reversing the entire PC1 sign.
They quantify relative squared loading magnitude, **not** a causal contribution
or each protein's final per-SNP GWAMA weight (which also depends on sample size).
The signed loading is provided separately to retain direction information.

Negative later PCs do not themselves fail PC1 QC. Existing matrix QC and explicit
`--negative_eigen_action error` still apply; PC1 PASS does not mean the full run
passed. Reports include FAIL if decomposition succeeds but PC1 is invalid;
upstream input errors may stop before these reports can be generated.

`Problematic_Trait_Combinations.csv` is no longer generated. Negative-eigenvalue
flags remain matrix diagnostics, not a protein-removal list. Historical files
from older runs are not deleted; use a fresh output directory. Existing weight,
eigenvalue and warning/removal audit files remain available.

## Package layout

```text
ldsc-gpca/
├── README.md
├── pyproject.toml
├── MANIFEST.in
└── src/ldsc_gpca/
    ├── cli.py              # Command dispatch
    ├── prepare.py          # VCF to GPCA and optional munging inputs
    ├── ldsc_cli.py         # LDSC options and orchestration
    ├── extraction.py       # VCF filtering and prevalence preparation
    ├── munging.py          # Summary-statistic munging
    ├── pairwise.py         # LDSC batches and retries
    ├── results.py          # Result validation and compilation
    ├── utils.py            # Shared helpers
    ├── gpca.py             # R launcher
    ├── genomicsem.py       # Native GenomicSEM GPCA dispatch
    ├── postprocess.py      # Automatic GWAMA combination and export
    └── r/
        ├── gpsca_gwama_python_ldsc.r # Thin Python-table entry point
        ├── gpsca_gwama_v2.r          # Thin native GenomicSEM entry point
        ├── load_modules.R          # Loads shared code plus one backend
        ├── shared/                 # Manifest, matrix checks, PCA, reports, GWAMA
        ├── python_ldsc/            # CLI, constants, reader, QC, covariance, reports, workflow
        └── genomicsem/             # CLI, native reader, QC, reports, workflow
```

Both backends use the same PCA, whole-vector sign orientation and GWAMA helpers.
GenomicSEM does not load the Python-LDSC reader or workflow. Backend-specific QC
and reports stay with their backend; shared numerical functions are defined once.
All R modules are included in the wheel. Keep the full `r/` directory together
when running R directly; neither thin entry point is a self-contained R script.
The local repository-root `gpsca_gwama_python_ldsc.r` is a compatibility launcher
for `src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r`, so existing local commands still work.

## Attribution

Method: [Anna Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
Upstream: [Bulik-Sullivan LDSC](https://github.com/bulik/ldsc).
Packaging follows [setuptools entry points](https://setuptools.pypa.io/en/latest/userguide/entry_point.html)
and [package data](https://setuptools.pypa.io/en/latest/userguide/datafiles.html).

## License

No redistribution license has been selected for this repository. Public availability
does not grant an open-source license. Third-party software retains its own licensing.
