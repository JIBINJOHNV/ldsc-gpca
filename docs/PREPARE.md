# Prepare GPCA and LDSC inputs

[README](../README.md) · [Input formats](INPUTS.md) · [Python LDSC](PYTHON_LDSC.md) · [GenomicSEM LDSC](GENOMICSEM_LDSC.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Choose what to prepare](#choose-what-to-prepare)
- [Files you need](#files-you-need)
- [Basic command](#basic-command)
- [Prepare shared raw LDSC inputs](#prepare-shared-raw-ldsc-inputs)
- [Produce munged files](#produce-munged-files)
- [What happens in order](#what-happens-in-order)
- [Choose IDs, layout and P-value floor](#choose-ids-layout-and-p-value-floor)
- [Outputs and next step](#outputs-and-next-step)
- [Legacy raw export](#also-write-raw-tables-for-later-munging)
- [All options](#all-options)

## When to use this command

Use `ldsc-gpca prepare` to turn GWAS summary-statistics VCFs into inputs for a
later analysis. You can make GPCA/GWAMA tables, shared raw LDSC tables, or
LDSC-ready munged files. Preparation accepts one or more traits.

**Raw munging inputs still need munging. Munged files are ready for the selected
LDSC backend. Neither is a completed LDSC result.** This command stops before
LDSC regression, PCA and GWAMA. No LD-score directory is needed here.

## Choose what to prepare

The default is GPCA-only preparation. `--mode ldsc` adds munging by default;
`--raw_only` stops before that step. `--mode both` also creates GPCA/GWAMA tables.

| What you need | Options after `ldsc-gpca prepare`¹ | Folders created |
| --- | --- | --- |
| GPCA/GWAMA inputs only | No mode option, or `--mode gpca` | `gpca_inputs/` |
| Shared raw LDSC tables only | `--mode ldsc --raw_only` | `munge_inputs/` |
| Raw LDSC tables and Python-munged files | `--mode ldsc` | `munge_inputs/`, `munged/` |
| Raw LDSC tables and GenomicSEM-munged files | `--mode ldsc --munge_backend genomicsem` | `munge_inputs/`, `munged/` |
| GPCA inputs and shared raw LDSC tables | `--mode both --raw_only` | `gpca_inputs/`, `munge_inputs/` |
| All outputs, with Python munging | `--mode both` | All three folders |
| All outputs, with GenomicSEM munging | `--mode both --munge_backend genomicsem` | All three folders |

¹ Also supply `--input`, `--outdir`, and `--hm3` for either LDSC mode.
Every mode writes QC reports. Use `--splitby_chr nosplit` for whole-genome GPCA
files unless your traits have valid records on every chromosome 1–22.

`--munge_backend` selects the program performing munging. The raw table format
is shared. This option belongs to standalone preparation. For `ldsc` or
`pipeline`, use `--ldsc_backend` instead; it selects both the matching munger and
regression implementation, with no separate munging override. Omit
`--munge_backend` with `--raw_only` or `--mode gpca`; irrelevant
munging options are rejected. **For later GenomicSEM LDSC reuse, select
`--munge_backend genomicsem`**, using the same `--hm3` allele reference for all
traits. Do not directly reuse Python-munged files in GenomicSEM: Python
HapMap matching can retain strand complements rather than rewriting all alleles
to the reference orientation. The GenomicSEM reuse check does not verify this
orientation or the preparation backend. See the [GenomicSEM reuse guide](GENOMICSEM_LDSC.md#reuse-munged-files).

## Files you need

Supply a CSV with `traitname,vcf_files`. Relative VCF paths resolve beside this
manifest. Each VCF must contain exactly one GWAS sample, not individual-level
genotypes. Plain VCF and gzip-compressed VCF are accepted.

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

Save this quantitative example as `/data/traits_vcf.csv`; replace the paths with
your files. All modes need variant coordinates, REF/ALT and FORMAT `AF,ES,SE,LP`.
LP means −log10(P); ES is the signed effect for ALT.

| Requested output | Additional fields and files |
| --- | --- |
| GPCA/GWAMA tables | FORMAT `NEF`, used for per-SNP GWAMA N. No INFO field or HapMap reference is required. |
| Shared raw LDSC tables | FORMAT `SI` for INFO; FORMAT `NEF` unless manifest N overrides it. A tab-separated `--hm3` file with a `SNP` header. Selected IDs must match this reference. |
| Munged LDSC files | The raw-input requirements above, plus `SNP,A1,A2` headers in `--hm3` and the chosen munging runtime. |

You need `bcftools` on PATH or `--bcftools /path/to/bcftools` for every mode.
Python munging uses the isolated LDSC environment configured by installation.
GenomicSEM munging needs Rscript and the pinned GenomicSEM package. GPCA-only
and raw-only preparation do not require either munging runtime.

### Sample size and binary traits

For new LDSC modes, a positive, finite manifest `N` replaces per-SNP NEF in
**LDSC tables only**. Without it, quantitative traits use NEF. In `--mode both`,
GPCA/GWAMA tables always keep NEF, so NEF remains required even with an N override.
The original VCF is never modified.

**Binary traits require explicit total N (cases + controls) and both prevalence
values in this shared preparation route.** The command does not infer total N
from NEF or case-count fields. Sample prevalence must be the case fraction
appropriate for that total N. Each prevalence must be strictly between 0 and 1.
For quantitative traits, leave both empty. For example:

```csv
traitname,vcf_files,N,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,50000,,
Trait_B,/data/Trait_B.vcf.gz,60000,,
Binary_C,/data/Binary_C.vcf.gz,40000,0.1,0.25
Binary_D,/data/Binary_D.vcf.gz,50000,0.05,0.4
```

These sample sizes and prevalences are illustrative. For binary meta-analysis
requiring an effective-N convention or varying per-SNP total N, supply suitable
raw/munged tables through the existing LDSC input routes instead. Merely renaming
an effective-N column to N does not make it total N.

## Basic command

The existing GPCA-only command keeps its behavior:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_gpca \
  --splitby_chr nosplit
```

It writes `gpca_inputs/Trait_A_GenomicPCA_inputs.tsv` and one file per other
trait, plus QC reports. Each file has nine tab-separated columns:

```text
SNPID CHR BP EA OA EAF N Z P
```

Pass the `gpca_inputs/` subdirectory to `--gpca_input_folder` when running GPCA
or pipeline, with the same `--splitby_chr` setting.

## Prepare shared raw LDSC inputs

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_raw \
  --mode ldsc \
  --raw_only \
  --hm3 /references/hm3_alleles.tsv
```

This writes `munge_inputs/{traitname}_munge_inputs.tsv` and a reusable
`Prepared_LDSC_Manifest.csv`. It does not create GPCA files or run munging.
`--splitby_chr` has no effect in LDSC-only mode: raw LDSC files are always
whole-autosome tables.

Both Python LDSC and GenomicSEM can read the shared, tab-separated schema:

| Column | Meaning / VCF source |
| --- | --- |
| `SNP` | Selected variant identifier; VCF ID by default. Must match HapMap. |
| `CHR` | Autosome 1–22, from CHROM. |
| `BP` | Positive integer position, from POS. |
| `A1` | Effect allele, ALT. |
| `A2` | Other allele, REF. |
| `EAF` | Effect-allele frequency, FORMAT/AF. |
| `BETA` | Signed effect for A1, FORMAT/ES. Binary effects must be on the log-odds scale. |
| `SE` | Effect standard error, FORMAT/SE. |
| `P` | 10^(-LP), with the requested floor. |
| `N` | Manifest N when supplied; otherwise quantitative FORMAT/NEF. |
| `INFO` | FORMAT/SI. |

Basic QC removes missing/invalid values and resolves duplicates. HapMap selection
matches IDs at this stage; it does not align alleles. Valid low INFO/MAF values
remain in the raw table for the selected munger to filter later. Missing or
invalid SI is removed and audited for LDSC; SI must be finite and in `[0,1]`.

Add `--mode both --raw_only --splitby_chr nosplit` instead of `--mode ldsc
--raw_only` to produce GPCA tables alongside these raw tables. Their QC is
separate: LDSC INFO validity and HapMap selection do not filter GPCA rows.

## Produce munged files

Omit `--raw_only` to run munging. Python is the default:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_python \
  --mode ldsc \
  --hm3 /references/hm3_alleles.tsv \
  --info_filter 0.9 \
  --maf_filter 0.01 \
  --n_cores 4
```

To make all three outputs using GenomicSEM munging:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_native \
  --mode both \
  --munge_backend genomicsem \
  --hm3 /references/hm3_alleles.tsv \
  --splitby_chr nosplit \
  --info_filter 0.9 \
  --maf_filter 0.01 \
  --n_cores 4
```

The new preparation modes default to INFO 0.9 and MAF 0.01 for either munger.
They do not apply Python `ldsc` VCF extraction's separate INFO/AF–EUR comparison
or optional MHC filtering. Those extraction settings are not options of this
command. Python munging uses BETA as the signed statistic here; its column
mappings are supplied internally.

The original [Python munging](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py)
and [GenomicSEM munging](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge_main.R)
implementations are unchanged. Python uses a strict MAF
`>` comparison and has its own low-N and ambiguous-SNP handling. GenomicSEM uses
`>=` and aligns output to reference alleles. Both reconstruct Z from P and effect
direction. Matching thresholds and a shared format therefore do not guarantee
identical SNP sets or estimates. Inspect each trait's munging log.

## What happens in order

1. Validate the chosen mode, input files, reference and any required runtime.
2. Extract the VCF once per trait, adding SI only for the new LDSC modes.
3. Apply basic QC and duplicate handling independently for requested GPCA and
   LDSC outputs. Preserve manifest trait order and audit affected original rows.
4. Write requested GPCA and raw LDSC tables after every trait succeeds.
5. For either LDSC mode, run the chosen munger unless `--raw_only` is set.
   GPCA-only mode ends after preparation.
6. Validate the munged output schema and usable N/Z rows; publish `munged/` only
   after all traits succeed. Save provenance and update the reusable manifest.

Failed workers receive one retry, for two attempts total. Failed extraction
publishes QC but no prepared tables. Failed munging leaves successfully prepared
raw/GPCA tables and logs available, and publishes no final `munged/` directory
when a trait fails. `Preparation_Settings.json` records `completed` only after
the directory and updated manifest have both been published. A later publication
failure records `failed` when possible; files already published are retained for
diagnosis and must not be treated as a completed preparation. Manifest and status
updates use atomic replacement, so an interrupted write cannot truncate them.
Existing named outputs are refused; use a fresh directory for a new run.

## Choose IDs, layout and P-value floor

`--gpca_id_source chr_pos_ref_alt` is the GPCA default. `--munge_id_source vcf_id`
is the LDSC default. Choose IDs that match across traits and the intended
reference. Neither option performs genome-build conversion or strand resolution.

Split GPCA output requires variants on all chromosomes 1–22 for every trait.
`--splitby_chr nosplit` writes one whole-autosome GPCA file per trait. LDSC tables
are always whole-autosome files.

Duplicates retain the largest valid LP, with the first source row breaking ties.
GPCA and LDSC ID choices are evaluated separately. `--p_min` defaults to `1e-300`.
It changes written P values but not the ES/SE-derived GWAMA Z. Munging derives Z
from P, so the floor can affect extreme munged Z. Adjusted rows are audited.
Values originate from bcftools numeric decoding; original textual precision is
not reconstructed.

## Outputs and next step

| Output | Meaning |
| --- | --- |
| `gpca_inputs/{traitname}_GenomicPCA_inputs.tsv` | Whole-genome GPCA/GWAMA input. Split mode inserts `_chr{CHR}` after the trait name. |
| `munge_inputs/{traitname}_munge_inputs.tsv` | Shared raw LDSC table with the eleven columns defined above; needs munging. |
| **`munged/{traitname}.sumstats.gz`** | Final LDSC-ready files. Core headers are `SNP,N,Z,A1,A2`; order can differ by backend. Python can include reference-only rows with missing statistics; later LDSC discards these. |
| `munged/{traitname}.prevalence.json` | Output checksum, sample-size convention, prevalence, backend and munging thresholds. Keep beside the munged file for Python reuse. |
| **`Prepared_LDSC_Manifest.csv`** | Trait order, `ref=yes`, raw/munged paths, N override and prevalences. Munged paths are filled only after successful munging. |
| `Preparation_Status.csv` | Per-trait `gpca_rows`, `munge_rows`, excluded non-autosomal count and success/error. |
| `Preparation_Settings.json` | Input paths, ID/layout/P rules and, when run, munging settings/status. |
| `GPCA_Input_QC_Summary.csv`, `GPCA_Input_QC_Issues.csv` | GPCA counts and original affected VCF records. Present when GPCA tables are requested. |
| `LDSC_Input_QC_Summary.csv`, `LDSC_Input_QC_Issues.csv` | Shared LDSC preparation counts, HapMap removals, N source/override and original affected records. Present in new LDSC modes. |
| `Preparation_Worker_Attempts.csv` | Extraction/preparation worker attempts: `Job,Attempt,Success,Error`. |
| `Preparation_Munging_Status.csv` | Munging success/error plus output, usable and missing row counts. Counts may be unavailable for failures. |
| `Preparation_Munging_Worker_Attempts.csv` | Munging attempts: `Job,Attempt,Success,Error`. |
| `munge_logs/` | Exact per-trait command JSON, console output and munging logs. Console output includes retries. |

QC summaries contain `traitname,input_rows,retained_rows,removed_rows,
p_adjusted_rows,excluded_non_autosomal_rows,success,error`. LDSC summaries also
contain `hapmap_removed_rows,N_source,N_override,p_min`. LDSC retained counts
include HapMap selection. QC issue files preserve original VCF fields and add
`traitname,QC_action,QC_reason`.

To run **Python LDSC** using the prepared munged files:

```bash
ldsc-gpca ldsc \
  --input /results/prepared_python/Prepared_LDSC_Manifest.csv \
  --ldsc_only \
  --munged_dir /results/prepared_python/munged \
  --ld_ref /references/eur_w_ld_chr \
  --outdir /results/python_ldsc
```

Python reuse verifies the checksum and N convention. Shared-preparation provenance
is reported; original Python VCF extraction filters are not reapplied.

To run **GenomicSEM LDSC**, reuse files produced with
`--mode ldsc --munge_backend genomicsem` or
`--mode both --munge_backend genomicsem`, using the same allele reference for
all traits. Use the GenomicSEM preparation directory, not the Python-munged output:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /results/prepared_native/Prepared_LDSC_Manifest.csv \
  --munged_dir /results/prepared_native/munged \
  --ld_ref /references/eur_w_ld_chr \
  --outdir /results/genomicsem_ldsc
```

Use appropriate LD references and the same N/prevalence convention. For shared
raw tables, their `sumstats_file` manifest can instead be passed to GenomicSEM LDSC
with `--hm3`; omit munged-reuse options so GenomicSEM munging runs there.

## Also write raw tables for later munging

The older `--write_munge_inputs --hm3 FILE` option remains available for existing
commands and automatic pipeline preparation. It keeps GPCA output plus the
legacy `munge_inputs/{traitname}_munge_inputs.txt` files, whose columns are:

```text
SNP CHR POS A1 A2 eaf_A1 beta se N p
```

That legacy export uses N=NEF and contains no INFO. It does not run munging.
For the new shared schema with INFO, use `--mode both --raw_only` instead.
Do not combine the legacy flag with `--mode ldsc`, `--mode both` or `--raw_only`.
This preserves existing pipeline behavior; the new standalone modes are explicit.

## All options

Required options have no default. Boolean switches are off unless supplied.

| Option | Default | Applies to / meaning |
| --- | --- | --- |
| `--input` | Required | CSV with `traitname,vcf_files`; LDSC metadata rules above. |
| `--outdir` | Required | Output directory; existing named outputs are protected. |
| `--mode` | `gpca` | `gpca`, `ldsc` or `both`. LDSC modes run munging unless raw-only. |
| `--raw_only` | Off | `ldsc`/`both`: stop before munging. |
| `--munge_backend` | `python` | `python` or `genomicsem`; select `genomicsem` for later GenomicSEM LDSC reuse. Select only when munging runs. |
| `--hm3` | Unset | Tab-separated HapMap reference; required for LDSC or legacy raw export. Munging also requires A1/A2. |
| `--splitby_chr` | `split` | GPCA only: `split` or `nosplit`; ignored in LDSC-only mode. |
| `--gpca_id_source` | `chr_pos_ref_alt` | GPCA IDs: `chr_pos_ref_alt` or `vcf_id`. |
| `--munge_id_source` | `vcf_id` | Raw LDSC IDs: `vcf_id` or `chr_pos_ref_alt`; must match HapMap. |
| `--p_min` | `1e-300` | P floor, finite and strictly between 0 and 1; affected rows audited. |
| `--info_filter` | `0.9` | Munging only; finite value in [0,1]. |
| `--maf_filter` | `0.01` | Munging only; finite value in [0,0.5]. Boundary comparison is backend-specific. |
| `--n_cores` | `4` | Positive worker count for preparation and munging; stages run in order. |
| `--bcftools` | `bcftools` | Local VCF query executable or path. |
| `--conda_executable` | `CONDA_EXE`, otherwise `conda` | Python munging only: Conda executable. |
| `--ldsc_env` | `ldsc-cbiit` unless a prefix is configured | Python munging: named environment; mutually exclusive with explicit prefix. |
| `--ldsc_env_prefix` | `LDSC_GPCA_LDSC_PREFIX`, otherwise unset | Python munging: environment directory. An explicit environment name bypasses the saved prefix. |
| `--rscript` | `Rscript` | GenomicSEM munging only: executable with the pinned package installed. |
| `--write_munge_inputs` | Off | Legacy GPCA-plus-raw export with its original schema. Requires --hm3; incompatible with new LDSC modes. |
| `--help` | Off | Show help, defaults and choices without analysis. |
