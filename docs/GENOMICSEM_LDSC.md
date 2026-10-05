# Run GenomicSEM LDSC

[README](../README.md) · [Inputs](INPUTS.md) · [Pipeline](PIPELINE.md) · [Next: GPCA/GWAMA](GENOMICSEM_GPCA.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Choose the input mode explicitly](#choose-the-input-mode-explicitly)
- [Files and references you need](#files-and-references-you-need)
- [Start from GWAS VCFs](#start-from-gwas-vcfs)
- [Basic command from raw tables](#basic-command-from-raw-tables)
- [What happens in order](#what-happens-in-order)
- [Filters, sample size and prevalence](#filters-sample-size-and-prevalence)
- [Reuse munged files](#reuse-munged-files)
- [Results and next step](#results-and-next-step)
- [All options](#all-options)

## When to use this command

`ldsc-gpca ldsc --ldsc_backend genomicsem` runs GenomicSEM's R LDSC workflow and produces
the final `LDSCoutput` object for GenomicSEM GPCA. It accepts raw GWAS text tables,
GWAS VCFs with `--vcf_input`, or existing munged files. VCF mode converts the
records to raw tables before calling the unchanged GenomicSEM munging/LDSC workflow.
This command stops after LDSC; it does not prepare GWAMA tables or run PCA.

The existing `ldsc-gpca genomicsem ldsc` command remains available and uses the
same workflow. With the unified `ldsc` command, explicitly select
`--ldsc_backend genomicsem`; omission defaults to Python. No separate
`--munge_backend` is accepted: raw inputs use GenomicSEM munging automatically.
Use `ldsc-gpca ldsc --ldsc_backend genomicsem --help` for this backend's options.

## Choose the input mode explicitly

The command uses one mode for all traits:

1. **Raw tables, the default:** put paths in manifest `sumstats_file`. The
   command runs `GenomicSEM::munge()` before LDSC.
2. **Munged directory:** set `--munged_dir`; the directory must contain exactly
   one `{traitname}.sumstats` or `{traitname}.sumstats.gz` per trait.
3. **Explicit munged paths:** set `--munged_input` and supply manifest
   `munged_file` paths.
4. **GWAS VCFs:** set `--vcf_input` and supply manifest `vcf_files` paths.
   The command extracts raw tables with INFO and allele frequency, then munges
   them and runs LDSC.

`--vcf_input`, `--munged_dir` and `--munged_input` are mutually exclusive. The
tool does not switch modes based on filenames. Put VCFs in `vcf_files` and
select `--vcf_input`; put raw tables in `sumstats_file` for the default mode.

For either munged reuse mode, produce the files with
`prepare --mode ldsc --munge_backend genomicsem` (or `--mode both`), using the
same `--hm3` allele reference for all traits. Do not directly reuse Python-munged
files. Raw-table and VCF modes run GenomicSEM munging themselves; a separate prepare
step is not required. See [reuse requirements and commands](#reuse-munged-files).

## Files and references you need

Use at least two traits. All manifests need
`traitname,population_prevalence,sample_prevalence` and the path column required
by the selected mode. Relative paths resolve beside the manifest. For binary
traits supply both prevalences between 0 and 1; for quantitative traits leave
both empty/NA. GenomicSEM LDSC does not infer a missing prevalence.

The [four-trait raw manifest](INPUTS.md#raw-table-manifest-for-genomicsem) and
[table schemas](INPUTS.md#raw-and-munged-tables) show the starting files. Raw
input needs SNP IDs, alleles, P, signed effects and appropriate sample sizes.
A positive manifest `N` overrides raw-file N during munging. Recognized INFO
and frequency columns are needed for their respective filters.

You need the installed R/GenomicSEM environment, LD-score reference files and,
for raw munging, a HapMap allele reference with `SNP,A1,A2`. Reuse requires
tab-separated `SNP,A1,A2,N,Z` files from GenomicSEM preparation and does not accept `--hm3`.
See [references](INPUTS.md#reference-files).

## Start from GWAS VCFs

Use one GWAS summary-statistics VCF per trait, with exactly one GWAS sample in
each file. Save this quantitative example as `/data/traits_genomicsem_vcf.csv`:

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/traits_genomicsem_vcf.csv \
  --vcf_input \
  --outdir /results/genomicsem_ldsc_vcf \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist \
  --info_filter 0.9 \
  --maf_filter 0.01 \
  --n_cores 4
```

The thresholds shown are GenomicSEM defaults. `--bcftools /path/to/bcftools` can
select a local executable. Use a fresh or empty output directory.

| VCF value | Raw table column / use |
| --- | --- |
| Variant ID | `SNP`; IDs must match the HapMap reference. |
| CHROM, POS | `CHR,POS`; retain valid autosomal positions. |
| ALT, REF | `A1,A2`: effect and other allele. |
| FORMAT/ES, FORMAT/SE | `beta,se`; ES supplies direction for GenomicSEM munging. Both must be finite and SE positive. |
| FORMAT/LP | `p = 10^(-LP)`, with a default floor of `1e-300`. GenomicSEM munging calculates Z from this P and the effect direction. |
| FORMAT/AF | `MAF = min(AF, 1−AF)`. GenomicSEM `--maf_filter` applies during munging. |
| FORMAT/SI | `INFO`. GenomicSEM `--info_filter` applies during munging. |
| FORMAT/NEF | Per-SNP `N` for quantitative traits, unless manifest `N` overrides it. |

The adapter uses the existing preparation QC: remove invalid records, retain
autosomes, and resolve duplicate IDs/coordinates by keeping the largest valid
LP (first source row on a tie). Missing or invalid SI is removed and audited;
finite SI must be in `[0,1]`. Low but valid INFO/MAF values reach GenomicSEM munging,
which applies your thresholds. Extraction does not apply Python's INFO/AF–EUR
comparison, MHC filter or Python-specific munging settings. These routes can
therefore retain different SNPs even when starting from the same VCFs.

**Sample size:** quantitative traits use NEF by default. A positive `N` column
in the manifest replaces it with that constant and makes NEF unnecessary for
this LDSC conversion. Binary VCFs require an explicit, appropriate manifest N
and both prevalences. Choose N and sample prevalence together for your study's
GenomicSEM convention; the command does not infer that convention from NEF or
NC/NCO. For binary data needing varying per-SNP N, supply suitable raw tables or
munged files instead. A mixture of quantitative and binary rows is supported
when each row satisfies these rules.

`--p_min` changes the VCF conversion floor and records adjusted rows. Because
GenomicSEM munging derives Z from P, flooring can affect extreme Z values. Values
come from bcftools' numeric decoding; the adapter does not reconstruct the
original VCF's textual precision.

After successful conversion and LDSC, the output includes:

| File or folder | Meaning |
| --- | --- |
| `vcf_input/{traitname}_munge_inputs.tsv` | Tab-separated raw inputs: `SNP,CHR,POS,A1,A2,MAF,beta,se,N,p,INFO`, in that order. These are not munged files. |
| `GenomicSEM_VCF_QC_Summary.csv` | Per-trait conversion counts, N source/override, P floor and success/error. Counts describe conversion before GenomicSEM INFO/MAF filtering. |
| `GenomicSEM_VCF_QC_Issues.csv` | Original VCF record fields plus `traitname,QC_action,QC_reason` for removed or P-adjusted rows. |
| `GenomicSEM_VCF_Worker_Attempts.csv` | `Job,Attempt,Success,Error`; each failed conversion gets one retry, for two total attempts. |
| `Resolved_Manifest.csv` | Resolved raw-table paths, trait order, N and prevalence passed to GenomicSEM munging. |
| `munge_output/` | GenomicSEM munged files and munging logs. |
| **`genomicsem_LDSC.RData`** | Final GenomicSEM LDSC result for GPCA, as in the other input modes. |
| `Selected_Traits.csv`, `GenomicSEM_LDSC_Trait_QC.csv` | Retained traits and GenomicSEM heritability QC. |

The conversion summary columns are `traitname,vcf_files,input_rows,retained_rows,
removed_rows,p_adjusted_rows,excluded_non_autosomal_rows,N_source,N_override,
p_min,INFO_source,success,error`. A failure before conversion can leave counts
unavailable. Failed conversion stops before R starts; converted tables are
published only after every trait succeeds. Original VCFs are unchanged.

## Basic command from raw tables

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/traits_genomicsem_raw.csv \
  --outdir /results/genomicsem_ldsc \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist
```

Use a fresh or empty output directory. The final GPCA input will be
`/results/genomicsem_ldsc/genomicsem_LDSC.RData`.

## What happens in order

1. Validate the manifest and selected input mode.
2. In VCF mode, extract and check the raw tables first. Munge raw tables, or
   use the explicitly selected existing munged files.
3. Run the first LDSC pass and save its raw result.
4. Check raw heritabilities. By default, audit and drop traits with nonfinite
   or nonpositive raw h2; `--invalid_h2_action error` stops instead.
5. Run the standardized LDSC pass for retained traits, check the matrices and
   preserve their trait order.
6. Save the final `LDSCoutput` object, retained-trait manifest and QC reports.

At least two valid traits must remain. The GenomicSEM workflow estimates the actual
`S,V,I,S_Stand,V_Stand` matrices; it does not reconstruct a GenomicSEM object from a
Python CSV.

## Filters, sample size and prevalence

This example uses stricter raw-munging thresholds, a fixed chi-square cutoff
and a policy that stops on invalid first-pass heritabilities:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/traits_genomicsem_raw.csv \
  --outdir /results/genomicsem_ldsc_filtered \
  --ld_ref /references/eur_w_ld_chr \
  --ld_weights /references/weights_hm3_no_hla \
  --hm3 /references/w_hm3.snplist \
  --info_filter 0.95 \
  --maf_filter 0.02 \
  --chisq_max 80 \
  --invalid_h2_action error \
  --n_cores 4
```

The default GenomicSEM munging thresholds are INFO `0.9` and MAF `0.01`, inclusive.
Filtering depends on recognized fields being present; a table with no INFO
cannot be quality-filtered by INFO. Reusing munged files skips these filters.
These settings do not filter the separate GWAMA input tables.

When `--chisq_max` is omitted, GenomicSEM uses its automatic
`max(80, 0.001 × maximum N)` rule in LDSC after relevant matching. To override,
supply a positive finite number. **Do not write `--chisq_max auto` here**: this
CLI does not accept that string. Python's managed option has different default
behavior and preprocessing; selecting an analogous cutoff does not guarantee
identical estimates.

To set a constant N during raw munging, add an `N` column to the manifest,
for example:

```csv
traitname,sumstats_file,N,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A_raw.tsv,50000,,
Trait_B,/data/Trait_B_raw.tsv,60000,,
Trait_C,/data/Trait_C_raw.tsv,70000,,
Trait_D,/data/Trait_D_raw.tsv,80000,,
```

These numbers are illustrative. This replaces per-variant N in each raw file;
it is not just missing-value filling. Omit manifest N to use appropriate
per-variant N already present. Manifest N does not replace N in reused munged
files. See [sample-size and prevalence rules](INPUTS.md#sample-size-and-prevalence).

## Reuse munged files

**Reused munged files must come from `prepare --munge_backend genomicsem`**
with `--mode ldsc` or `--mode both`, using the same `--hm3` allele reference for
all traits. Do not directly reuse Python-munged files, even when they used the
same HapMap reference. Python munging accepts strand complements and can retain
their coding; GenomicSEM LDSC assumes compatible allele orientation and can reverse
their cross-trait products incorrectly. The upstream munging and regression
implementations are unchanged. See [Python munging](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py)
and [GenomicSEM alignment](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R).

Prepare GenomicSEM inputs and pass the generated manifest and `munged/` directory:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_native \
  --mode ldsc \
  --munge_backend genomicsem \
  --hm3 /references/hm3_alleles.tsv

ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /results/prepared_native/Prepared_LDSC_Manifest.csv \
  --outdir /results/genomicsem_ldsc_reused \
  --ld_ref /references/eur_w_ld_chr \
  --munged_dir /results/prepared_native/munged
```

Alternatively, remunge shared raw tables by using the manifest's `sumstats_file`
paths in default raw mode with `--hm3`. If assembling a reuse manifest separately,
use the [GenomicSEM reuse schema](INPUTS.md#reusing-munged-files).

For files stored in different locations, save
`/data/traits_genomicsem_munged_paths.csv`:

```csv
traitname,munged_file,population_prevalence,sample_prevalence
Trait_A,/data/genomicsem_munged/Trait_A.sumstats.gz,,
Trait_B,/data/genomicsem_munged/Trait_B.sumstats.gz,,
Trait_C,/data/genomicsem_munged/Trait_C.sumstats.gz,,
Trait_D,/data/genomicsem_munged/Trait_D.sumstats.gz,,
```

Then run:

```bash
ldsc-gpca ldsc --ldsc_backend genomicsem \
  --input /data/traits_genomicsem_munged_paths.csv \
  --outdir /results/genomicsem_ldsc_explicit_paths \
  --ld_ref /references/eur_w_ld_chr \
  --munged_input
```

The current reuse check verifies the required headers and at least one data
row. It does not verify the preparation backend, allele orientation, sidecar
checksums or N/prevalence provenance. The preparation requirement above is a
documented input requirement, not an automatic compatibility check.

## Results and next step

Use **`genomicsem_LDSC.RData`** as the input to
[`genomicsem gpca`](GENOMICSEM_GPCA.md). It contains `LDSCoutput` with
`S,V,I,S_Stand,V_Stand`. Use `Selected_Traits.csv` to carry forward the retained
traits in order.

Update scripts that refer to the previous output filename. Existing RData files
remain readable through an explicit `--ldsc_results` path when they contain the
required `LDSCoutput` object and matrices; filenames do not determine compatibility.

`genomicsem_LDSC_raw.RData` is the first-pass audit result containing
`LDSCoutput_raw`; it is **not** the final GPCA input. Keep it, QC events, logs,
resolved manifest and session information for reproducibility. See the
[output catalog](REFERENCE.md#genomicsem-ldsc-outputs).

Munging workers retry each failed trait once; an exhausted failure stops before
LDSC. `--n_cores` controls munging, not GenomicSEM LDSC regression. There is no GenomicSEM
whole-run `--restart` flag. After a failure, correct the cause and use a fresh
output directory; explicitly reuse suitable munged files when appropriate.

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--ldsc_backend` | `python` on `ldsc` | Set `genomicsem` for this workflow. The existing `genomicsem ldsc` command implies GenomicSEM and does not take this selector. |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Fresh or empty GenomicSEM LDSC output directory. |
| `--ld_ref` | Required | Directory of chromosome LD scores and M reference files. |
| `--ld_weights` | Use `--ld_ref` | Optional separate directory of regression-weight LD scores. |
| `--hm3` | Unset | Whitespace SNP/A1/A2 allele reference, required for raw munging; omit for munged reuse. |
| `--vcf_input` | Off | Read `vcf_files` from the manifest, convert VCFs to raw tables and run GenomicSEM munging. Mutually exclusive with both munged input modes. |
| `--munged_dir` | Unset | Reuse files from `prepare --munge_backend genomicsem`: exactly one {trait}.sumstats or .sumstats.gz per trait. Mutually exclusive with `--munged_input` and `--vcf_input`. |
| `--munged_input` | Off | Reuse manifest `munged_file` paths from `prepare --munge_backend genomicsem`; mutually exclusive with `--munged_dir` and `--vcf_input`. |

### LDSC filtering and estimation

| Option | Default | What it changes |
| --- | --- | --- |
| `--info_filter` | `0.9` | GenomicSEM raw munging: inclusive INFO threshold from 0 to 1, when a recognized INFO column exists. Inactive for munged reuse. |
| `--maf_filter` | `0.01` | GenomicSEM raw munging: inclusive MAF threshold from 0 to 0.5, when recognized frequency data exist. Inactive for munged reuse. |
| `--chisq_max` | GenomicSEM automatic rule | Positive finite number to override GenomicSEM automatic max(80, 0.001 × maximum matched N). Omit for automatic; the string auto is not accepted. |
| `--chromosomes` | `22` | GenomicSEM LDSC: use reference chromosomes 1 through this integer (1–22). Does not change GWAMA split-file requirements. |
| `--n_blocks` | `200` | GenomicSEM LDSC: requested jackknife blocks, integer ≥2. Pinned GenomicSEM overrides the count for more than 18 traits; inspect its log. |

### VCF conversion

| Option | Default | What it changes |
| --- | --- | --- |
| `--bcftools` | `bcftools` on PATH | VCF query executable; used only with `--vcf_input`. |
| `--p_min` | `1e-300` | P floor for VCF LP conversion, finite and strictly between 0 and 1. Used only with `--vcf_input`; adjusted rows are audited. |

### Failures and reuse

| Option | Default | What it changes |
| --- | --- | --- |
| `--invalid_h2_action` | `drop` | GenomicSEM first-pass h2 QC: `drop` audits/removes nonfinite or nonpositive h2 traits; `error` stops. At least two must remain. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--n_cores` | `1` | Positive integer; VCF extraction and GenomicSEM munging workers. Regression is sequential; Windows munging is sequential. |
| `--rscript` | `Rscript` | GenomicSEM LDSC R executable/path. In pipeline this affects LDSC only; GPCA still needs Rscript on PATH. |

### Help

| Option | Default | What it changes |
| --- | --- | --- |
| `--help` | Off | Show command usage, defaults and choices; exit without analysis. |
