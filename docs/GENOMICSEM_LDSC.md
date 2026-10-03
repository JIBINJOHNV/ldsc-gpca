# Run native GenomicSEM LDSC

[README](../README.md) · [Inputs](INPUTS.md) · [Pipeline](PIPELINE.md) · [Next: GPCA/GWAMA](GENOMICSEM_GPCA.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Choose the input mode explicitly](#choose-the-input-mode-explicitly)
- [Files and references you need](#files-and-references-you-need)
- [Basic command from raw tables](#basic-command-from-raw-tables)
- [What happens in order](#what-happens-in-order)
- [Filters, sample size and prevalence](#filters-sample-size-and-prevalence)
- [Reuse munged files](#reuse-munged-files)
- [Results and next step](#results-and-next-step)
- [All options](#all-options)

## When to use this command

`ldsc-gpca genomicsem ldsc` runs GenomicSEM's native R LDSC workflow and produces
the final `LDSCoutput` object for native GPCA. It accepts raw GWAS text tables or
explicitly selected existing munged files. It does not directly read VCFs,
prepare GWAMA tables or run PCA.

## Choose the input mode explicitly

The command uses one mode for all traits:

1. **Raw tables, the default:** put paths in manifest `sumstats_file`. The
   command runs `GenomicSEM::munge()` before LDSC.
2. **Munged directory:** set `--munged_dir`; the directory must contain exactly
   one `{traitname}.sumstats` or `{traitname}.sumstats.gz` per trait.
3. **Explicit munged paths:** set `--munged_input` and supply manifest
   `munged_file` paths.

The two reuse options are mutually exclusive. The tool does not identify raw
versus munged data from filenames or automatically switch modes. Do not put a
munged file under `sumstats_file` and expect munging to be skipped.

## Files and references you need

Use at least two traits. All manifests need
`traitname,population_prevalence,sample_prevalence` and the path column required
by the selected mode. Relative paths resolve beside the manifest. For binary
traits supply both prevalences between 0 and 1; for quantitative traits leave
both empty/NA. Native LDSC does not infer a missing prevalence.

The [four-trait raw manifest](INPUTS.md#raw-table-manifest-for-genomicsem) and
[table schemas](INPUTS.md#raw-and-munged-tables) show the starting files. Raw
input needs SNP IDs, alleles, P, signed effects and appropriate sample sizes.
A positive manifest `N` overrides raw-file N during munging. Recognized INFO
and frequency columns are needed for their respective filters.

You need the installed R/GenomicSEM environment, LD-score reference files and,
for raw munging, a HapMap allele reference with `SNP,A1,A2`. Reuse requires
compatible tab-separated `SNP,A1,A2,N,Z` files and does not accept `--hm3`.
See [references](INPUTS.md#reference-files).

## Basic command from raw tables

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/traits_native_raw.csv \
  --outdir /results/native_ldsc \
  --ld_ref /references/eur_w_ld_chr \
  --hm3 /references/w_hm3.snplist
```

Use a fresh or empty output directory. The final GPCA input will be
`/results/native_ldsc/genomicPCA_LDSC.RData`.

## What happens in order

1. Validate the manifest and selected input mode.
2. Munge raw tables, or use the explicitly selected munged files.
3. Run the first LDSC pass and save its raw result.
4. Check raw heritabilities. By default, audit and drop traits with nonfinite
   or nonpositive raw h2; `--invalid_h2_action error` stops instead.
5. Run the standardized LDSC pass for retained traits, check the matrices and
   preserve their trait order.
6. Save the final `LDSCoutput` object, retained-trait manifest and QC reports.

At least two valid traits must remain. The native workflow estimates the actual
`S,V,I,S_Stand,V_Stand` matrices; it does not reconstruct a native object from a
Python CSV.

## Filters, sample size and prevalence

This example uses stricter raw-munging thresholds, a fixed chi-square cutoff
and a policy that stops on invalid first-pass heritabilities:

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/traits_native_raw.csv \
  --outdir /results/native_ldsc_filtered \
  --ld_ref /references/eur_w_ld_chr \
  --ld_weights /references/weights_hm3_no_hla \
  --hm3 /references/w_hm3.snplist \
  --info_filter 0.95 \
  --maf_filter 0.02 \
  --chisq_max 80 \
  --invalid_h2_action error \
  --n_cores 4
```

The default native munging thresholds are INFO `0.9` and MAF `0.01`, inclusive.
Filtering depends on recognized fields being present; a table with no INFO
cannot be quality-filtered by INFO. Reusing munged files skips these filters.
These settings do not filter the separate GWAMA input tables.

When `--chisq_max` is omitted, native GenomicSEM uses its automatic
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

Use the [native reuse manifest](INPUTS.md#reusing-munged-files), with prevalence
columns and no raw-file paths needed:

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/traits_native_munged.csv \
  --outdir /results/native_ldsc_reused \
  --ld_ref /references/eur_w_ld_chr \
  --munged_dir /data/native_munged
```

For files stored in different locations, save
`/data/traits_native_munged_paths.csv`:

```csv
traitname,munged_file,population_prevalence,sample_prevalence
Trait_A,/data/native_munged/Trait_A.sumstats.gz,,
Trait_B,/data/native_munged/Trait_B.sumstats.gz,,
Trait_C,/data/native_munged/Trait_C.sumstats.gz,,
Trait_D,/data/native_munged/Trait_D.sumstats.gz,,
```

Then run:

```bash
ldsc-gpca genomicsem ldsc \
  --input /data/traits_native_munged_paths.csv \
  --outdir /results/native_ldsc_explicit_paths \
  --ld_ref /references/eur_w_ld_chr \
  --munged_input
```

The initial structural check verifies required headers and early records; it
does not certify the inputs' previous QC, allele handling or N conventions.

## Results and next step

Use **`genomicPCA_LDSC.RData`** as the input to
[`genomicsem gpca`](GENOMICSEM_GPCA.md). It contains `LDSCoutput` with
`S,V,I,S_Stand,V_Stand`. Use `Selected_Traits.csv` to carry forward the retained
traits in order.

`genomicPCA_LDSC_raw.RData` is the first-pass audit result containing
`LDSCoutput_raw`; it is **not** the final GPCA input. Keep it, QC events, logs,
resolved manifest and session information for reproducibility. See the
[output catalog](REFERENCE.md#genomicsem-ldsc-outputs).

Munging workers retry each failed trait once; an exhausted failure stops before
LDSC. `--n_cores` controls munging, not native LDSC regression. There is no native
whole-run `--restart` flag. After a failure, correct the cause and use a fresh
output directory; explicitly reuse suitable munged files when appropriate.

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Fresh or empty native LDSC output directory. |
| `--ld_ref` | Required | Directory of chromosome LD scores and M reference files. |
| `--ld_weights` | Use `--ld_ref` | Optional separate directory of regression-weight LD scores. |
| `--hm3` | Unset | Whitespace SNP/A1/A2 allele reference, required for raw munging; omit for munged reuse. |
| `--munged_dir` | Unset | Select native munged reuse: exactly one {trait}.sumstats or .sumstats.gz per trait. Mutually exclusive with --munged_input. |
| `--munged_input` | Off | Native reuse from manifest `munged_file` paths; mutually exclusive with `--munged_dir`. |

### LDSC filtering and estimation

| Option | Default | What it changes |
| --- | --- | --- |
| `--info_filter` | `0.9` | Native raw munging: inclusive INFO threshold from 0 to 1, when a recognized INFO column exists. Inactive for munged reuse. |
| `--maf_filter` | `0.01` | Native raw munging: inclusive MAF threshold from 0 to 0.5, when recognized frequency data exist. Inactive for munged reuse. |
| `--chisq_max` | Native automatic rule | Positive finite number to override native automatic max(80, 0.001 × maximum matched N). Omit for automatic; the string auto is not accepted. |
| `--chromosomes` | `22` | Native LDSC: use reference chromosomes 1 through this integer (1–22). Does not change GWAMA split-file requirements. |
| `--n_blocks` | `200` | Native LDSC: requested jackknife blocks, integer ≥2. Pinned GenomicSEM overrides the count for more than 18 traits; inspect its log. |

### Failures and reuse

| Option | Default | What it changes |
| --- | --- | --- |
| `--invalid_h2_action` | `drop` | Native first-pass h2 QC: `drop` audits/removes nonfinite or nonpositive h2 traits; `error` stops. At least two must remain. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--n_cores` | `1` | Positive integer; native munging workers only. Regression is sequential; Windows munging is sequential. |
| `--rscript` | `Rscript` | Native LDSC R executable/path. In pipeline this affects LDSC only; GPCA still needs Rscript on PATH. |

### Help

| Option | Default | What it changes |
| --- | --- | --- |
| `--help` | Off | Show command usage, defaults and choices; exit without analysis. |
