# Prepare per-trait GWAMA inputs

[README](../README.md) · [Input formats](INPUTS.md) · [Pipeline](PIPELINE.md)

## In this guide

- [When to use this command](#when-to-use-this-command)
- [Files you need](#files-you-need)
- [Basic command](#basic-command)
- [What happens in order](#what-happens-in-order)
- [Choose IDs, layout and P-value floor](#choose-ids-layout-and-p-value-floor)
- [Also write raw tables for later munging](#also-write-raw-tables-for-later-munging)
- [Outputs and next step](#outputs-and-next-step)
- [All options](#all-options)

## When to use this command

Use `ldsc-gpca prepare` when you have GWAS VCFs and need the per-SNP tables for
GWAMA. It creates those tables and QC reports. It does not run LDSC, PCA or
GWAMA. Preparation alone accepts one or more traits.

## Files you need

Supply a CSV with `traitname,vcf_files`. The four-trait
[VCF manifest](INPUTS.md#vcf-manifest-for-pipeline) works without changes.
Relative VCF paths resolve beside that manifest. Each VCF needs one GWAS sample
and FORMAT `AF,ES,SE,LP,NEF`, plus coordinates and alleles. See
[the field definitions](INPUTS.md#vcf-fields).

You need `bcftools` on PATH or `--bcftools /path/to/bcftools`.
You do not need LD scores or HapMap for GWAMA-only preparation.

## Basic command

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared \
  --splitby_chr nosplit
```

This writes one file per trait, for example
`/results/prepared/gpca_inputs/Trait_A_GenomicPCA_inputs.tsv`, containing:

```text
SNPID CHR BP EA OA EAF N Z P
```

The files are tab-separated. Pass the **`gpca_inputs` subdirectory** to
`--gpca_input_folder` in a subsequent GPCA or pipeline command, using the same
`--splitby_chr` setting.

## What happens in order

1. Read the manifest and extract the required VCF fields.
2. Keep valid autosomal records: finite values, positive N and SE, positive
   integer positions, EAF in `[0,1]`, nonnegative LP and valid unequal alleles.
3. Resolve duplicates, keeping the record with the largest LP; ties keep the
   first source record. Repeated chromosome–position–REF–ALT combinations
   are checked, as are repeated IDs when VCF IDs are selected.
4. Set EA=ALT, OA=REF, EAF=AF, N=NEF, Z=ES/SE and P from `10^(-LP)`, applying
   the requested P floor. Sort and write the selected layout.
5. Publish the prepared outputs after all traits succeed, retaining QC/audit
   information. Existing named outputs are refused.

There is no INFO threshold or constant-N override in preparation. The filters
used by an LDSC command do not automatically apply to these GWAMA inputs.

## Choose IDs, layout and P-value floor

This example preserves VCF IDs, floors very small P values at `1e-250` and
prepares chromosome-split files with eight workers:

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_split \
  --splitby_chr split \
  --gpca_id_source vcf_id \
  --p_min 1e-250 \
  --n_cores 8
```

Use `vcf_id` only when IDs are present and compatible across traits. The default
`chr_pos_ref_alt` builds IDs from position and alleles. Neither choice aligns
strand flips or converts genome builds. Split mode requires retained variants
on **every chromosome 1–22 for every trait**; use `nosplit` otherwise.

The P floor changes the written P value, not the ES/SE-derived GWAMA Z. If the
optional raw tables below are later munged using P, the floor can affect the
munged Z.

## Also write raw tables for later munging

```bash
ldsc-gpca prepare \
  --input /data/traits_vcf.csv \
  --outdir /results/prepared_with_raw_tables \
  --splitby_chr nosplit \
  --write_munge_inputs \
  --hm3 /references/hm3_snp_list.tsv \
  --munge_id_source vcf_id
```

Here HapMap is a tab-separated file with a `SNP` header. Its IDs must match the
selected munging IDs. The additional space-separated files are
`munge_inputs/{traitname}_munge_inputs.txt`, with:

```text
SNP CHR POS A1 A2 eaf_A1 beta se N p
```

These files are **not munged LDSC files**. They use N=NEF, contain no INFO column
and have been selected by HapMap ID without allele alignment. To use them with
native LDSC, create a manifest with `sumstats_file` paths and run
[`genomicsem ldsc`](GENOMICSEM_LDSC.md). The quantitative VCF route in
[`pipeline`](PIPELINE.md#native-genomicsem-from-quantitative-vcfs) handles this
handoff automatically. Do not assume the generated N convention is appropriate
for a binary trait.

## Outputs and next step

Check preparation QC reports and the expected per-trait files before proceeding.
If a trait fails, inspect its preparation audit and fix the data; existing
successful named outputs are not silently overwritten. A fresh output directory
is the simplest way to rerun with different choices.

Use `gpca_inputs/` with [Python GPCA](PYTHON_GPCA.md) or
[native GPCA](GENOMICSEM_GPCA.md) once LDSC results exist. See the
[preparation output reference](REFERENCE.md#preparation-outputs-and-filtering)
for the audit-file catalog.

## All options

Required options have no default. “Off” means omit the flag; include it alone to enable. The tables cover this command, including wrapper preparation/export options.

### Files and input modes

| Option | Default | What it changes |
| --- | --- | --- |
| `--input` | Required | Manifest CSV. Use the schema for your selected input mode; preserve trait order. |
| `--outdir` | Required | Output directory. Use a fresh directory for a separate analysis. |
| `--hm3` | Unset | Tab-separated file with SNP header; needed only for optional --write_munge_inputs. Preparation selects IDs without allele alignment. |

### VCF preparation and file layout

| Option | Default | What it changes |
| --- | --- | --- |
| `--splitby_chr` | `split` | `split`: all chromosomes 1–22 per trait. `nosplit`: one autosomal table per trait. Preparation and GWAMA must use matching layouts. |
| `--gpca_id_source` | `chr_pos_ref_alt` | `chr_pos_ref_alt` or `vcf_id` for prepared GWAMA SNPIDs. Choose compatible IDs across traits. |
| `--p_min` | `1e-300` | Preparation P floor, finite and strictly between 0 and 1. Does not change the ES/SE-derived GWAMA Z. |
| `--write_munge_inputs` | Off | Also prepare HapMap-selected raw tables, using NEF and no INFO column. Does not itself run munging. |
| `--munge_id_source` | `vcf_id` | `vcf_id` or `chr_pos_ref_alt` for optional raw-table SNP IDs; must match HapMap IDs. |

### Execution

| Option | Default | What it changes |
| --- | --- | --- |
| `--n_cores` | `4` | Positive integer; parallel VCF preparation workers. |
| `--bcftools` | `bcftools` on PATH | Executable name/path for VCF extraction or preparation; unnecessary when neither runs. |

### Help

| Option | Default | What it changes |
| --- | --- | --- |
| `--help` | Off | Show command usage, defaults and choices; exit without analysis. |
