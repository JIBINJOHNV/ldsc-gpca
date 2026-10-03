# Input files, references and sample sizes

[README](../README.md) · [Pipeline](PIPELINE.md) · [Full reference](REFERENCE.md)

## In this guide

- [Know which file you have](#know-which-file-you-have)
- [Manifests](#manifests)
- [VCF fields](#vcf-fields)
- [Raw and munged tables](#raw-and-munged-tables)
- [GWAMA tables](#gwama-tables)
- [Reference files](#reference-files)
- [Sample size and prevalence](#sample-size-and-prevalence)
- [Export metadata](#export-metadata)

## Know which file you have

A **manifest** is a small CSV with one row per trait. It names the traits and,
when needed, points to their data files. A **summary-statistics table** has one
row per variant. A **munged file** is a filtered per-trait table ready for LDSC.
**Completed LDSC results** describe trait heritabilities, correlations and
intercepts. **GWAMA input tables** contain the variant statistics to combine.

Use at least two traits for analysis. Preparation alone can process one. Trait
names must be unique and nonempty; use simple names without spaces or path
characters, such as `Trait_A`. Manifest row order defines the analysis order.
Extra annotation columns do not change that order.

Examples use four traits, input files under `/data`, LD references under
`/references` and outputs under `/results`. Replace these paths with your own.
Relative file paths inside manifests resolve beside the manifest, except for
standalone Python `ldsc` VCF paths, which resolve from the working directory.
Absolute paths work consistently across commands.

## Manifests

### VCF manifest for pipeline

Save as `/data/traits_vcf.csv`. Empty prevalence cells here mean quantitative
traits. The Python pipeline sets `ref=yes` automatically when that column is
absent; if supplied, all values must be `yes`.

```csv
traitname,vcf_files,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,,
Trait_B,/data/Trait_B.vcf.gz,,
Trait_C,/data/Trait_C.vcf.gz,,
Trait_D,/data/Trait_D.vcf.gz,,
```

`prepare` needs only `traitname,vcf_files`; the extra columns above are accepted.
For **standalone Python LDSC**, use `/data/traits_python.csv` with `ref` included:

```csv
traitname,vcf_files,ref,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A.vcf.gz,yes,,
Trait_B,/data/Trait_B.vcf.gz,yes,,
Trait_C,/data/Trait_C.vcf.gz,yes,,
Trait_D,/data/Trait_D.vcf.gz,yes,,
```

Setting every `ref=yes` generates the complete self/pair coverage required by
GPCA. Reference-versus-target screening with some `ref=no` values need not
produce a complete GPCA input.

### Raw-table manifest for GenomicSEM

Save as `/data/traits_native_raw.csv`:

```csv
traitname,sumstats_file,population_prevalence,sample_prevalence
Trait_A,/data/Trait_A_raw.tsv,,
Trait_B,/data/Trait_B_raw.tsv,,
Trait_C,/data/Trait_C_raw.tsv,,
Trait_D,/data/Trait_D_raw.tsv,,
```

An optional positive `N` column sets one sample size per trait during raw
munging, **replacing N in that trait's raw file**. Add `vcf_files` only if you
also want automatic GWAMA preparation; otherwise supply existing GWAMA tables.

### Reusing munged files

For native `--munged_dir`, save `/data/traits_native_munged.csv`:

```csv
traitname,population_prevalence,sample_prevalence
Trait_A,,
Trait_B,,
Trait_C,,
Trait_D,,
```

The directory must contain exactly one `Trait_A.sumstats` or
`Trait_A.sumstats.gz`, and the corresponding file for each other trait. For
native `--munged_input`, add a `munged_file` column with one explicit path per
trait. Choose one mode for the entire manifest.

For standalone Python `--ldsc_only`, use `/data/traits_python_munged.csv`:

```csv
traitname,ref,population_prevalence,sample_prevalence
Trait_A,yes,,
Trait_B,yes,,
Trait_C,yes,,
Trait_D,yes,,
```

Python reuse expects `{traitname}.sumstats.gz`. Preserve associated provenance
files; see [Python reuse rules](PYTHON_LDSC.md#reuse-munged-files-or-restart).
The Python pipeline accepts this manifest too.

### Selecting traits from completed LDSC results

Both GPCA commands need only `/data/selected_traits.csv` unless they will also
prepare VCFs:

```csv
traitname
Trait_A
Trait_B
Trait_C
Trait_D
```

The selected traits must be present in completed results unless you explicitly
allow audited removal. Use a retained-trait manifest after upstream QC drops
traits. Add `vcf_files` for automatic GWAMA preparation.

## VCF fields

These are GWAS summary-statistics VCFs with **one GWAS sample per file**, not
individual-level genotype VCFs. A `.vcf.gz` extension alone does not establish
compatibility.

| Processing step | Required VCF fields | Values used |
| --- | --- | --- |
| GWAMA preparation (`prepare`, or automatic preparation) | Coordinates, REF, ALT; FORMAT `AF,ES,SE,LP,NEF` | EA=ALT; OA=REF; EAF=AF; Z=ES/SE; P from LP; N=NEF. |
| Python LDSC extraction | Coordinates, alleles and IDs; INFO `AF,EUR`; FORMAT `SI,AF,EZ,LP,NEF` | SI for INFO filtering; FORMAT/AF for MAF; INFO/AF versus INFO/EUR for frequency difference; EZ and LP for association statistics. |
| Python binary-trait extraction with population prevalence | The Python fields above, plus FORMAT `NC,NCO` | Uses NC+NCO for N; can infer sample prevalence from case fractions. |

A Python VCF pipeline that prepares GWAMA files needs the union of both field
sets. The preparation path and Python extraction path have separate filtering
rules; passing an LDSC INFO filter does not certify all GWAMA rows.

Preparation handles autosomes 1–22. Its default split layout requires retained
variants on every chromosome for every trait. Use `--splitby_chr nosplit` when
you want one autosomal table per trait, including for data limited to some
autosomes. Neither mode converts genome builds.

## Raw and munged tables

Native GenomicSEM raw input is a whitespace-delimited GWAS table. A clear schema
is `SNP,A1,A2,P,BETA,N`, with a signed `BETA` (or supported signed Z column),
valid P values and positive N. Optional recognized INFO/frequency columns enable
native munging filters. For example, tab-separated:

```text
SNP	A1	A2	P	BETA	N	INFO	MAF
rs101	A	G	0.0455	0.02	50000	0.98	0.20
rs102	C	T	0.0027	-0.03	48000	0.96	0.15
```

These rows illustrate the format, not a sufficient LDSC dataset. Native munging
uses P and effect direction to derive Z. Supply a constant manifest N only when
it is appropriate for that trait.

Compatible munged files contain `SNP,A1,A2,N,Z`. Native reuse requires tab
separation. Prior filtering, allele handling and sample-size conventions remain
your responsibility when reusing files. The initial header check is not a
complete review of how they were created.

## GWAMA tables

Canonical tab-separated columns and a format example:

```text
SNPID	CHR	BP	EA	OA	EAF	N	Z	P
1_100_A_G	1	100	G	A	0.20	50000	2	0.0455003
1_200_C_T	1	200	T	C	0.15	48000	-3	0.0026998
```

Preparation writes exactly these nine columns. The reader also accepts extra
or reordered columns and the aliases `A1`→`EA`, `A2`→`OA`, `p`→`P`.
Use these exact filenames:

- `split`: `{traitname}_chr{CHR}_GenomicPCA_inputs.tsv`, for chromosomes 1–22.
- `nosplit`: `{traitname}_GenomicPCA_inputs.tsv`.

All traits must use compatible SNPIDs, genome builds and allele conventions.
The bundled GWAMA combines the union of SNPIDs and uses available traits at
each SNP. It aligns allele swaps for shared IDs, but does not resolve strand
flips or perform liftover. Matching coordinates do not reconcile different IDs.

## Reference files

Supply an ancestry-appropriate chromosome LD-score directory through `--ld_ref`.
It contains files named `1.l2.ldscore.gz` … `22.l2.ldscore.gz`, with accompanying
`1.l2.M_5_50` … `22.l2.M_5_50` files. `--ld_weights` can point to a separate
regression-weight directory; otherwise the reference directory supplies weights.
Pass the directory, not an individual chromosome file. Native `--chromosomes`
can limit LDSC to files 1 through the requested number; this does not change the
GWAMA split layout.

Python LDSC munging, native raw munging and pipeline use a whitespace-separated
HapMap allele reference with `SNP,A1,A2` headers, supplied through `--hm3`.
Standalone preparation's optional raw-table export accepts a **tab-separated
SNP list** with a `SNP` header; it selects IDs without aligning alleles. A SNP-only
preparation list is not a substitute for the allele reference used by munging.
Reference IDs must match the SNP IDs supplied to munging.

## Sample size and prevalence

There is no universal `--sample_size` flag. Sample sizes enter at different
stages:

| Stage | Source of sample size | How to change it |
| --- | --- | --- |
| Python VCF LDSC, population prevalence absent | FORMAT/NEF; sample prevalence alone is ignored | Supply correctly encoded VCF data. |
| Python VCF LDSC, population prevalence supplied | FORMAT/NC + FORMAT/NCO | Supply valid case/control counts and population prevalence. Missing sample prevalence is inferred from the median per-SNP case fraction. |
| Native raw munging | Raw-table N, unless positive manifest N is supplied | Add per-trait manifest N to replace raw N. |
| Either LDSC backend, munged reuse | N already in each munged file | Prepare appropriate munged inputs; an export override does not change them. Native manifest N is ignored in this mode. |
| Prepared GWAMA inputs | FORMAT/NEF becomes N | Supply correct VCF NEF, or correctly prepared nine-column tables. |
| Optional raw tables from `prepare` | FORMAT/NEF becomes N | These are unmunged tables and contain no INFO. |

For a binary trait, `sample_prevalence` is the fraction of cases in the GWAS
sample; `population_prevalence` is the population prevalence used for liability
scaling. Supplied values must be between 0 and 1. Native LDSC requires **both or
neither**; it does not infer one. Quantitative traits leave both empty/NA.
Choose appropriate N/prevalence conventions before analysis. Python's NEF versus
total-N extraction rule and native N inputs are not interchangeable settings.

## Export metadata

Both GPCA commands automatically export successful GWAMA results. The full
combined result retains the original GWAMA columns. A separate summary contains:

```text
SNPID CHR BP EA OA EAF N_eff BETA SE Z PVAL INFO
```

`--gwama_output_info` supplies a constant INFO between 0 and 1 for **that summary**.
The bundled GWAMA function does not produce INFO, so its full export requires
this explicitly chosen constant. Use your justified downstream policy; this
option does not estimate imputation quality or filter SNPs. Examples use
`${GWAMA_INFO:?...}` so they stop until you supply a value.

`--gwama_output_n_eff` optionally supplies a positive constant N_eff for that
same summary. Omit it to preserve GWAMA's reported values. It does not change
input N, weights, Z, P, BETA or SE. Neither override changes the full combined
result or recalibrates inherited GWAMA formulas. See
[scientific interpretation](REFERENCE.md#gwama-results-and-final-exports).
