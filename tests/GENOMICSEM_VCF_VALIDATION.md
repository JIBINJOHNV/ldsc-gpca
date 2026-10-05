# GenomicSEM VCF input validation

Validated on 2026-10-04 using Python 3.12.14, Polars 0.20.31, pandas 2.3.3,
NumPy 1.26.4, real bcftools and R 4.4.2 on macOS ARM64.

## Scope and results

The new `genomicsem ldsc --vcf_input` mode converts GWAS VCFs to raw tables,
then uses the existing native munging/LDSC runner. Pipeline can select this
mode explicitly. Existing raw-table, munged-file and legacy pipeline modes
remain available. Original Python LDSC, GenomicSEM R and GWAMA sources were
unchanged; hashes of all 53 other existing package source files matched the
pre-change baseline.

**75 Python tests passed, with no skips.** New conversion tests use synthetic
VCFs with known values and real bcftools, rather than mocked extraction. They
check quantitative NEF, explicit binary N without a NEF header, signed effects,
alleles, MAF, INFO, P flooring, missing/invalid fields, duplicate variants,
compressed input, manifest order, untouched input files and two-attempt
failures. CLI tests check the converted manifest passed to R; R execution is
stubbed only in those CLI tests. Pipeline tests check explicit VCF routing,
preflight rejection, mutually exclusive modes and validation-only behavior.
The selected suite also checks existing preparation, raw/munged interfaces,
I/O, retries and thread handling.

The existing R tests passed all ten matrix acceptance cases and the native
two-pass handoff fixture. The latter injects estimates to check orchestration;
it is not a fresh LDSC regression.

## Checks against pinned native munging

The unmodified [munge_main.R](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge_main.R)
and [utils.R](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/utils.R)
were sourced directly in R. Converted fixture tables were compared with
independently specified raw tables, using SNPs rs1/rs2, P values 1e-4/1e-300,
negative/positive effects, MAF 0.2 and INFO 0.95/0.8.

| Input | INFO threshold | Expected retained SNPs and N | Observed |
| --- | --- | --- | --- |
| Quantitative, N from NEF | 0.9 | rs1; N=20000 | Exact munged-table match |
| Quantitative, N from NEF | 0.8 | rs1/rs2; N=20000/22000 | Exact munged-table match |
| Binary, explicit N | 0.9 | rs1; N=35000 | Exact munged-table match |
| Binary, explicit N | 0.8 | rs1/rs2; N=35000/35000 | Exact munged-table match |

This checks compatibility with native munging, including effect direction and
INFO selection. It does not establish an appropriate N/prevalence convention
for an arbitrary binary study. Native Z is reconstructed from P and direction;
the adapter's P floor can therefore affect extreme Z. bcftools decoding can
also round VCF numeric fields; the adapter does not restore textual precision.

## Reproduction and limits

Run from the repository root with bcftools and Rscript on PATH:

```bash
PYTHONDONTWRITEBYTECODE=1 POLARS_MAX_THREADS=4 PYTHONPATH=src:tests \
  python -m unittest test_genomicsem_vcf test_pipeline \
  test_interface_aliases.CliNamesTests test_interface_aliases.ManifestNamesTests \
  test_io_acceleration test_worker_retries test_threads_and_workers
Rscript tests/test_failure_handling.R "$PWD"
Rscript tests/test_pipeline_native.R "$PWD"
```

Documentation checks passed for 319 local links/anchors, five GenomicSEM LDSC
examples, nine pipeline examples and all 19 GenomicSEM LDSC long options. Both
backends' pipeline options are covered in the pipeline guide.

The full GenomicSEM package was unavailable in this environment. No fresh
VCF-to-LDSC regression or original-user-VCF run was performed. These results
validate conversion, native munging compatibility and orchestration, not a
complete real-data analysis, whole-genome performance or cross-platform use.
