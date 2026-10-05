# Package cleanup validation

Baseline: `e4d2adedac47671f291e732612069365bd345df0` (2026-10-03).

This change optimizes package orchestration, parsing, QC and PCA reporting. It
does not edit the original LDSC, GenomicSEM or vendored GWAMA implementations.
Thirteen protected vendor/backend/adapter files matched the baseline SHA-256
hashes. Regression settings, normalization formulas, CTI, filtering rules,
trait ordering, retry policy and output schemas remain unchanged.

## Changes

- Use literal AWK patterns for the two INFO tags in the existing AF audit.
- Read LDSC trait identifiers as character, preserving names such as `001`.
- Share R row preparation between public QC functions; internal self-pair checks
  reuse prepared rows. Public entrypoints retain validation. Use reference-safe
  metadata setters and avoid reconverting plain double columns.
- Index Python self-pairs and aggregate failed-pair counts once; retain the
  original first-row/mean distinction, missing values, reason order and manifest
  order. Reuse the trait summary when publishing a retained subset.
- Share PC1 loading, sign and numerical-check results between computation and
  reporting. Compute the loading vector directly from the first eigenvector
  and eigenvalue; retain the full eigendecomposition and failure reports.
- Publish existing R audit tables from named lists and reuse the existing
  content-hashing helper for munging provenance.

The implementation has 117 fewer physical lines under `src/`; tests and this
documentation are additional. Binary TSV serialization, file-publication
semantics, status-update frequency, CLI defaults and upstream algorithms were
not rewritten. Those broader I/O changes require separate profiling and tests.

## Real-data equivalence

The original saved LDSC tables were read from
`/Users/JJOHN41/Downloads/python_gmenomicsem_comparison/python/`:

| Dataset | Traits | Directed LDSC rows | GWAMA subset SNPs |
|---|---:|---:|---:|
| T3_C0_C1_C2 | 43 | 1,849 | 521 |
| T3_C0_C1_NOISE | 49 | 2,401 | 259 |

For each dataset, run the baseline and changed package with all four combinations
of correlation/covariance PCA and pair/trait-wide normalization. Result
compilation, matrices, loadings, eigenvalues and 20–22 CSV audit tables per run
matched exactly using round-trip float parsing and exact frame comparison.
Every GWAMA output column also matched exactly, including Z, P, BETA, SE and N.

Saved GenomicSEM LDSC RData files were tested with both
correlation and covariance PCA on each dataset. All 19 compared CSV reports per
run and every GWAMA output column matched exactly. Thus the comparisons cover
12 analysis configurations and 244 diagnostic CSV tables, plus compilation and
GWAMA results. Run-status fields containing different output paths and binary
R session/settings artifacts were not used for numerical equality checks.

The complete original LDSC tables were used. The GWAMA inputs are existing
reconstructed real Z/N subsets with shared reference EAF, not the unavailable
original full GWAS input files. This validates computational equivalence on
those subsets; it does not independently establish BETA/SE calibration or
replace a new whole-genome regression run. Input hashes were unchanged.

## Performance and boundary checks

On a 300,000-record synthetic compressed VCF, one warm-up followed by three
measured runs per version gave median extraction times of 13.001 seconds for
the baseline and 0.528 seconds for the cleanup (about 25 times faster).
Execution order alternated. Output TSV hashes and audit text were identical.
Other local validation jobs were running, so these are illustrative local
stage timings, not an isolated system benchmark or whole-pipeline speed claim.

The audit fixture includes missing INFO/AF, AF disagreement, low SI and values
near filtering boundaries. The extraction tests additionally cover binary
traits, multiallelic missing values, missing headers, stale audits, empty input,
MHC exclusion and two-attempt worker failures.

Instrumentation on both real LDSC tables confirmed:

| Operation per validation workflow | Before | After |
|---|---:|---:|
| Heritability normalization | 5 | 3 |
| Numeric conversion entrypoints | 4 | 2 |
| Self-pair QC | 2 | 2 |
| Duplicate validation | 2 | 2 |

Selection and final validation keep their independent checks. The cleanup avoids
repreparing the whole table inside each self-pair check.

New regression cases cover numeric-looking identifiers across CSV/gzip chunk
boundaries, input immutability, invalid numeric values and missing columns,
duplicate self rows with missing values, empty summaries, trait order, PC1 sign
orientation, negative eigenvalues and failure-report creation for invalid PC1.

## Reproduction

Local evidence and scripts:
`../package_cleanup_validation_20261003/` beside the checkout. This contains
the `e4d2ade` baseline snapshot, original-input hashes, before/after results,
`summary.json`, `native_gpca_summary.json`, `audit_regex_benchmark.json`,
`profile_before/`, `profile_after/`, and protected/validated source hashes.
The evidence directory is outside the repository and is not committed.

Focused portable regression command from the repository root:

```sh
PYTHONPATH=src:tests LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  python -m unittest test_cleanup test_extraction_af test_failure_handling \
  test_ldsc_drop_traits test_self_pair_zero_se test_normalization
```

The local `run_tests.py` selects tracked Python test modules plus the new cleanup
module, avoiding unrelated ignored tests in the workspace. The final suite log
and machine-readable result are `tests_final.log` and `tests_summary.json`.
The final run executed 187 tests: 184 passed, 3 skipped, no failures or errors.
The skipped checks require an explicitly configured isolated Python LDSC runtime
(native table export, three-trait regression and separate-weight regression).
R tests and native munging integration were enabled and passed.

Environment: macOS ARM64, Python 3.12.14, pandas 2.3.3, Polars 0.20.31,
R 4.4.2, data.table 1.17.8, bcftools/htslib 1.23 and GNU Awk 5.3.1.
Performance on the pinned deployment environment has not been measured.
