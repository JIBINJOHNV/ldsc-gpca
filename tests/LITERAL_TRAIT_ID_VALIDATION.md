# Literal trait identifiers at R boundaries

The shared GPCA manifest reader, Python LDSC table reader and GenomicSEM LDSC
manifest reader preserve literal `NA`, `NaN`, `nan` and `001` identifiers.
Missing tokens are interpreted in numeric fields. Empty/whitespace identifiers
remain invalid; malformed numeric text remains a structural error.

## Reproduction and regression coverage

`test_literal_trait_ids.py` runs the real R readers and validation functions.
All three test groups failed on the preceding code and pass with this change.
The fixtures contain five traits in a deliberately nonalphabetical order:
`NaN`, `001`, `Trait_A`, `NA`, `nan`.

- Quoted and unquoted manifests preserve identifiers and text annotations.
  Empty/whitespace names and duplicate names fail.
- A complete 15-row LDSC triangle is shuffled and read in two-row chunks.
  All 24 combinations of quoting, CSV/TSV, plain/gzip and observed/liability/mixed
  heritability retain every row and the exact manifest order. Correlation, CTI,
  SE matrices, self heritabilities and PC1 loadings match the in-memory baseline.
  Enabling the missing-trait policy does not remove these present traits.
- Numeric `NA`, `NaN`, `nan` and empty fields in rg, SE and h2 remain reportable
  estimation failures: strict mode rejects them; explicit exclusion removes only
  the affected trait. Nonnumeric text fails under both policies. These checks run
  with both quoting styles, for 30 input cases. Empty/whitespace p1/p2 values fail
  before row selection, with the field and file row identified.
- The real GenomicSEM wrapper runs both LDSC passes with an injected estimator.
  Both quoting styles preserve exact names/order at the estimator handoff, in
  the selected manifest and in the saved RData. Missing prevalence/N fields and
  valid numeric metadata retain their meaning without conversion warnings.
  Invalid empty names fail before the estimator is called.

The GenomicSEM test validates ingestion and orchestration, not an installed
GenomicSEM regression. No production cohort analysis was rerun, and no vendor
LDSC/PCA/GWAMA calculation was changed.

## Running the checks

From the repository root, set `LDSC_GPCA_TEST_RSCRIPT` to an Rscript with the
project's R dependencies, and use a Python environment with its dependencies:

```bash
PYTHONPATH=src:tests python -m unittest test_literal_trait_ids -v
```

Related regression suites are `test_failure_handling`, `test_self_duplicates`,
`test_self_pair_zero_se`, `test_interface_aliases`, `test_normalization`,
`test_pipeline` and `test_ldsc_drop_traits`.

On 2026-10-04, the combined eight-module run passed all 97 tests with no skips,
using Python 3.12 in `.venv-ldsc-normalization` and R 4.4.2. The three new test
groups also passed after adding final metadata and self-heritability assertions.
The installed data.table emitted its existing R-build-version warning; this did
not fail the checks. Logs for this checkout are `/tmp/literal_trait_ids_before.log`,
`/tmp/literal_trait_ids_after.log` and `/tmp/literal_trait_ids_regression.log`.

Reader semantics: [data.table fread](https://rdatatable.gitlab.io/data.table/reference/fread.html)
and [R read.table/read.csv](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html).
