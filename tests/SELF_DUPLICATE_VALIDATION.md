# Repeated self-estimate validation

Validated on 2026-10-04 with Python 3.12 and R 4.4.2.

## Change

Python result collection checks duplicates before scale mapping, derived
annotations, trait removal and publication, under `error`, `report` and
`drop_traits`. Conflicts identify the pair, fields, source files and CSV rows.
R GPCA input QC applies the same self-only field checks before aggregation or
trait removal. Original observed/liability columns remain available through
normalization so an unselected scale cannot hide a conflict.

Repeated self-pairs must agree in heritability, heritability SE, heritability
intercept and intercept SE. Supplied observed and liability scales are compared
separately. Existing common pair-field checks still apply; off-diagonal reverse
rows may report different target-trait heritabilities/intercepts.

Each original self estimate is checked for validity. Consistent but invalid
duplicates follow the existing estimate-failure policy; conflicting duplicates
are always structural errors. Positive means cannot rescue negative source
heritabilities, including mixtures whose differences fall within tolerance.

## Executed checks

From the repository root:

```bash
PYTHONPATH=src:tests \
LDSC_GPCA_TEST_RSCRIPT=/Users/JJOHN41/miniconda3/bin/Rscript \
../.venv-ldsc-normalization/bin/python -m unittest \
  test_self_duplicates test_self_pair_zero_se test_ldsc_drop_traits \
  test_failure_handling test_results_csv test_normalization -v

PYTHONPATH=src:tests \
LDSC_GPCA_TEST_RSCRIPT=/Users/JJOHN41/miniconda3/bin/Rscript \
../.venv-ldsc-normalization/bin/python -m unittest \
  test_interface_aliases test_io_acceleration -v
```

Observed: **62 core tests and 23 interface/I/O tests passed; no skips**.
The Python suite executes the R regression scripts as subprocesses. Local
interpreter paths above should be replaced with equivalent environments elsewhere.

| Representative input | Expected and observed result |
| --- | --- |
| Each of six self fields changed separately, including non-selected liability fields | Conflict under all Python policies and both R selection policies; no final Python result published |
| Conflicting copies in two CSV files | Error identifies both source files and row numbers |
| Finite versus missing, positive infinity or negative infinity | Duplicate conflict |
| Difference exactly 0.001 / just above 0.001 | Accepted / rejected; R custom tolerance respected |
| Identical self duplicates; entirely empty alternative scale | Accepted; compiled numerical output unchanged |
| Self h2 values -0.2 and +0.8 | Structural conflict, never averaged into a usable result |
| Self h2 values -0.0001 and +0.0003 | Invalid source detected despite positive mean and within-tolerance difference; error/report/drop policies behave as specified |
| Identical invalid self estimates or SEs | Invalid estimates remain invalid; R cannot rescue them by averaging |
| Reverse off-diagonal rows with different h2/intercept fields | Accepted |
| Repeated self-rg SE=0, Z=+Inf, P=0 | Existing exception preserved and audited |
| Shuffled complete three-trait data in manifest order C, A, B | Identical correlation, CTI, SE matrices, covariance matrices and PC1 loadings with consistent duplicates |
| Three aligned SNPs across those traits | Bundled modified GWAMA executed for correlation and covariance PCA; all result columns identical with and without consistent duplicates |
| Conflicting rows read in separate R input chunks | Conflict survives normalization and is rejected |

The collection test previously expecting conflicting rows to reach R now expects
an immediate error. An unrelated stale test regex was aligned with the existing
`Incomplete Python LDSC` message. Documentation and duplicate-tolerance help were
updated. `git diff --check` passed.

## Evidence boundary

These are representative synthetic regression fixtures, including execution of
the bundled GWAMA function. No production cohort was rerun, and this validation
does not establish whether earlier production results contained conflicting
duplicates. LDSC estimators and GWAMA weighting code were not modified.
