# LDSC failure handling and permissive input columns

Validated on 2026-10-01. Statistical regression code and the supplied GWAMA
implementation were not changed. The change separates numerical-result
collection from eligibility for genomicPCA, preserving strict defaults.

## Reproduction

Use Python 3.11 with the declared dependency ranges and R with argparse,
data.table and glue:

```bash
PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 \
  LDSC_GPCA_TEST_RSCRIPT="$(command -v Rscript)" \
  python -m unittest discover -s tests -v
```

Test environment: Python 3.11.15, NumPy 1.26.4, pandas 2.3.3, Polars 0.20.31;
R 4.4.2, data.table 1.17.8, argparse 2.3.1, glue 1.8.0. The default host
Python/dependency versions were newer than the package supports, so final
validation used a separate temporary environment with supported versions.

Observed results: the complete suite ran 100 tests in 179.964 seconds, with
97 passing and three optional native-runtime tests skipped. The additional
independent-batch-retention test, added after that suite started, also passed
in the same environment (98 passing tests in total). No production code
changed after the full suite started.

## Expected and observed behavior

- Three simulated traits, complete pair coverage, one negative self h2 and five
  directed missing-rg rows: strict mode saves diagnostics then fails; report
  mode preserves all nine rows, the five missing correlations and every valid
  numerical estimate at binary64 precision. It does not publish a success CSV.
- Real Python-to-R handoff: strict R mode records failure before GWAMA; opt-in
  removal retains the two valid traits in original order and passes both
  correlation- and covariance-PCA validation. No correlation is imputed.
- Matrix checks: complete upper triangle, reversed pairs, consistent/conflicting
  duplicates, missing pairs, unknown traits, shuffled input, intercept mapping,
  unchanged PC1 when only diagnostic columns change, and both GWAMA naming modes.
- Boundary/error checks: p=0/1, tiny positive SE, non-positive/non-finite SE,
  non-finite rg/h2, non-numeric text, all-missing estimates, missing fields/files,
  duplicate headers, fewer than two retained traits and non-positive-definite
  CTI. Structural failures remain errors in report mode.
- Failed process and missing/stale export fixtures: bounded retries, per-batch
  status, retained successful files. Failed recompilation archives previous
  successful outputs rather than presenting them as current.
- Additional and reordered columns in manifests, numerical LDSC tables and
  split/whole-genome GWAMA inputs: accepted; canonical fields remain authoritative.
  The GWAMA reader supplies only its nine required columns in canonical order.
  Missing/ambiguous required fields still fail. Former manifest aliases are
  annotations, not substitutes for required names.
- Observed/liability alternatives and mixed-scale rows: unused empty scale fields
  do not become false failures. Existing scale mapping remains unchanged.
- One existing installer test fixture was updated to answer the installer's
  `env create --help` capability probe. The installer itself was unchanged.

## Saved failed-run replay

The locally supplied 19-trait failed run was processed without rerunning LDSC:
all 76 current numerical exports yielded 361 directed rows; all 37 missing-rg
rows were retained and audited. Explicit R trait removal excluded the one
unestimable trait and validated the complete 18-trait subset. The original
inputs were not edited. Existing diagnostics remained visible: three traits
with self h2/SE below 2, two unique off-diagonal rg estimates above 1, and six
negative genetic-correlation eigenvalues. CTI/PC1 validation passed under the
existing warning policy; the input matrix was not replaced or clipped.

Real input files and derived results are deliberately not committed. GWAMA was
not run on this biological dataset; the replay used `--validate_only`.

## Limits and sources

The three optional fresh-native-LDSC integration tests require an explicitly
configured compatible LDSC runtime and were skipped. Saved native numerical
results and real R validation were exercised; no claim is made that a fresh
native regression or a full biological GWAMA analysis was rerun.

- [Pinned CBIIT LDSC negative-heritability handling](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py#L697)
- [LDSC heritability and correlation FAQ](https://github.com/bulik/ldsc/wiki/FAQ#genetic-correlation-and-heritability)

The existing deterministic trait-removal heuristic is reused. It does not
guarantee the largest complete subset and changes the analyzed trait set;
exclusion audits and the retained manifest must accompany interpretation.

## Restart validation (2026-10-01)

`tests/test_restart.py` adds 16 tests in the same supported environment. All 16
passed, including a real CLI/filter/compilation path with a simulated native
regression command and runtime identity. Native numerical estimates in these
tests are explicit fixtures, not newly estimated biological results.

The full regression suite after this change ran 117 tests in 203.299 seconds:
114 passed and the same three optional native-LDSC tests were skipped. Existing
R matrix, interface and worker validation tests passed. `git diff --check` passed.

Expected behavior matched observed results:

- Four completed fixture batches: restart launched zero regression commands,
  preserved result ordering and returned four `reused` statuses. Without
  `--restart`, all four ran again.
- One process interrupted after writing its result: only the three independently
  completed batches had checkpoints. Restart launched exactly one regression.
- Changed sumstats contents with the original modification time restored:
  three affected batches reran and the independent B-B result was reused.
- Changed LD scores, M files, regression weights, chi-square threshold, extraction
  filters, either prevalence, or runtime source identity: affected results reran.
- Identical decompressed sumstats with different gzip timestamps and extra
  manifest annotations: all four results were reused. The real CLI re-filtered
  ten-SNP inputs twice and produced byte-identical compiled results on restart.
- Missing, empty or modified results; absent/corrupt checkpoints; changed batch
  layout; malformed numerical text; missing required columns; incomplete or wrong
  pair coverage: no incorrect completed result was reused. A truncated CSV with
  a matching recorded checksum still failed pair-coverage validation.
- Numerical missing-rg estimates remained unchanged and reportable during reuse.
- Inputs changed during execution: no completion checkpoint was saved. Missing
  references and invalid gzip inputs stopped before regression commands launched.
- Runtime probe parsing, failure handling, native-source changes and CLI argument
  propagation were tested. The Python LDSC runtime itself remains unavailable
  for fresh integration testing on this host.

Hashing and reference-file selection were checked against the
[Python hashlib documentation](https://docs.python.org/3/library/hashlib.html)
and the pinned
[CBIIT LDSC reference parser](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/parse.py).
Restart relies on the existing numerical reader and exact required batch
coverage, and does not change regression, PCA or GWAMA calculations.

## LDSC-stage trait exclusion (2026-10-01)

The LDSC command now accepts `--result_failure_action drop_traits`. Its Python
selection helper matches the existing R heuristic without requiring R to run
LDSC. Every manifest trait must have `ref=yes`; at least two traits with usable
self estimates and complete unordered pair coverage must remain. Outputs include
the retained manifest and an ordered exclusion report. Extra manifest columns
and the original trait order are preserved.

All 14 tests in `test_ldsc_drop_traits.py` passed. The full suite then ran
131 tests in 194.502 seconds: 128 passed and the same three optional fresh-Python
LDSC tests were skipped. The R parity test exercised seven selection scenarios
and strict matrix validation. Expected and observed behavior agreed:

- A failed self estimate excludes only that trait; every retained numerical
  estimate is unchanged. Healthy traits, low positive h2/SE and finite rg outside
  [-1,1] remain eligible, with the existing diagnostic warnings.
- Missing self-pairs, missing unordered pairs, failed-pair degree and both
  tie-breaks produce the same retained order as the existing R implementation.
  One valid observed orientation covers a missing opposite orientation.
- Consistent duplicate estimates are accepted. Conflicting estimates, including
  finite-versus-missing conflicts, remain fatal even if exclusion could hide
  them. Absolute duplicate tolerances match R, including for large z values.
- Zero successful batches, all traits failing or fewer than two retained traits
  save diagnostic/exclusion reports and fail without publishing a success CSV.
  Malformed numeric fields, unknown traits and unsupported reference manifests
  remain errors. Mixed observed/liability scales select the appropriate h2 Z.
- CLI tests use real filtering, scheduling and compilation with simulated native
  commands. Recoverable execution failures are retried, then missing comparisons
  are audited before selection. Malformed numerical exports remain fatal.
- Restart reuses all completed numerical-failure fixtures without launching a
  regression. In an execution-failure fixture it retries the five failed batches
  and reuses the four completed batches. Previous retained manifests are archived
  before recompilation, so stale selections cannot masquerade as current output.

The saved 19-trait biological run was replayed through the new LDSC policy using
its 76 existing numerical exports. It excluded only
`AMDHD2_Q9Y303_OID30366_v1_Cardiometabolic_II`, retained 18 traits and wrote all
324 directed retained rows. The diagnostic CSV preserved all 361 original rows,
including the 37 missing correlations. An exact DataFrame comparison confirmed
that the selected estimates equal the corresponding diagnostic rows.

Strict R `--validate_only`, using the new retained manifest and result file with
the default error policy, passed. The correlation matrix, CTI matrix and PC1
weights exactly matched the earlier R-stage trait-removal outputs. Existing
low-h2-Z, out-of-range-rg and negative-eigenvalue warnings remained visible.
Original biological inputs and R/GWAMA implementation files were unchanged.

This validates saved native outputs, the selection policy, restart behavior and
the real R handoff. A fresh native regression and full biological GWAMA run were
not executed. The heuristic does not guarantee the largest possible subset, and
a failed batch is recorded as an unavailable comparison rather than evidence of
a biological defect in a trait.
