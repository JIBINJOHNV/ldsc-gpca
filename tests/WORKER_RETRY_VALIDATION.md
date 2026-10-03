# Worker retry validation — 2026-10-02

Package-managed parallel jobs default to two total attempts: the initial attempt
and one retry. Successful jobs are retained; only failed/unconfirmed jobs retry.
Exhaustion stops downstream analysis, including when Python LDSC collection uses
`--result_failure_action drop_traits`. That policy still applies to numerical
estimation failures from completed commands. The existing `--ldsc_retries`
override is preserved; other managed stages have one retry.

## Automated checks

Run from the repository root with its Python dependencies and an R installation
containing the package's R dependencies:

```bash
PYTHONPATH=src LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  python -m unittest discover -s tests -v
```

Observed: **167 tests ran; 164 passed, 3 skipped** (166.660 seconds). The skipped
tests require an explicitly configured native Python LDSC child runtime, which
was unavailable in this test environment. Python 3.12.14, pandas 2.3.3, Polars
0.20.31 and R 4.4.2 were used. An additional check of recovered GenomicSEM munging
was then added inside the existing R test; the seven retry tests passed again
(10.581 seconds):

```bash
PYTHONPATH=src LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  python -m unittest discover -s tests -p 'test_worker_retries.py' -v
```

After updating CLI help, all ten tests in `test_threads_and_workers.py` also
passed (93.701 seconds), including the real R parsers and subprocess settings.

| Case | Expected and observed |
| --- | --- |
| Successful job alongside a transient failure | Successful job runs once; failed job runs twice; manifest order retained. |
| Persistent Python worker error | Two attempts, named job/error in audit, fatal exception. |
| Python child actually exits via `os._exit` | Broken pool replaced; job succeeds in fresh pool on attempt 2. |
| R child returns `NULL`, `try-error`, malformed data, or an ordinary error | Explicit failed status followed by successful retry; status reporting does not crash. |
| R child actually killed with SIGKILL | Controller survives; only failed job retries successfully. |
| Whole-genome R worker repeatedly fails | `whole_genome` status retained with two attempts and failure details. |
| All 22 chromosome workers fail in either GPCA backend | Both attempts recorded; both workflows stop with chromosome/error details. |
| Python munging interruption | Only failed trait retries; persistent failure stops instead of silently dropping the trait. |
| Unchanged old munging output | Rejected after both attempts; no stale success accepted. |
| Native GenomicSEM munging orchestration | Failed traits retry; permanent failures stop before LDSC; recovered traits reach LDSC in manifest order. |
| Audit cannot be written | Audit error propagates; analysis cannot silently report success. |
| Empty job list / invalid Python worker count | Empty list returns no jobs; invalid counts fail clearly. |
| LDSC execution failure with `drop_traits` | Five failed batches each attempted twice; four successful batches run once; compilation and selected-manifest publication do not run. |
| Existing matrix, ordering, preparation and export checks | Passed in the full suite. |

GenomicSEM is not installed in this environment. Its munging tests inject a
function with the upstream API, including a successful gzip output, and check
scheduling, arguments, freshness, order and the LDSC handoff. They do not test
the native GenomicSEM statistical implementation. The GenomicSEM GPCA failure
test also replaces the read/QC boundary with a small validated matrix fixture;
it exercises the actual GWAMA orchestration and fatal error handling.

## Real-data GWAMA checks

Reused the existing native LDSC tables and reconstructed real Z/N GWAMA subsets
from the earlier normalization/covariance validation. No new LDSC regressions
were fitted for this check.

| Dataset | Traits | SNPs | Normal run | First child killed | Persistent worker failure |
| --- | ---: | ---: | --- | --- | --- |
| T3_C0_C1_C2 | 43 | 521 | Success, 1 attempt | Success, 2 attempts | Stopped after 2 attempts; export skipped |
| T3_C0_C1_NOISE | 49 | 259 | Success, 1 attempt | Success, 2 attempts | Stopped after 2 attempts; export skipped |

For both normal and recovered runs, SNP IDs, Z scores and P values matched the
previous results exactly. The correlation matrix, CTI and PC1 weights also
matched exactly. SHA-256 checks confirmed all source GWAMA input files were
unchanged. The retry fixture kills the first R child before calling the bundled
estimator, then calls the unchanged estimator on attempt 2. The permanent-failure
fixture raises an error on both attempts. These checks ran the complete Python
`gpca` launcher, including automatic export for successful runs.

The reconstructed inputs lack INFO. Successful export checks therefore explicitly
used `--gwama_output_info 0.9` as **test-only metadata**, not as an inferred INFO
score or a recommendation for downstream scientific use. The original missing-INFO
run correctly failed both export read attempts. This override changes neither
the raw GWAMA result nor LDSC, PCA or GWAMA calculations. These are small real-data
subsets; whole-genome performance and reconstruction-dependent BETA/SE/N_eff were
not validated by this comparison.

Local evidence is retained beside the checkout in
`worker_retry_validation_20261002/`: `validate.py`, `summary.json`, per-dataset
normal/retry/failure logs, result files and worker audits. Full and targeted test
logs are `worker-retries-suite.log`, `worker-retries-new-tests.log` and
`worker-retries-cli-tests.log` in the same parent directory. Private datasets and
generated outputs are not bundled in Git.

## Behavior and limits

The stage collects submitted worker outcomes before reporting final failure;
it does not immediately kill other running jobs. Retried jobs outside pairwise
Python LDSC run one at a time to reduce concurrent memory demand. Unix/macOS R
jobs run in child processes even for a single job; Windows uses sequential R
execution. Terminating the controller itself cannot be recovered within that
process. Preflight errors, sequential publication and archival are not worker
retries. No upstream LDSC, GenomicSEM or GWAMA estimator was modified.

Native GenomicSEM munging is now called per trait with upstream `parallel=FALSE`;
the package owns scheduling. This can repeat the upstream HapMap reference read
per trait. Its performance has not been benchmarked with the native package.

Verified API and failure semantics against the official sources:

- [Python concurrent.futures](https://docs.python.org/3/library/concurrent.futures.html)
  describes broken process pools and future exceptions.
- [R mclapply](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)
  describes `try-error`/missing child results and prescheduling.
- [R mcparallel/mccollect](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mcparallel.html)
  describes isolated child execution and result collection.
- [Pinned GenomicSEM munge source](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge.R)
  supplies the single-trait API used by the retry wrapper.
