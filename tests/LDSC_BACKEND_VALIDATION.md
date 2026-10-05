# LDSC backend selector validation

Validated on 2026-10-04. `ldsc-gpca ldsc --ldsc_backend` selects `python`
(default) or `genomicsem`, reusing the existing parsers and analysis workflows.
Standalone `prepare` retains `--munge_backend`; `ldsc` and `pipeline` reject it.

## Checks and observed results

| Check | Expected and observed behavior |
| --- | --- |
| Default and explicit Python selection | Forward the same arguments to Python LDSC, including `--chisq_max auto`, and preserve its exit status. |
| GenomicSEM raw, VCF and both reuse modes | Dispatch to the existing GenomicSEM workflow with unchanged input-mode arguments and numerical settings. |
| Existing `genomicsem ldsc` entry point | Continue to dispatch to the same GenomicSEM workflow. |
| Empty invocation, selected backend alone, `-h` and `--help` | Show the selected backend's options and outputs without starting analysis. |
| Invalid, missing, empty, abbreviated or conflicting backend choices | Exit with an argument error before either analysis starts. Consistent repeated choices are accepted. |
| Options from another backend or a separate munging override | Reject before either analysis starts. GenomicSEM also rejects incompatible input-mode combinations. |
| Parser reuse | Backend defaults are unchanged after repeated selection; the standalone Python parser and pipeline parsers are not mutated. |
| Preparation | Omission still selects Python munging; explicit GenomicSEM selection remains available; preparation rejects `--ldsc_backend`. |
| Python reuse handoff | Two-trait fixture preserves trait names/order, sample-size conventions and separate LD weights through the unified CLI. Regression/runtime and collection are stubbed. |
| GenomicSEM reuse handoff | Two-trait fixture preserves paths, trait order, prevalence, worker count and separate weights in the R invocation. R execution is stubbed. |
| GenomicSEM VCF handoff | Real bcftools conversion reaches the existing R entry point in `munge` mode with the selected HapMap reference. Only the R execution is stubbed. |

The following suite ran **74 tests, all passed**, with the R parser and
orchestration checks enabled:

```bash
PYTHONPATH=src:tests LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  python -m unittest test_ldsc_dispatch test_pipeline test_interface_aliases \
  test_genomicsem_vcf test_threads_and_workers
```

Standalone preparation ran **13 tests: 12 passed, 1 skipped**. The skipped test
requires an installed GenomicSEM package. It was run without
`LDSC_GPCA_TEST_RSCRIPT` because that package is unavailable in this environment:

```bash
PYTHONPATH=src:tests python -m unittest test_prepare_modes
```

Additional checks passed for 11 actual CLI help invocations, all 18 documented
LDSC analysis commands (9 Python and 9 GenomicSEM, parser checks only), and 340
local documentation links/anchors. `git diff --check` passed.

## Evidence boundary

These checks validate selection, argument handling, existing preparation,
orchestration and file handoff. They do not establish full installed GenomicSEM
regression correctness or equivalence between the estimators. No statistical
implementation was changed for this selector. Reused files still skip munging;
backend selection does not harmonize alleles or add enforcement of GenomicSEM
reuse provenance. Those input requirements remain documented separately.

The dispatcher uses Python's documented
[`parse_known_args`](https://docs.python.org/3.12/library/argparse.html#partial-parsing)
to select a backend, then validates the full argument list with a copy of its
existing parser before forwarding the remaining arguments.
