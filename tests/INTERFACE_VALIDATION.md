# Single-interface validation — v0.6.0

Validated on 2026-09-29. This release removes alternative managed option and input
manifest spellings, uses underscores for multiword managed options, and disables
R option abbreviation. Variant headers, upstream raw-script options, GWAMA and
the statistical calculations are unchanged.

## Tests and observed results

| Check | Data / expected behavior | Observed |
| --- | --- | --- |
| 16 interface tests (15 Python tests and one real R test harness) | Two-trait manifests with names `002`,`001`; single names accepted; previous spellings rejected even alongside valid names; malformed/duplicate/missing input rejected; row order/leading zeros retained | Passed |
| Native GenomicSEM modes | Unmunged paths, munged paths and munged directory; quantitative and binary prevalence rows; partial prevalence rejected; canonical reference/weight/worker options forwarded to a captured R invocation | Passed; R estimation itself not run in this invocation test |
| Both real R parsers | Current GPCA flags parse; removed flags and abbreviations fail; manifest names and unsupported headers checked | Passed |
| 10 worker/thread tests | Defaults, overrides, early initialization, child environment inheritance, raw passthrough, real R worker values and invalid inputs | Passed |
| 39 existing precision/export/runtime tests | Existing result reading, precise numeric export, compilation, invocation and failure cases; real bcftools extraction fixture | Passed |
| 17 scientific regression checks | Complete/reversed/duplicate pairs, missing/conflicting input failures, CTI mapping, ordering, diagnostics-only fields, PCA and analytic GWAMA behavior | Passed |
| Isolated built-wheel installation | Version/help for all managed modes; additional help sections; explicit rejection of removed flags; packaged source/resources match tested source | Passed |

The first Python interface test run exposed an incorrect test expectation for
macOS's `/var` versus `/private/var` temporary-path symlink. The test now compares
resolved paths, consistent with the existing implementation; its rerun passed.

The 17 scientific checks include two reproductions of known inherited edge cases
(invalid EAF and zero-loading-only SNPs). Passing these checks establishes unchanged
behavior; those pre-existing GWAMA issues were not repaired in this naming task.

One unrelated installer mock test (`test_setup_routes_and_stops_on_failure`) was
excluded. Its failure was previously reproduced on unmodified v0.5.3; installer
behavior was not changed by this release.

## Scope and reproducibility

The following modules match the v0.5.5 baseline byte-for-byte: pairwise LDSC,
munging, extraction, native GenomicSEM LDSC pipeline, shared matrix validation,
shared PCA, shared GWAMA and the original vendor GWAMA function. The vendor SHA256
remains `99c041b1e31366e62f1650088e750615d2284047163a55edceeafdf00fcb34d2`.
No variant table was renamed or rewritten. The existing 15 templates retain their
schemas; input guides and command examples use only current names.

Real native CBIIT regressions were already tested for v0.5.5, including separate
weights. They were not rerun for this name-removal update, which does not change
those modules. Full native GenomicSEM estimation, Nextflow execution and the
180-trait production analysis were not rerun. No speedup claim is made.

Run with declared Python dependencies:

```bash
python -m unittest discover -s tests -p test_interface_aliases.py -v
python -m unittest discover -s tests -p test_threads_and_workers.py -v
```

The interface test filename is retained; its tests now assert rejection of
alternative spellings. Set `LDSC_GPCA_TEST_RSCRIPT` to R with `argparse` and
`data.table` to run its R checks. Local tests used Python 3.10.21, R 4.3.2 for
parsers and R 4.4.3 for the scientific audit. Wheel testing used an isolated
`--target` installation; the user's analysis environment was not reinstalled.

Exact option spelling uses the documented
[Python argparse `allow_abbrev=False` behavior](https://docs.python.org/3.10/library/argparse.html#allow-abbrev)
through both Python and [R argparse](https://trevorld.r-universe.dev/argparse/doc/argparse.html).
Logs, the wheel, templates and change patch are retained in the v0.6.0 artifacts.
