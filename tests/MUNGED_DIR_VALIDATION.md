# Python LDSC supplied-directory protection

Validated on 2026-10-04.

## Behavior

Standalone Python LDSC rejects `--munged_dir` without `--ldsc_only` immediately
after argument parsing. It exits with code 2 before reading the manifest,
checking runtimes, creating outputs, archiving prior results or starting workers.
The error explains that fresh runs must omit `--munged_dir`.

Explicit reuse remains supported. `--ldsc_only` alone still uses
`<outdir>/ldsc_input`; fresh runs continue writing to that default directory.
GenomicSEM directory reuse and the existing pipeline guard are unchanged.
README, CLI help and the Python guides now describe the reuse-only contract.

## Regression tests

From the repository root:

```bash
PYTHONPATH=src:tests \
LDSC_GPCA_TEST_RSCRIPT=/Users/JJOHN41/miniconda3/bin/Rscript \
../.venv-ldsc-normalization/bin/python -m unittest \
  test_ldsc_dispatch test_pipeline test_restart test_interface_aliases -v
```

Observed: **74 tests passed, no skips**. Substitute equivalent Python/R paths
on other machines. The new rejection test failed against the original code and
passed after the guard was added.

The test covers 12 rejection combinations: direct Python entry, public default
entry and explicit Python backend selection; presence/absence of HapMap; and
new/existing output directories. Existing-output cases also supply `--restart`.
Every case verifies SHA-256 preservation of both supplied trait files and their
sidecars, unrelated files and existing results. No downstream stage is called.
Separate controls cover explicit/default-directory reuse and fresh-run routing.
Existing tests cover GenomicSEM backend dispatch, pipeline and restart behavior.

## Independent real extraction/munging check

The earlier disposable overwrite reproduction was repeated against the fix,
using real bcftools extraction and installed CBIIT munging on the same synthetic
VCF/HapMap fixture. A direct test launcher selected the existing environment;
Conda environment selection was not tested. Only later regression and result
compilation were stubbed.

| Case | Expected and observed |
| --- | --- |
| Supplied directory without `--ldsc_only` | Exit 2; no output directory; all supplied hashes unchanged; regression not reached |
| Explicit-directory reuse | Reuse gate passed; all supplied hashes unchanged |
| Default-directory reuse | Reuse gate passed; all supplied hashes unchanged |
| Fresh run | Output under `<outdir>/ldsc_input`; uncompressed munged file identical to the pre-fix run |

The fresh munged payload SHA-256 was
`282046a33d6d1529a4311a502750d42b5bef5e5d04ceb18f4ebe1e603e3ded50`
before and after the fix. Gzip-container bytes were not compared across fresh
runs because their timestamps can differ. Reused files were compared byte for
byte using SHA-256, including their gzip containers and metadata sidecars.

Local reproduction scripts, logs and all before/after hashes are retained in
`../munged_dir_verification_20261004/` (`verify_fix.py`, `verify_fix.log`,
`fixed_summary.json`). Disposable inputs/results are not committed.

`git diff --check` passed. No production inputs were used or modified, and these
checks do not claim a fresh end-to-end LDSC regression run.
