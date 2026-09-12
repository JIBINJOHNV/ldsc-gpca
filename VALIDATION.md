# Local validation — 2026-09-12

Environment: macOS, Python 3.10.21, pandas 2.0.3, Polars 0.20.31,
setuptools 84.0.0. A temporary virtual environment reused the existing Python
dependencies; a completely fresh dependency installation was not tested.
Polars 0.20.31 was installed into the temporary environment because the host copy
lacked pip metadata. `pip check` then reported no broken requirements.

* Built a wheel and installed it into a temporary environment: succeeded.
* Built a source distribution including both compatibility scripts and tests: succeeded.
* Ran 25 unit/regression tests against the installed wheel: all passed.
* Exercised top-level help/version, LDSC help, and the legacy Python CLI: succeeded.
* Simulated GPCA execution: arguments including spaces remained intact; the R
  exit status was propagated; missing Rscript returned 127 with a clear message.
* Compared bundled and original R script bytes: identical.
* Tested process-spawn imports, normal/missing/invalid prevalence data, mixed-scale
  result compilation, incomplete/malformed results, filter provenance, quoted paths,
  simulated command failures, retry recovery and exhaustion: expected outcomes matched.
* P-value calculations and default analysis settings were preserved.

Limitations:

* No real Docker LDSC analysis was run. Command execution failures were simulated;
  the awk transformations and shell failure propagation were exercised locally.
* Real R help failed before analysis because the local `data.table` installation
  could not load `data_table.dylib`. No GPCA/GWAMA scientific run was performed.
* Python 3.11/3.12 and GitHub Actions are configured but have not been run locally.
* No repository license has been selected and no PyPI release has been made.
