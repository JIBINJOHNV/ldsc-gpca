# Numerical threads and workers — v0.6.0

The current interface uses `--n_cores` in all managed commands. Standalone
preparation uses this option for its preparation workers. GPCA uses it for GWAMA;
its optional earlier VCF preparation stage has the separate `--prepare_workers`
setting. Only these spellings are accepted.

Ten tests in `test_threads_and_workers.py` passed on 2026-09-29, including real R
parsers for both backends, default/0/1/22/100 worker settings where supported,
negative/invalid/missing values, Python worker validation, independent GPCA
preparation workers, raw-script forwarding and subprocess return codes.
The tests also verify that numerical defaults are set before numerical imports,
existing environment values survive, and spawned Conda/R processes inherit them.

The unchanged startup helper supplies `1` only when these variables are unset:

- `OPENBLAS_NUM_THREADS`
- `OMP_NUM_THREADS`
- `MKL_NUM_THREADS`
- `VECLIB_MAXIMUM_THREADS`
- `NUMEXPR_NUM_THREADS`

Importing the CLI alone does not modify the environment or import numerical
libraries. Runtime commands set the defaults before lazy imports. Caller settings
are preserved, including values which may be invalid for the actual library.
Conda activation hooks or the libraries themselves can override these values.
An environment variable is not a runtime measurement of the active thread pool.

Per-mode defaults are unchanged: Python LDSC 5 workers, preparation 4, native
GenomicSEM munging 1, split GPCA 0/auto (up to 22 chromosomes). Whole-genome GPCA
and Windows GWAMA remain sequential. Native GenomicSEM LDSC remains serial;
its worker setting applies to munging. No production speedup is claimed.

For earlier threadpool measurements and full historical logs, see the v0.5.4
artifacts retained in the conversation. The current source no longer supports
alternative option spellings. See `INTERFACE_VALIDATION.md` for v0.6.0 checks.
