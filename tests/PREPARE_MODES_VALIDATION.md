# Preparation modes validation

Validated on 2026-10-04 with Python 3.12.14, Polars 0.20.31, pandas 2.3.3,
NumPy 1.26.4, SciPy 1.13.1, installed bcftools and R 4.4.2.

## Results

- Full Python suite: **241 tests, OK; 6 optional integration tests skipped**.
  Two additional tests were added afterward; the final targeted suite below
  includes those additions and the final native runner.
- Final `test_prepare_modes`: **13 tests, all passed**, including actual
  bcftools extraction, Python munging, pinned GenomicSEM munging and reuse checks.
- R `test_failure_handling.R .`: all ten matrix acceptance cases passed,
  including trait order, duplicate pairs, intercept mapping and split/whole-genome
  GWAMA input validation.
- Before/after preparation comparison: default GPCA files, legacy raw files,
  QC reports, settings and status files were byte-identical on the same fixture
  (six files for default preparation; eight for legacy raw export).
- Documentation: all **20** preparation options documented, **6** preparation
  and reuse command examples accepted by their parsers, and **332** local links
  checked (external source links excluded). CLI help displays defaults and choices.

## Data and expected behavior

The integration fixture contains two traits (`002`, `001`, in that order), each
with seven VCF records. It covers duplicate IDs, a P value below the configured
floor, missing/invalid INFO, invalid frequency and zero SE. These are **synthetic
VCF fixtures**, not the user's original GWAS datasets.

Expected and observed results:

| Check | Expected and observed |
| --- | --- |
| Shared raw schema | Eleven columns: `SNP,CHR,BP,A1,A2,EAF,BETA,SE,P,N,INFO`; two retained SNPs per trait. |
| GPCA versus LDSC QC | Five GPCA rows versus two raw LDSC rows; LDSC INFO/HapMap QC does not remove GPCA rows. |
| N override | Shared LDSC N uses the manifest value; GPCA N remains the original NEF. LDSC-only preparation can omit NEF with an override. |
| Binary metadata | Explicit positive total N and both valid prevalences required; incomplete/invalid declarations fail. |
| Python munging | Output equals a direct call to the existing precision-preserving Python munger. At INFO 0.9, one usable SNP remains; reference-only placeholders are preserved. |
| Native munging | Output equals direct pinned `GenomicSEM::munge`; one SNP remains, with Z sign changed appropriately when reference A1 is swapped. |
| Paths with spaces | Preparation and native munging succeed; a path-free temporary trait name avoids upstream removal of spaces in output prefixes. |
| Reuse | Real checksum/N-convention validation accepts the generated manifest and sidecars. Both LDSC manifest readers preserve trait order. A modified checksum fails. Regression itself is mocked in this handoff test. |
| Failure handling | Missing SI and failed munging receive two attempts. Preparation failure publishes no tables; munging failure retains raw inputs/logs and publishes no `munged/`. |
| Existing outputs/options | Existing files remain unchanged. Irrelevant munging options fail before outputs are created. Legacy filter-provenance mismatch still fails. |

The Python success fixture uses BETA −0.08 so that the one retained SNP passes
the original LDSC median-effect check. An initial fixture with BETA −0.123457
correctly failed that upstream check; it was not disabled or altered.

## Native integration boundary

Python used the installed pinned CBIIT LDSC runtime through the existing
`ldsc_munge.py` wrapper. The test substitutes the launch command to use that
Python environment directly; Conda environment selection is covered by existing
runtime tests, not a new full Conda installation.

The complete GenomicSEM package is not installed locally. Native integration
used a temporary R package namespace containing the **unchanged pinned**
`munge.R`, `munge_main.R`, `utils.R` and `utils_sanitychecks.R`, importing their
installed dependencies. This exercises the actual munging functions, the new
R entry point and the Python worker/CLI path. It does **not** establish a full
GenomicSEM installation or validate a new native LDSC regression.

Reference sources:

- [Pinned Python munging](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py)
- [Pinned GenomicSEM public munging](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge.R)
- [Pinned GenomicSEM munging calculation](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge_main.R)

No original LDSC, GenomicSEM or GWAMA calculation source was changed. New modes
use existing extraction, QC, GPCA writing, parallel retry and munging helpers.
Original VCFs were unchanged in the tests. Original user VCF datasets were not
available in the inspected workspace, so original-data preparation and a full
installed-GenomicSEM run remain unverified.

Local audit logs: `/private/tmp/prepare-modes-2wl50z73/`.

## Reproduce

```bash
PYTHONPATH=src:tests python -m unittest discover -s tests -p 'test_*.py'
Rscript tests/test_failure_handling.R .
```

The Python integration test needs bcftools and the pinned Python LDSC runtime.
For native integration, set `LDSC_GPCA_TEST_RSCRIPT` to an Rscript executable
whose library contains GenomicSEM. Without that explicit setting the native
test is skipped. The temporary munging-only harness described above was used
for this validation; it is not distributed as a replacement GenomicSEM package.
