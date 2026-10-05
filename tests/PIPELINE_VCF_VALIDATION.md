# VCF-only pipeline validation (v0.8.0)

This update makes `pipeline` accept GWAS summary-statistics VCFs as its only
study-data input. Standalone preparation, LDSC and GPCA commands retain their
existing input contracts. The Python/GenomicSEM estimators, mungers, VCF
converters, PCA implementation and bundled GWAMA function were not edited.

## Contract tests

`test_pipeline.py` checks full and validation-only VCF workflows for both
backends, automatic GenomicSEM VCF conversion, retained-trait ordering,
backend defaults, filter forwarding, optional raw-table export, required
references, malformed/missing inputs, invalid binary N/prevalences, and
failure-status preservation.

Removed input-mode flags are rejected before stages run or output directories
are created; supplied munged-file bytes remain unchanged. Populated raw/munged
manifest paths are rejected even when VCFs are also supplied. Blank metadata
columns do not select an alternative input route. Standalone parser and
interface tests confirm that their input modes remain available.

## Checks run on 2026-10-05

With Python 3.12, Polars, bcftools and Rscript available:

```bash
PYTHONPATH=src:tests python -m unittest discover -s tests -p 'test_*.py'
```

Observed: **305 tests run, 297 passed, 8 optional checks skipped**. These include
real bcftools extraction/preparation and R matrix/QC tests, alongside tests
that inject stage execution. The R interface test also passed in the earlier full-suite
run with an explicit Rscript setting.

An earlier full-suite attempt enabled installed-GenomicSEM integration tests
through `LDSC_GPCA_TEST_RSCRIPT`. Two such tests failed because that R runtime
could not find the pinned GenomicSEM package. They are not claimed as passing;
the final run left those optional integrations disabled. Fresh pinned Python
LDSC regression integration also requires a separately configured runtime.

Documentation checks accepted all **16 pipeline command examples**, checked
**56 Bash blocks** and resolved **221 local links** across the six updated user
guides. Syntax/argument checks do not execute
analysis examples against their placeholder paths.

## Before/after routing and input comparison

The baseline pipeline source was from commit `e3bb43101fc0cd4fb1502c80921f49594b0d2d74`.
A disposable three-trait fixture used two variants with positive and negative
effects, N values 1,000/2,000 and INFO 0.95/0.8. Real bcftools preparation was
used. LDSC fitting and downstream GPCA/GWAMA execution were injected.

| Route | Observed comparison |
| --- | --- |
| Python VCF, full run | Identical stage arguments and three GWAMA input tables |
| Python VCF, validation only | Identical stage arguments; no GWAMA tables requested |
| GenomicSEM VCF, full run | Identical stage arguments and six tables: three GWAMA inputs and three INFO-preserving LDSC raw inputs |
| GenomicSEM VCF, validation only | Identical stage arguments and three INFO-preserving LDSC raw inputs |

Output-directory paths were normalized and the same dataset identifier was used
when comparing arguments. All **12 generated-table byte comparisons passed**.
For GenomicSEM, the baseline explicitly supplied `--vcf_input`; the new interface
selects that same adapter automatically.

The old implicit GenomicSEM VCF fallback is intentionally not equivalent: it
exported no INFO. Its removal allows INFO filtering and can change retained
SNPs. That migration is documented in the pipeline guide.

This validation does **not** establish newly fitted estimates/SEs, a fresh
whole-genome run, runtime performance, or cross-backend equivalence. It verifies
the new input boundary and reuse of the existing supported VCF execution paths.
