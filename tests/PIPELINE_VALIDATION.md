# Pipeline integration validation — v0.7.0

Validated on 2026-10-03 against main commit
`f913277c113f3756c8de9517fb0036de9373a54f`.

## Scope

The existing unpublished pipeline draft was integrated into the current package.
It orchestrates preparation, Python LDSC or GenomicSEM, GPCA/GWAMA and export using
existing commands. Both backends now have a README contents entry, examples,
input contracts, grouped CLI help, defaults and choices.

Integration fixes include Python `--chisq_max auto`, forwarding
`--rg_normalization pair|trait_wide`, and handing Python's explicitly retained
traits to GPCA in original manifest order. GenomicSEM continues to use its own
retained-trait audit. Missing/invalid audits stop the workflow. Pipeline requires
a fresh output directory; restart remains a standalone LDSC operation.

The original local draft's 10 files were unchanged. All 35 checked existing R
and analysis files matched the baseline byte for byte, including GWAMA,
GenomicSEM orchestration, LDSC adapters, preparation and export. No estimator,
PCA formula or existing stage default was changed.

## Automated and interface checks

- Full suite: **217 tests run, 214 passed, three optional Python LDSC integration
  tests skipped** because no usable isolated child runtime was configured.
- All **23 pipeline tests passed**, including a real bcftools/preparation
  subprocess on a one-variant quantitative VCF. This is a preparation check,
  not a real LDSC regression.
- Native R two-pass handoff fixture passed: an injected estimator gives one
  trait negative heritability; both passes receive cutoff 80 and the second pass
  retains the other traits in their original order. The injected matrices are
  test fixtures only; production calls the native estimator.
- Coverage includes both backends, VCF/table/munged inputs, validation-only mode,
  fixed/automatic cutoffs, normalization forwarding, retained-trait handoff,
  leading-zero names, default/parser isolation, invalid arguments, missing
  files/metadata, duplicate traits, binary-input restrictions, subprocess failure,
  missing expected outputs, and preservation of existing output directories.
- README checks passed for **58 local links/anchors, 35 Bash blocks, 23 analysis
  examples and eight real R argument-parser examples**. All 61 Python-pipeline
  and 42 GenomicSEM-pipeline long options appear in the README and show defaults
  or required-input markers in help. Eight pipeline help routes passed.
- The built **0.7.0 wheel** contains the pipeline and shared formatter. Seven
  fresh-process help/version checks passed with analysis tools absent from PATH.

## Saved real-data replay

LDSC-stage output was supplied from the user's saved estimates. Pipeline then ran
real GPCA, the unchanged bundled GWAMA function, and export. Results were compared
with separate GPCA command runs using the same saved inputs and settings.

| Dataset | Traits | Available GWAMA SNP subset | Methods checked |
| --- | ---: | ---: | --- |
| T3_C0_C1_C2 | 43 | 521 | Python pair, Python trait-wide, GenomicSEM |
| T3_C0_C1_NOISE | 49 | 259 | Python pair, Python trait-wide, GenomicSEM |

All six runs matched exactly for six matrix/PC reports each: correlation, CTI,
PCA matrix, PC1 weights, selected eigenvalues and all-PC variance. Both exported
tables in every run were byte-identical after decompression. Shuffled retained
trait audits were restored to manifest order. Source input hashes were unchanged.

These are saved-estimate handoff checks, **not new LDSC regressions or complete
VCF-to-result runs**. GWAMA used reconstructed real Z/N subsets with shared
reference EAF. INFO=1 was an explicit test override, not an imputation-quality
estimate. The checks establish orchestration equivalence; they do not establish
BETA/SE calibration, whole-genome performance or cross-backend equivalence.
GenomicSEM was not installed in the validation R environment, and a usable
isolated Python LDSC runtime was not configured for fresh regressions.

## Reproduction

```bash
PYTHONDONTWRITEBYTECODE=1 POLARS_MAX_THREADS=4 \
  LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  LDSC_GPCA_TEST_MUNGE_SCRIPT=/path/to/pinned/munge_sumstats.py \
  PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src:tests python -m unittest test_pipeline -v
Rscript tests/test_pipeline_native.R "$PWD"
```

Validation used Python 3.12.14, Polars 0.20.31, pandas 2.3.3, NumPy 1.26.4 and
R 4.4.2 on macOS ARM64. This does not test every supported platform/version.
Local scripts, logs, input hashes, saved-data replay outputs, wheel and summaries
are in `../pipeline_integration_20261003/`; private data are not committed.
