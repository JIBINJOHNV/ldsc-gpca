# Separate PC1/CTI validation evidence

Date: 2026-10-05. This change separates matrix assessment from the action taken
on a failed assessment, adds PC1-only execution and explicit exploratory CTI
selection, and preserves source estimates. The bundled GWAMA algorithm and
upstream LDSC estimators were not changed.

## Executed checks

The new CLI tests exercise real R subprocesses for Python LDSC CSV and a
GenomicSEM RData fixture. They cover:

- PC1 passing while CTI fails: strict validation fails; PC1-only succeeds
  without deleting traits or invoking SNP preparation/export.
- Explicit CTI selection, manifest order, identical retained estimates, both
  matrices rechecked after selection, and a fresh strict run of the retained CSV.
- Exhausted exclusion limits, conflicting/invalid options, and refusal to
  overwrite an existing matrix assessment.
- Covariance-mode selection retaining scale and source metadata.
- GenomicSEM selection consistently subsetting the full V/V_Stand dimensions.
- Shared numerical tests for positive, indefinite, singular and ill-conditioned
  CTI; invalid PC1; missing values/asymmetry; later negative genetic eigenvalues;
  and separation of PC1 computability from the optional whole-matrix error policy.
- Pipeline option routing and skipping GWAMA preparation/export in PC1-only mode.
- An invalid matrix blocking automatic VCF preparation before that work begins.

The test fixtures do not constitute a fresh GenomicSEM LDSC regression or
whole-genome GWAMA validation. Pipeline routing tests use mocked upstream stages.
Existing tests separately cover input validity, duplicate handling, matrix
mapping, worker retries, GWAMA and export behavior.

Commands (with the project's numerical environment active):

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m unittest discover -s tests -v
python scripts/update_cli_help.py --check
git diff --check
```

Final local regression result: **317 tests discovered; 309 passed, 8 skipped**.
The skipped R interface test was then run separately with an explicit Rscript
path and passed. The remaining seven opt-in native munging/LDSC integration
checks were not enabled. Additional CSV, TSV and whitespace export tests with
plain/gzipped input passed, including preservation of original numerical text.
Generated-help consistency and `git diff --check` passed. A wheel build
succeeded and included the new shared R module and both backend help files.

## Real-data CSV replay

Input: a supplied complete Python LDSC table with **1,095 traits and 1,199,025
ordered rows**, from the user's extracted LDSC analysis archive. It contains
observed-scale h2, self intercepts and cross-trait intercepts. No synthetic
GenomicSEM sampling covariance was constructed from this input.

Input SHA256:
`c58f0f151124389efc951b896496d603ba0f3657fd02e7a154f8aff99352b1ad`.

All runs used the actual package CLI, `--rg_normalization pair` and
`--pca_matrix correlation` defaults. A manifest was derived in source trait
order. The original input was preserved. No LDSC regression or SNP-level GWAMA
was run in this real-data replay.

| Run | Expected and observed outcome | Traits at final assessment | PC1 eigenvalue | CTI minimum eigenvalue |
| --- | --- | ---: | ---: | ---: |
| Original, `--validate_only` | Exit 1; PC1 PASS, CTI INDEFINITE; both reports saved | 1,095 | 194.568102299312 | -0.030361155606344 |
| Original, `--pc1_only` | Exit 0; all traits retained, CTI failure reported, no GWAMA | 1,095 | 194.568102299312 | -0.030361155606344 |
| Original, `--validate_only --cti_action explore_drop --max_cti_drop_fraction 0.01` | Exit 0; four exploratory exclusions, both checks pass | 1,091 | 193.570428758045 | 0.005801521585967 |
| Original, same exploration with fraction `0.002` | Exit 1 after the allowed two exclusions; no usable export | 1,093 | 194.014923220542 | -0.010354137658777 |
| Exported subset and matching manifest, `--validate_only` | Exit 0 without further exclusions | 1,091 | 193.570428758045 | 0.005801521585967 |

A separate NumPy reconstruction from the **original CSV** matched both exported
matrices within absolute/relative tolerance 1e-12, with identical trait order.
Sign-aligned PC1 loadings differed by at most **3.997e-15**. Independent Cholesky
passed on both retained-subset runs.

Byte-level subset hashing confirmed every retained CSV row, column and numerical
text was unchanged. The PC1-only export matched the original CSV checksum, and
the source checksum was unchanged after all runs. The failing original/limited
runs published no final `GenomicPCA_Retained_Traits.csv` or retained-results file.

The retained CTI condition number was approximately 25,047.54. The PC1 loading
RMS change after selection was approximately 0.000729805 and the maximum
absolute change approximately 0.002160881, after overall sign alignment.
These describe this exclusion sensitivity, not sampling SEs or proof of a
scientifically preferable phenotype. The other genetic eigenvalues remain
reported; they were not clipped in the exported matrix or used for automatic
trait removal.

The private input and large output matrices are not committed. Full commands,
logs, manifests, exports, independent-check results and the replay script are
saved with the user's input analysis folder in `cli_validation_20261005/`.

## Interpretation boundary

Passing establishes the tested numerical/input contracts and CLI behavior.
It does not establish why the original CTI was indefinite, that excluded traits
were faulty, a globally minimal exclusion set, main-analysis suitability, or
correct calibration of GWAMA BETA/SE/N_eff. No universal exclusion-fraction or
PC1-sensitivity threshold is inferred from this dataset.
