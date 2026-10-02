# Native covariance and pair heritability export validation

Validated on 2026-10-02 before committing, against parent
`5cc60b1ccb5cddddfd3a1c364e239dd26867095f`.

## Implementation

The existing in-process export hook copies `RG.gencov.tot`, `RG.hsq1.tot`, and
`RG.hsq2.tot` to `gcov_native_obs`, `h2_p1_pair_obs`, and `h2_p2_pair_obs`.
The original table and readable rendering are preserved. These fields always
retain the unconverted observed scale, including when the original table reports
liability h2. With effective sample size, this is the unconverted LDSC output,
not an assertion about population observed-scale interpretation.

Compilation retains the fields and calculates annotations in memory with Polars.
Native covariance is normalized using observed self-pair h2; older CSVs retain
the existing reconstruction fallback. `gcov_pair` retains its reported-scale
meaning. Default PCA/GWAMA uses original rg; `--rg_normalization trait_wide`
selects the additional column. Original SE/Z/P and CTI remain unchanged.
No original Python LDSC, GenomicSEM, or bundled GWAMA estimator was edited.
No extra fit or final-CSV reread is added to a fresh analysis. Existing restart
fingerprinting includes the exporter, so old checkpoints require rerunning to
obtain the new native fields.

## Full original datasets

The actual filtered summary statistics and matching LD scores from the saved
`python_gmenomicsem_comparison` analyses were used. Every directed pair was fitted
again, including both orientations and all self-pairs; no synthetic values were
substituted for native fit attributes.

| Check | C2 | NOISE |
| --- | ---: | ---: |
| Traits | 43 | 49 |
| Native directed fits | 1,849 | 2,401 |
| Unique pairs including self | 946 | 1,225 |
| SNPs per pair, minimum–maximum | 1,175,671–1,176,294 | 1,175,335–1,176,280 |
| Exported native fields equal fit attributes exactly | Pass | Pass |
| Readable table identical with/without export hook | Pass | Pass |
| Maximum refit vs saved rg difference | 2.96e-13 | 2.02e-13 |
| Maximum native vs reconstructed covariance difference | 1.67e-14 | 1.13e-14 |
| Maximum native vs legacy trait-wide rg difference | 6.10e-14 | 7.78e-14 |
| Maximum native vs legacy trait-wide PC1 loading difference | 3.71e-14 | 2.50e-14 |
| Original CSV cells preserved when native fields are attached | Exact | Exact |
| Default/explicit pair PCA matrices and PC1 vs previous version | Exact | Exact |
| CTI and diagnostic SE matrices after normalization | Unchanged | Unchanged |
| Native triangular input without reverse rows | Pass | Pass |
| Median annotation time, 12 repetitions | 4.48 ms | 5.29 ms |

The production compiler was exercised on all 92 native batch exports. Every
original/native numerical field survived compilation exactly. All 4,250 rows
received `normalization_status=calculated`. Actual Python CLI calls tested
default, explicit pair, trait-wide, triangular, and newly refitted input paths.
PC1 loadings were also checked independently using NumPy eigen decomposition.
Saved source CSVs were not modified.

The refits use the unchanged [pinned CBIIT regression source](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py)
with 200 jackknife blocks, `twostep=None`, M=1,173,569, the original per-trait
chi-square filter of 80, and identical reference/weight LD scores. The exact
upstream `_get_rg_table` function and allele-matching constants were loaded from
[sumstats.py](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/sumstats.py)
using AST selection, avoiding unrelated parsing dependencies. Regression,
jackknife, IRWLS, and sumstats source hashes were verified against live upstream
downloads; all four files matched byte for byte.

Four workers completed alignment preparation and 4,250 refits in 577.83 seconds.
This is validation cost, not added production overhead. Native refits used
Python 3.13.5, NumPy 2.2.6, SciPy 1.15.2, pandas 2.3.3; this differs from the
managed pinned LDSC environment. Package compilation/annotations/tests used
Python 3.12.14, NumPy 1.26.4, pandas 2.3.3, Polars 0.20.31. R was 4.4.2.

Refitting across numerical environments is not bit-identical to saved estimates.
Point estimates were checked with absolute tolerance 1e-10; SE/Z/P used absolute
1e-10 plus relative 1e-5. Extremely small self-pair SEs amplify tiny differences
in their Z ratios: maximum absolute Z differences were 23.69 and 752.82, against
saved Z values of approximately 49.8 million and 244 million. Maximum relative
Z differences were 6.73e-7 and 3.08e-6. This does not represent an export change:
within a fit, values are copied exactly. To isolate the feature from refitting
roundoff, default-output preservation was tested by attaching actual native
fields to the original saved cell text.

## GWAMA checks

The unmodified bundled GWAMA ran on the available reconstructed subsets of real
munged Z/N data: 521 C2 SNPs and 259 NOISE SNPs. Z scores matched an independent
`sum(sqrt(N)*loading*Z) / sqrt(w' CTI w)` calculation within 1e-11. Both modes'
Z/P outputs matched the preceding implementation within 1e-10. The original
whole-genome nine-column GWAMA inputs are unavailable; whole-genome GWAMA and
BETA/SE/N_eff accuracy were not validated. Shared reference EAF in the reconstructed
inputs is not treated as original per-trait EAF.

## Automated tests and reproduction

The full suite ran 145 tests: 141 passed, four skipped. A focused runtime run
then passed all ten tests, including the previously skipped real-bcftools test:
**142 unique tests passed**. Three optional full native-CLI environment integration
tests remain unrun. The full native-data check above separately validates the
unchanged estimators, native table builder, export hook, compiler and downstream
CLI. The nine exporter tests passed again after strengthening the liability test.

Tests include existing matrix/ordering/duplicate acceptance cases, complete and
triangular inputs, mixed scales, negative/zero covariance, non-positive/missing/
infinite h2, partial or malformed native schemas, inconsistent native rg or scale
conversion, conflicting duplicates, failed fit objects, missing self-pairs,
atomic write failures, unchanged logs, and one input read per batch.

```bash
PYTHONPATH=src LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
  python -m unittest discover -s tests -v
```

The local sibling folder `native_covariance_validation_20261002/` contains
`refit_native.py`, `validate_outputs.py`, source/input hash manifests,
`refit_provenance.json`, `refit_roundoff.json`, `summary.json`, native batch
exports, compiled/enriched CSVs, CLI matrix/PC1 audits and GWAMA subset outputs.
The scripts require the private original datasets at the recorded locations.
Large/private data and validation caches are not committed.

Original result CSV SHA-256:

- C2: `9beafb3579f852189be3b762bcfbe74aa3dfd5356c46bd3f75b552e3667fce52`
- NOISE: `5588a3bffa88977bfbd57cb49bcd90449761727173fb2fc5e9db4ce83bdefee5`

These checks establish correct export and agreement with the preceding
normalization on these datasets. They do not establish superiority over native
rg, remove regression-weight differences, or supply a sampling covariance matrix.
