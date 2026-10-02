# Trait-wide normalization validation

Validated on 2026-10-02 against baseline commit
`8ed2e6367d6de7a6efe10839fd0752178e2125f0`, in an isolated worktree.
Python 3.12.14, NumPy 1.26.4, pandas 2.3.3, Polars **0.20.31**, R 4.4.2.
Original Python LDSC, GenomicSEM and the bundled GWAMA estimator were not modified.

## Automated tests

The complete suite ran 138 tests: 134 passed and four optional integrations skipped.
A subsequent focused run passed nine tests, including the previously skipped R
interface test and a new compiler integration test. Across these runs, **136 unique
tests passed**; the three optional native-LDSC integration tests were not run because
no runtime was configured for that suite. Native LDSC regressions were not rerun.
After the final reader adjustment, the plain/gzip/in-place CSV test passed again.

Coverage includes a four-trait independent numerical oracle; signed and out-of-range
rg; both orientations; row order; unmodified original columns; automatic compiler
annotation; duplicate/conflicting h2; missing reverse/self-pairs; non-positive,
missing and non-finite estimates; liability/observed and mixed trait scales;
consistent scale conversion; compressed filenames containing brackets; malformed
headers; atomic in-place replacement; repeated annotation; and R matrix selection,
duplicate validation, unchanged defaults and SE/Z/P independence.
The existing R tests cover all ten matrix/input acceptance cases.

Run from the checkout using a supported Python environment:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Set `LDSC_GPCA_TEST_RSCRIPT` to an Rscript executable with argparse/data.table/glue
for the optional R interface test. Tests require no private data.

## Actual saved LDSC results and PCA

Inputs were the existing `ldsc_results.csv` files in
`python_gmenomicsem_comparison/python/{dataset}_python_ldsc/`, with trait order
from each saved `GenomicPCA_PC1_Weights_Used.csv`. All rows had positive observed
h2 and both orientations/self-pairs. Native comparison files came from
`genomicsem/{dataset}/gpca_gwama/`. Input SHA-256 hashes and every command/output
are retained with the local artifacts; original input files were unchanged.

| Check | C2 | NOISE |
| --- | ---: | ---: |
| Traits | 43 | 49 |
| Directed rows, including self-pairs | 1,849 | 2,401 |
| Rows successfully annotated | 1,849 | 2,401 |
| Default and explicit `pair` vs baseline matrices/PC1 | Exactly equal | Exactly equal |
| Original CSV cells vs enriched CSV | Exactly equal | Exactly equal |
| CTI and original SE matrices across modes | Exactly equal | Exactly equal |
| Maximum normalized matrix error vs independent formula | 1.12e-15 | 7.78e-16 |
| Maximum absolute rg change | 0.02463625 | 0.04720400 |
| Maximum absolute PC1 loading change | 0.00540557 | 0.00735974 |
| PC1 eigenvalue: original → trait-wide | 10.724812 → 10.699627 | 14.975089 → 14.848527 |
| Median annotation time, 12 calls, seconds | 0.0062 | 0.0163 |

Timing covers in-memory annotation only, including string-to-number conversion,
on these two tables. It excludes imports, disk I/O, R startup, PCA and GWAMA;
it is not an estimate of LDSC runtime savings or a large-scale benchmark.

The actual managed Python CLI was used for default, `pair`, and `trait_wide`
validation-only runs. The prior baseline R entry point independently established
unchanged defaults. A separate NumPy eigendecomposition reproduced new PC1
loadings within 1e-13. The C2 covariance mode was also checked against reconstructed
pair covariances, with self h2 on the diagonal.

## Comparison with saved GenomicSEM results

Correlation RMSE uses unique off-diagonal pairs (903 C2; 1,176 NOISE), excluding
self-pairs. PC1 comparisons align the sign of the entire vector once. Lower
RMSE here means closer agreement with this saved reference, not proven greater
statistical accuracy.

| Metric | C2 original | C2 trait-wide | NOISE original | NOISE trait-wide |
| --- | ---: | ---: | ---: | ---: |
| Correlation RMSE | 0.00522415 | 0.00484607 | 0.00521918 | 0.00383258 |
| PC1 loading RMSE | 0.00525206 | 0.00503005 | 0.00409820 | 0.00347933 |

Correlation disagreement decreased by **7.24% / 26.57%**; PC1 loading disagreement
by **4.23% / 15.10%**, respectively. Regression weighting and other differences
remain. Negative eigenvalues were reported, not repaired.

Two illustrative pairs with the largest absolute normalization changes:

| Dataset / pair | Original rg | Trait-wide rg | Saved GenomicSEM |
| --- | ---: | ---: | ---: |
| C2: FABP2 OID20200 – FABP6 OID20076 | 0.44036815 | 0.41573191 | 0.41573471 |
| NOISE: IDO1 OID30563 – IDO1 OID31474 | 1.00734225 | 0.96013825 | 0.96212663 |

The IDO1 result is a calculated rescaling, not clipping to [-1,1].

## GWAMA check on a reconstructed actual-data subset

The original per-trait nine-column GWAMA inputs were unavailable. Therefore,
**this is not a whole-genome rerun of those original inputs**. Existing cached,
allele-aligned munged Z/N arrays supplied 521 C2 and 259 NOISE complete-case SNPs.
These SNPs were previously selected from large differences/threshold crossings;
they are not a random or genome-representative sample.

For each mode, the actual R workflow and unchanged bundled GWAMA function ran
on reconstructed nine-column files. Alleles/positions and a shared EAF came from
the saved comparison metadata. Shared EAF does not affect the tested GWAMA Z
formula; reconstructed BETA/SE/N_eff were not validated or treated as original
analysis results. No INFO value or final harmonisation export was fabricated.

The output Z scores matched the independent calculation
`sum(sqrt(N) * loadings * Z) / sqrt(w' CTI w)` within 1e-11 in both modes.
Only PC1 loadings changed between modes; CTI and SNP inputs were held fixed.

| Subset | SNPs | Maximum absolute change in Z | RMS change in Z |
| --- | ---: | ---: | ---: |
| C2 | 521 | 0.00451589 | 0.00137674 |
| NOISE | 259 | 0.00380917 | 0.00125790 |

These checks demonstrate correct propagation of the selected correlation through
PCA and GWAMA. They do not establish whole-genome association accuracy or supply
new uncertainty estimates for normalized rg. Original LDSC SE/Z/P remain attached
to original rg only.

## Local reproducibility artifacts

The sibling directory `trait_wide_validation_20261002/` contains
`validate_real_data.py`, `summary.json`, `summary.csv`, input hashes, command logs,
enriched LDSC copies, baseline/default/pair/trait-wide PCA audits, per-pair and
per-trait comparison CSVs, and both reconstructed GWAMA subset outputs. Private
inputs and generated results are excluded from Git.

The mathematical reconstruction follows the
[pinned LDSC ratio definition](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py#L702-L713).
The unchanged GWAMA weighting follows the
[Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
