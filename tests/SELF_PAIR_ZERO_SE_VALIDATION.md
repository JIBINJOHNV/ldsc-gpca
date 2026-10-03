# Self-correlation zero-SE validation

Validated on 2026-10-03 with Python 3.12.14 and R 4.4.2, against baseline
commit `a309e972489adc7358028053e9061e9c643573b1`.

## Rule and implementation

Only self-pairs with finite rg within the existing tolerance of 1 and the
combination `se=0, z=+Inf, p=0` receive an exception. Every other required
estimate still passes its existing checks. No numerical values are imputed.
This combination is consistent with the pinned LDSC
[`p_z_norm` implementation](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py#L29-L37).

Python CSV collection and legacy log recovery share the same predicate. R row
selection, final validation and self-pair auditing share numerical checks.
Duplicate checks distinguish finite, positive-infinite, negative-infinite and
missing values, including before optional trait removal.

The Python pair-status warning and R self-pair QC note identify the exception.
R preserves `Self_Z=Inf` and `Self_RG_SE=0`. Literal finite/positive audit flags
remain false for this case; new validity flags include the allowed exception.
Original Python LDSC, PCA and GWAMA calculation code is unchanged.

## Automated tests

From the repository root, with dependencies installed:

```bash
PYTHONPATH=src:tests LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript python -m unittest \
  test_results_csv test_results_precision test_failure_handling \
  test_ldsc_drop_traits test_normalization test_self_pair_zero_se -v
```

All **67 tests passed**, including R tests and all ten matrix/input acceptance
cases. Tests cover native CSV and legacy logs, strict/report/drop modes,
preserved native values, zero-SE audit round-trips, default/custom self-rg
tolerance, consistent infinite duplicates, and conflicting finite/infinite/
missing duplicates. Negative/missing SEs, wrong Z/P combinations, between-trait
zero SE, invalid h2/intercept SE, missing pairs and conflicting traits still
fail. Existing tiny positive SE and rounded-log precision recovery tests pass.

## Real-data regression checks

Both saved real LDSC tables were recompiled with the baseline and updated code.
Each was run through genomicPCA and the unmodified bundled GWAMA function using
both `pair` and `trait_wide` normalization.

| Dataset | Traits | Directed LDSC rows | GWAMA subset SNPs | Observed before/after difference |
| --- | ---: | ---: | ---: | --- |
| T3_C0_C1_C2 | 43 | 1,849 | 521 | Exactly zero |
| T3_C0_C1_NOISE | 49 | 2,401 | 259 | Exactly zero |

The expected and observed result was exact equality of the compiled LDSC table,
correlation and intercept matrices, PC1 loadings, eigenvalues, and every GWAMA
output column. Original self-pair SEs were positive (minimum approximately
`2.01e-8` and `4.08e-9`, respectively).

Separate stress-test copies set all self-pair SE/Z/P values to `0/Inf/0` while
preserving the remaining real estimates. All 92 self-pairs were retained and
audited. Matrices, loadings, eigenvalues and GWAMA outputs still matched the
baseline exactly in both normalization modes. These copies are test fixtures,
not newly fitted LDSC estimates. Original data files were verified unchanged by
SHA-256 hashes.

Local reproduction script, logs, complete comparisons and input hashes are in
`../self_pair_zero_se_validation_20261003/` relative to the repository root.
The GWAMA subsets reuse previously reconstructed real Z/N inputs with shared
reference EAF. They validate numerical invariance, not independent BETA/SE
calibration. This check did not rerun the underlying LDSC regressions or the
whole-genome GWAMA analysis.
