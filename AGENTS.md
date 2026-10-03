# Self-pair LDSC validation exception

The following narrow exception supersedes the parent instructions requiring all
LDSC standard errors to be positive and all numeric diagnostics to be finite:

- Accept a self-pair (`p1 == p2`) with finite `rg` within the existing tolerance
  of 1, `se == 0`, `z == +Inf`, and `p == 0`, provided all other checks pass.
- Preserve these values and report the exception in the QC audit. Do not impute
  a positive SE or replace the infinite Z with a finite value.
- Between-trait correlations, heritability and intercept estimates retain their
  existing validation rules. Missing values and negative SEs remain invalid.
- Apply the same rule during Python result collection and R trait selection and
  final validation. Consistent infinite duplicate values may agree; infinity,
  negative infinity and missing values must not be treated as interchangeable.
- Test both the accepted case and invalid neighboring cases. Existing valid
  inputs must retain identical matrices, PCA loadings and GWAMA results.
