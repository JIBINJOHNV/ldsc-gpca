# Bundled modified GWAMA source

`N_weighted_GWAMA.function.1_2_6.R` is copied unchanged from the file supplied
by the project maintainer. Its header credits Hill Ip and Bart Baselmans,
version 1.2.6 (11 July 2019). The genomicPCA weighting modification follows
[Anna Fürtjes' tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html):
`w = t(t(sqrt(N))*h2)` and `sqrt_W <- w`, with signed PC1 loadings supplied
through the argument named `h2`.

Both GPCA backends use this file by default; `--source_path` overrides it.
Do not mistake the argument `h2` or the inherited log wording for actual
heritabilities. The source does not output INFO. Automatic export retains
its existing missing-INFO error policy; a scientifically justified explicit
`--gwama-output-info` override is needed with this source. No INFO value is
invented by bundling it. The inherited N_eff/BETA/SE formulas are unchanged.
