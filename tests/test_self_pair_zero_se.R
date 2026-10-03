root <- commandArgs(trailingOnly = TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern = "non-finite|non-positive|Self-pair") {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
traits <- c("C", "A", "B")
rows <- data.table(p1 = c("C", "C", "C", "A", "A", "B"),
                  p2 = c("C", "A", "B", "A", "B", "B"),
                  rg = c(1, .2, .3, 1, .4, 1), se = .1,
                  h2_obs = .2, h2_obs_se = .03,
                  h2_int = 1.05, h2_int_se = .01,
                  gcov_int = c(1.05, .1, .2, 1.05, .3, 1.05), gcov_int_se = .01)
rows[, `:=`(z = rg/se, p = 2*pnorm(-abs(rg/se)))]
check <- function(x, ...) suppressWarnings(canonicalize_and_validate_ldsc(x, traits, ...))
baseline <- check(rows)
zero <- copy(rows)
zero[p1 == p2, `:=`(se = 0, z = Inf, p = 0)]
actual <- check(zero)
stopifnot(identical(actual$S_Stand, baseline$S_Stand), identical(actual$I, baseline$I),
          identical(compute_pc1(actual$S_Stand, traits)$loadings,
                    compute_pc1(baseline$S_Stand, traits)$loadings),
          all(diag(actual$rg_se_matrix) == 0),
          all(actual$self_pair_qc$Self_QC_Pass),
          all(actual$self_pair_qc$Self_RG_Zero_SE_Exception),
          all(actual$self_pair_qc$Self_Z == Inf),
          all(!actual$self_pair_qc$Required_Numeric_Finite),
          all(!actual$self_pair_qc$Required_SE_Positive),
          all(actual$self_pair_qc$Required_Numeric_Valid),
          all(actual$self_pair_qc$Required_SE_Valid),
          all(nzchar(actual$self_pair_qc$Self_QC_Note)))
for (action in c("error", "drop_traits")) {
  selection <- resolve_incomplete_ldsc_traits(zero, traits, action)
  stopifnot(identical(selection$trait_order, traits), nrow(selection$excluded_traits) == 0)
}
reversed <- copy(zero); reversed[2, `:=`(p1 = "A", p2 = "C")]
stopifnot(identical(check(reversed[.N:1])$S_Stand, actual$S_Stand),
          identical(check(rbind(zero, zero[p1 == p2]))$S_Stand, actual$S_Stand))
for (z_value in c(-Inf, NA_real_, 100)) {
  duplicate <- copy(zero[1]); duplicate[, z := z_value]
  fail(resolve_incomplete_ldsc_traits(rbind(zero, duplicate), traits, "drop_traits"),
       "Conflicting duplicate")
}
for (field in c("se", "z", "p", "rg", "h2_obs_se", "h2_int_se", "gcov_int_se")) {
  values <- switch(field, se = c(-.1, NA, Inf), z = c(-Inf, NA, 100),
                   p = c(.05, NA, -1), rg = c(.98, 1.02, NA, Inf), c(0, NA))
  for (value in values) {
    bad <- copy(zero); set(bad, i = 1L, j = field, value = value)
    fail(check(bad))
    self_qc <- evaluate_self_pair_qc(bad, traits)
    stopifnot(!self_qc$Self_QC_Pass[1])
    selection <- suppressWarnings(resolve_incomplete_ldsc_traits(bad, traits, "drop_traits"))
    stopifnot(identical(selection$trait_order, c("A", "B")))
  }
}
bad <- copy(zero); bad[2, `:=`(se = 0, z = Inf, p = 0)]
fail(check(bad))
boundary <- copy(zero); boundary[1, rg := 1.01]
stopifnot(all(check(boundary)$self_pair_qc$Self_QC_Pass))
fail(check(boundary, self_rg_tolerance = .005))
selection <- suppressWarnings(resolve_incomplete_ldsc_traits(boundary, traits, "drop_traits",
                                                             self_rg_tolerance = .005))
stopifnot(identical(selection$trait_order, c("A", "B")))
boundary[1, rg := 1.02]
stopifnot(all(check(boundary, self_rg_tolerance = .03)$self_pair_qc$Self_QC_Pass),
          identical(resolve_incomplete_ldsc_traits(boundary, traits, "error",
                    self_rg_tolerance = .03)$trait_order, traits))
tmp <- tempfile(); dir.create(tmp)
fwrite(zero, file.path(tmp, "input.csv"))
roundtrip <- read_python_ldsc_selected(file.path(tmp, "input.csv"), traits)
invisible(resolve_incomplete_ldsc_traits(roundtrip, traits, "error", audit_outdir = tmp))
audit <- fread(file.path(tmp, "Python_LDSC_Self_Pair_QC.csv"))
stopifnot(all(audit$Self_Z == Inf), all(audit$Self_RG_Zero_SE_Exception), all(audit$Self_QC_Pass))
unlink(tmp, recursive = TRUE)
cat("Self-pair zero-SE acceptance, rejection, duplicates, selection, matrices and audit passed\n")
