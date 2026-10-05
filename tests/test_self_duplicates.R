root <- commandArgs(trailingOnly = TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
  error <- tryCatch({suppressWarnings(force(expr)); NULL}, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error), perl = TRUE))
    stop("Expected error matching ", pattern, "; observed: ",
         if (is.null(error)) "success" else conditionMessage(error))
}
traits <- c("C", "A", "B")
rows <- data.table(p1 = c("C", "C", "C", "A", "A", "B"),
                  p2 = c("C", "A", "B", "A", "B", "B"),
                  rg = c(1, .2, .3, 1, .4, 1), se = .1,
                  h2_obs = .2, h2_obs_se = .03, h2_liab = .3, h2_liab_se = .04,
                  h2_int = 1.05, h2_int_se = .01,
                  gcov_int = c(1.05, .1, .2, 1.05, .3, 1.05), gcov_int_se = .01)
rows[, `:=`(z = rg/se, p = 2*pnorm(-abs(rg/se)))]
check <- function(x, ...) suppressWarnings(canonicalize_and_validate_ldsc(
  x, traits, heritability_scale = "observed", ...))
select <- function(x, action, ...) suppressWarnings(resolve_incomplete_ldsc_traits(
  x, traits, action, heritability_scale = "observed", ...))
self_qc <- function(x) evaluate_self_pair_qc(x, traits, heritability_scale = "observed")
baseline <- check(rows)
fields <- c("h2_obs", "h2_obs_se", "h2_liab", "h2_liab_se", "h2_int", "h2_int_se")

for (field in fields) {
  for (value in c(.8, NA_real_, Inf, -Inf)) {
    extra <- copy(rows[1]); set(extra, j = field, value = value)
    bad <- rbind(rows, extra)
    pattern <- paste0("Conflicting duplicate[\\s\\S]*", field)
    fail(check(bad), pattern)
    fail(self_qc(bad), pattern)
    for (action in c("error", "drop_traits")) fail(select(bad, action), pattern)
  }
  for (difference in c(.001, .00100001)) {
    extra <- copy(rows[1]); set(extra, j = field, value = extra[[field]] + difference)
    candidate <- rbind(rows, extra)
    if (difference == .001) {
      invisible(check(candidate))
      for (action in c("error", "drop_traits"))
        stopifnot(identical(select(candidate, action)$trait_order, traits))
    } else {
      fail(check(candidate), "Conflicting duplicate")
    }
  }
}

# User-specified tolerances still apply, including after file normalization.
extra <- copy(rows[1]); extra[, h2_int := h2_int + .005]
stopifnot(all(check(rbind(rows, extra), tolerance = .01)$self_pair_qc$Self_QC_Pass))
for (action in c("error", "drop_traits"))
  stopifnot(identical(select(rbind(rows, extra), action, duplicate_tolerance = .01)$trait_order, traits))

# Large conflicts are structural failures, even if a negative source could
# otherwise remove the trait. A small mixture within tolerance is still invalid.
for (values in list(c(-.2, .8), c(-.0001, .0003))) {
  bad <- copy(rows); bad[1, h2_obs := values[1]]
  extra <- copy(bad[1]); extra[, h2_obs := values[2]]
  bad <- rbind(bad, extra)
  if (diff(values) > .001) {
    fail(check(bad), "Conflicting duplicate")
    for (action in c("error", "drop_traits")) fail(select(bad, action), "Conflicting duplicate")
  } else {
    qc <- self_qc(bad)
    stopifnot(qc$Self_H2[1] > 0, !qc$Self_H2_Positive[1], !qc$Self_QC_Pass[1])
    fail(check(bad), "self-pair")
    fail(select(bad, "error"), "missing/invalid")
    stopifnot(identical(select(bad, "drop_traits")$trait_order, c("A", "B")))
  }
}

# Invalid source SEs and non-finite estimates remain invalid even if all
# repeated copies agree. The unused alternative scale can be entirely empty.
for (field in fields) {
  for (value in if (endsWith(field, "_se")) c(0, -.0001, NA, Inf) else c(NA, Inf)) {
    bad <- copy(rows); set(bad, i = 1L, j = field, value = as.numeric(value))
    bad <- rbind(bad, bad[1])
    stopifnot(!self_qc(bad)$Self_QC_Pass[1])
    fail(check(bad), "non-finite|non-positive")
    stopifnot(identical(select(bad, "drop_traits")$trait_order, c("A", "B")))
  }
}
empty <- copy(rows); empty[, `:=`(h2_liab = NA_real_, h2_liab_se = NA_real_)]
stopifnot(identical(check(rbind(empty, empty[1]))$I, baseline$I))
# Even an unselected but supplied scale must not hide a negative source value.
bad <- copy(rows); bad[1, h2_liab := -.01]
bad <- rbind(bad, bad[1])
stopifnot(!self_qc(bad)$Self_H2_Positive[1])
fail(check(bad), "Non-positive self-pair heritability in a source row")
stopifnot(identical(select(bad, "drop_traits")$trait_order, c("A", "B")))

# Per-source self-rg validity cannot be rescued by averaging near its boundary.
bad <- copy(rows); bad[1, rg := 1.0101]
extra <- copy(bad[1]); extra[, rg := 1.0095]
bad <- rbind(bad, extra)
stopifnot(!self_qc(bad)$Self_RG_Within_Tolerance[1])
fail(check(bad), "Self-pair rg differs")
stopifnot(identical(select(bad, "drop_traits")$trait_order, c("A", "B")))

# Reverse off-diagonal h2 values and h2 intercepts belong to different p2 traits.
reverse <- copy(rows[2]); reverse[, `:=`(p1 = "A", p2 = "C")]
for (field in fields) set(reverse, j = field, value = .8)
consistent <- rbind(rows, rows[p1 == p2], reverse)
actual <- check(consistent[.N:1])
for (name in c("S_Stand", "I", "rg_se_matrix", "intercept_se_matrix"))
  stopifnot(identical(actual[[name]], baseline[[name]]))
stopifnot(identical(actual$self_pair_qc$Self_H2, baseline$self_pair_qc$Self_H2),
          identical(compute_pc1(actual$S_Stand, traits)$loadings,
                    compute_pc1(baseline$S_Stand, traits)$loadings))
covariance <- function(result) construct_genetic_covariance_matrix(
  result$S_Stand, result$self_pair_qc$Self_H2, traits, "observed")
stopifnot(identical(covariance(actual), covariance(baseline)),
          identical(compute_pc1(covariance(actual), traits)$loadings,
                    compute_pc1(covariance(baseline), traits)$loadings))
for (action in c("error", "drop_traits"))
  stopifnot(identical(select(consistent, action)$trait_order, traits))

# Reader normalization must preserve conflicts in the non-selected scale.
tmp <- tempfile(); dir.create(tmp)
extra <- copy(rows[1]); extra[, h2_liab := .8]
fwrite(rbind(rows, extra), file.path(tmp, "conflict.csv"))
read <- read_python_ldsc_selected(file.path(tmp, "conflict.csv"), traits,
                                 chunk_size = 2L, heritability_scale = "observed")
fail(check(read), "Conflicting duplicate[\\s\\S]*h2_liab")
for (action in c("error", "drop_traits")) fail(select(read, action), "Conflicting duplicate")

# Already normalized callers (without raw columns) still receive the same guard.
internal <- prepare_python_ldsc_rows(rows, traits, "observed")
internal[, (unlist(heritability_column_sets)) := NULL]
extra <- copy(internal[1]); extra[, h2 := .8]
fail(check(rbind(internal, extra)), "Conflicting duplicate")

# Exercise the bundled modified GWAMA on aligned, representative SNPs, using
# matrices/loadings from both the original table and accepted duplicates.
dat <- setNames(lapply(seq_along(traits), function(i) data.frame(
  SNPID = paste0("rs", 1:3), CHR = 1L, BP = 101:103, EA = "A", OA = "G",
  EAF = c(.2, .3, .4), N = 1000 + 100*i, Z = c(1, -2, 3) + i/10,
  P = 2*pnorm(-abs(c(1, -2, 3) + i/10)))), traits)
gwama <- load_modified_gwama(bundled_gwama_path)
for (mode in c("correlation", "covariance")) {
  outputs <- lapply(list(baseline, actual), function(result) {
    out <- tempfile(tmpdir = tmp); dir.create(out)
    matrix <- if (mode == "correlation") result$S_Stand else covariance(result)
    invisible(gwama(x = dat, cov_Z = result$I,
      h2 = compute_pc1(matrix, traits)$loadings, out = out, name = "parity",
      output_gz = FALSE, check_columns = FALSE))
    fread(file.path(out, "parity.N_weighted_GWAMA.results.txt"))
  })
  stopifnot(nrow(outputs[[1]]) == 3L, identical(outputs[[1]], outputs[[2]]))
}
unlink(tmp, recursive = TRUE)
cat("Self duplicate fields, policies, source validity, scales, tolerances, matrix/PC1/GWAMA parity passed\n")
