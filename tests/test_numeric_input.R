root <- commandArgs(TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
  error <- tryCatch({suppressWarnings(force(expr)); NULL}, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error)))
    stop("Expected ", pattern, "; got ", if (is.null(error)) "success" else conditionMessage(error))
}
out <- tempfile(); dir.create(out)
traits <- c("A", "B", "C")
rows <- CJ(p1 = traits, p2 = traits)
rows[, `:=`(rg = ifelse(p1 == p2, 1, .2), se = .1, h2_obs = .2,
  h2_obs_se = .03, h2_int = 1.05, h2_int_se = .01,
  gcov_int = ifelse(p1 == p2, 1.05, .1), gcov_int_se = .01)]
rows[, `:=`(z = rg/se, p = 2*pnorm(-abs(rg/se)))]
base <- canonicalize_and_validate_ldsc(rows, traits)
for (field in c("h2_obs", "h2_obs_se", "h2_int", "h2_int_se")) {
  bad <- copy(rows); set(bad, j = field, value = as.character(bad[[field]]))
  set(bad, i = 1L, j = field, value = "not-a-number")
  for (action in c("error", "drop_traits"))
    fail(resolve_incomplete_ldsc_traits(bad, traits, action), paste0(field, ".*input row 1"))
  file <- file.path(out, "bad.csv"); fwrite(bad, file)
  fail(read_python_ldsc_selected(file, traits, chunk_size = 2L), paste0(field, ".*file row 2"))
}
# Guard runs during normalization itself, including an unused supplied scale.
for (scale in c("observed", "liability", "mixed")) {
  bad <- copy(rows); bad[, `:=`(h2_liab = NA_real_, h2_liab_se = NA_real_)]
  bad[, h2_obs := as.character(h2_obs)]; bad[1, h2_obs := "broken"]
  fail(normalize_heritability_columns(bad, scale), "h2_obs.*input row 1")
}
# Missing estimates remain removable; scientific notation and whitespace parse.
for (token in c("", "NA", "NaN", "nan")) {
  missing <- copy(rows); missing[, h2_obs := as.character(h2_obs)]
  missing[1, h2_obs := token]
  fail(resolve_incomplete_ldsc_traits(missing, traits, "error"), "missing/invalid")
  kept <- resolve_incomplete_ldsc_traits(missing, traits, "drop_traits")
  stopifnot(identical(kept$trait_order, c("B", "C")))
}
valid <- copy(rows); set(valid, j = "h2_obs", value = rep(" 2e-1 ", nrow(valid)))
set(valid, j = "h2_obs_se", value = rep("3e-2", nrow(valid)))
stopifnot(identical(canonicalize_and_validate_ldsc(valid, traits)$S_Stand, base$S_Stand),
          identical(canonicalize_and_validate_ldsc(valid, traits)$I, base$I))
liab <- copy(valid); setnames(liab, c("h2_obs", "h2_obs_se"), c("h2_liab", "h2_liab_se"))
stopifnot(all(normalize_heritability_columns(liab)$H2_Scale == "liability"))
mixed <- copy(valid); mixed[, `:=`(h2_liab = NA_character_, h2_liab_se = NA_character_)]
mixed[p2 == "B", `:=`(h2_liab = h2_obs, h2_liab_se = h2_obs_se,
                       h2_obs = NA_character_, h2_obs_se = NA_character_)]
stopifnot(identical(normalize_heritability_columns(mixed, "mixed")$h2, rep(.2, 9L)))
# File row numbers include skipped traits and earlier chunks.
file <- file.path(out, "offset.csv")
skipped <- copy(rows[1]); skipped[, `:=`(p1 = "other", p2 = "other")]
bad <- rbind(skipped, rows); bad[, h2_obs := as.character(h2_obs)]; bad[6, h2_obs := "bad"]
fwrite(bad, file)
fail(read_python_ldsc_selected(file, traits, chunk_size = 2L), "h2_obs.*file row 7")
cat("Early numeric parsing and genuine missing-value policy tests passed\n")

unlink(out, recursive = TRUE)
