root <- commandArgs(trailingOnly = TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
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
check <- function(x, order = traits) suppressWarnings(canonicalize_and_validate_ldsc(x, order))
baseline <- check(rows)
stopifnot(identical(rownames(baseline$S_Stand), traits),
          identical(baseline$S_Stand, t(baseline$S_Stand)),
          identical(baseline$I, t(baseline$I)),
          all(diag(baseline$I) == 1.05), baseline$I["C", "A"] == .1,
          baseline$I["C", "B"] == .2, baseline$I["A", "B"] == .3)
reversed <- copy(rows); reversed[2, `:=`(p1 = "A", p2 = "C")]
stopifnot(identical(check(reversed)$I, baseline$I),
          identical(check(reversed)$S_Stand, baseline$S_Stand),
          identical(check(rows[c(6,3,1,5,2,4)])$S_Stand, baseline$S_Stand))
dup <- rbind(rows, reversed[2])
stopifnot(identical(check(dup)$S_Stand, baseline$S_Stand))
dup[7, rg := .8]
fail(check(dup), "Conflicting duplicate")
fail(resolve_incomplete_ldsc_traits(dup, traits, "drop_traits"), "Conflicting duplicate")
fail(check(rows[-2]), "C <-> A")

fail(resolve_ldsc_trait_order(rows, c(traits, "UNKNOWN")), "UNKNOWN")
changed <- copy(rows); changed[, `:=`(p = .8, se = se*2, z = z/3,
                                      h2_obs_se = .05, h2_int_se = .02, gcov_int_se = .03)]
stopifnot(identical(compute_pc1(check(changed)$S_Stand, traits)$loadings,
                    compute_pc1(baseline$S_Stand, traits)$loadings))
bad <- copy(rows)
bad[p1 == "A" | p2 == "A", `:=`(rg = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_)]
bad[p2 == "A", h2_obs := -.03]
fail(resolve_incomplete_ldsc_traits(bad, traits, "error"), "missing/invalid")
retained <- suppressWarnings(resolve_incomplete_ldsc_traits(bad, traits, "drop_traits"))
stopifnot(identical(retained$trait_order, c("C", "B")),
          identical(retained$excluded_traits$Trait, "A"))
stopifnot(identical(check(retained$ldsc_rows, retained$trait_order)$S_Stand,
                    baseline$S_Stand[c("C", "B"), c("C", "B")]))
bad[p1 == "B" | p2 == "B", `:=`(rg = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_)]
fail(resolve_incomplete_ldsc_traits(bad, traits, "drop_traits"), "at least two")
bad <- copy(rows); bad[p1 != p2, gcov_int := 2]
fail(check(bad), "positive definite")

tmp <- tempfile(); dir.create(tmp)
manifest <- file.path(tmp, "traits.csv")
fwrite(data.table(note = "annotation", traitname = traits, gwas_name = "ignored"), manifest)
stopifnot(identical(read_trait_manifest(manifest), traits))
ldsc <- file.path(tmp, "ldsc.csv")
annotated <- copy(rows); annotated[, `:=`(extra = "text", h2 = "ignored metadata")]
setcolorder(annotated, rev(names(annotated)))
fwrite(annotated, ldsc)
stopifnot(identical(check(read_python_ldsc_selected(ldsc, traits))$S_Stand, baseline$S_Stand))

gwama <- data.frame(note = "annotation", P = .05, Z = 2, N = 1000,
                    EAF = .2, OA = "G", EA = "A", BP = 123, CHR = 1, SNPID = "rs1")
for (layout in c("split", "nosplit")) {
  chromosomes <- if (layout == "split") as.list(1:22) else list(NULL)
  for (chr in chromosomes) {
    paths <- build_gwama_input_paths(traits, tmp, layout, chr)
    for (path in paths) fwrite(gwama, path, sep = "\t")
  }
  preflight_gwama_inputs(traits, tmp, layout)
  dat <- read_gwama_data(traits, tmp, layout, if (layout == "split") 1 else NULL)
  stopifnot(identical(names(dat), traits), identical(names(dat[[1]]), gwama_required_columns),
            dat[[1]]$SNPID == "rs1", dat[[1]]$Z == 2)
}
alias <- gwama; names(alias)[names(alias) == "EA"] <- "A1"
names(alias)[names(alias) == "OA"] <- "A2"; names(alias)[names(alias) == "P"] <- "p"
stopifnot(all(gwama_required_columns %in% validate_gwama_columns(names(alias), "fixture", "A")))
both <- c(names(gwama), "A1", "A2", "p")
stopifnot(!anyDuplicated(validate_gwama_columns(both, "fixture", "A")))
fail(validate_gwama_columns(setdiff(names(gwama), "N"), "fixture", "A"), "missing required columns: N")
fail(validate_gwama_columns(c(names(gwama), "N"), "fixture", "A"), "ambiguous required columns: N")
unlink(file.path(tmp, paste0(traits[1], "_chr22_GenomicPCA_inputs.tsv")))
fail(preflight_gwama_inputs(traits, tmp, "split"), "Missing required GWAMA")
unlink(tmp, recursive = TRUE)
cat("R failure handling, all ten matrix acceptance cases, and extra-column checks passed\n")
