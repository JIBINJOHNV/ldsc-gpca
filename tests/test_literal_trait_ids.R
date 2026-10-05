args <- commandArgs(TRUE)
root <- args[1L]
group <- args[2L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
source(file.path(root, "src/ldsc_gpca/r/genomicsem/reader.R"))
source(file.path(root, "src/ldsc_gpca/r/genomicsem/ldsc_pipeline.R"))

fail <- function(expr, pattern) {
  error <- tryCatch({suppressWarnings(force(expr)); NULL}, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error)))
    stop("Expected error matching ", pattern, "; observed: ",
         if (is.null(error)) "success" else conditionMessage(error))
}

main <- function() {
  temporary <- tempfile("literal-trait-ids-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE))
  traits <- c("NaN", "001", "Trait_A", "NA", "nan")
  path <- file.path(temporary, "input.csv")

  if (group == "manifest") {
    for (quoted in c(FALSE, TRUE)) {
      manifest <- data.table(traitname = traits, note = "NA")
      fwrite(manifest, path, quote = quoted)
      stopifnot(identical(read_trait_manifest(path), traits),
                identical(read_gpca_manifest(path)$traitname, traits),
                identical(read_gpca_manifest(path)$note, manifest$note))
      for (invalid in c("", "   ")) {
        bad <- copy(manifest); bad[1, traitname := invalid]
        fwrite(bad, path, quote = quoted)
        fail(read_trait_manifest(path), "missing or empty")
      }
      duplicate <- rbind(manifest, manifest[1])
      fwrite(duplicate, path, quote = quoted)
      fail(read_trait_manifest(path), "unique")
    }
  }

  rows <- CJ(i = seq_along(traits), j = seq_along(traits))[i <= j]
  rows[, `:=`(p1 = traits[i], p2 = traits[j], rg = ifelse(i == j, 1, .2),
    se = .1, h2_obs = .2, h2_obs_se = .03, h2_int = 1.05, h2_int_se = .01,
    gcov_int = ifelse(i == j, 1.05, .1), gcov_int_se = .01)]
  rows[, `:=`(i = NULL, j = NULL, z = rg/se, p = 2*pnorm(-abs(rg/se)))]
  check <- function(x, scale = "auto") canonicalize_and_validate_ldsc(
    x, traits, heritability_scale = scale)
  select <- function(x, action) suppressWarnings(resolve_incomplete_ldsc_traits(
    x, traits, action, heritability_scale = "observed"))

  if (group == "ldsc") {
    baseline <- check(rows)
    baseline_weights <- compute_pc1(baseline$S_Stand, traits)$loadings
    for (scale in c("observed", "liability", "mixed")) {
      scaled <- copy(rows)
      if (scale == "liability") setnames(scaled, c("h2_obs", "h2_obs_se"), c("h2_liab", "h2_liab_se"))
      if (scale == "mixed") {
        scaled[, `:=`(h2_liab = NA_real_, h2_liab_se = NA_real_)]
        scaled[p2 %chin% traits[3:5], `:=`(h2_liab = h2_obs, h2_liab_se = h2_obs_se,
          h2_obs = NA_real_, h2_obs_se = NA_real_)]
      }
      for (quoted in c(FALSE, TRUE)) for (sep in c(",", "\t")) for (compressed in c(FALSE, TRUE)) {
        # Shuffled rows and small chunks exercise inference at every boundary.
        fwrite(scaled[.N:1], path, quote = quoted, sep = sep)
        input <- path
        if (compressed) {
          input <- paste0(path, ".gz")
          connection <- gzfile(input, "wt")
          writeLines(readLines(path), connection); close(connection)
        }
        parsed <- read_python_ldsc_selected(input, traits, chunk_size = 2L)
        stopifnot(nrow(parsed) == nrow(rows), identical(attr(parsed, "selected_traits_seen"), traits))
        for (allow in c(FALSE, TRUE))
          stopifnot(identical(resolve_ldsc_trait_order(parsed, traits, allow)$trait_order, traits))
        actual <- check(parsed)
        for (matrix in c("S_Stand", "I", "rg_se_matrix", "intercept_se_matrix"))
          stopifnot(identical(actual[[matrix]], baseline[[matrix]]))
        stopifnot(identical(actual$self_pair_qc$Self_H2, baseline$self_pair_qc$Self_H2),
                  identical(compute_pc1(actual$S_Stand, traits)$loadings, baseline_weights))
      }
    }

    # Missing tokens belong to estimates, never to p1/p2. Bad text stays fatal.
    for (quoted in c(FALSE, TRUE)) for (field in c("rg", "se", "h2_obs")) {
      for (token in c("NA", "NaN", "nan", "", "not-a-number")) {
        bad <- copy(rows)
        set(bad, j = field, value = as.character(bad[[field]]))
        set(bad, i = 1L, j = field, value = token)
        fwrite(bad, path, quote = quoted)
        parsed <- read_python_ldsc_selected(path, traits, heritability_scale = "observed")
        stopifnot(identical(attr(parsed, "selected_traits_seen"), traits), nrow(parsed) == nrow(rows))
        if (token == "not-a-number") {
          for (action in c("error", "drop_traits")) fail(select(parsed, action), "non-numeric")
        } else {
          fail(select(parsed, "error"), "missing/invalid")
          stopifnot(identical(select(parsed, "drop_traits")$trait_order, traits[-1L]))
        }
      }
    }
    for (quoted in c(FALSE, TRUE)) for (field in c("p1", "p2")) for (invalid in c("", "   ")) {
      bad <- copy(rows); set(bad, i = 1L, j = field, value = invalid)
      fwrite(bad, path, quote = quoted)
      fail(read_python_ldsc_selected(path, traits), paste0(field, ".*empty"))
    }
  }

  if (group == "genomicsem") {
    paths <- file.path(temporary, paste0(traits, ".sumstats"))
    for (file in paths) writeLines(c("SNP\tA1\tA2\tN\tZ", "rs1\tA\tG\t1000\t2"), file)
    manifest <- data.table(traitname = traits, source_file = paths,
      sampleprevalence = c("NA", "NaN", "nan", "", ".2"),
      populationprevalence = c("NA", "NaN", "nan", "", ".1"),
      N = c("NA", "NaN", "nan", "", "1000"))
    for (quoted in c(FALSE, TRUE)) {
      fwrite(manifest, path, quote = quoted)
      out <- file.path(temporary, paste0("genomicsem-", quoted)); dir.create(out)
      calls <- 0L
      estimator <- function(traits, sample.prev, population.prev, trait.names, stand, ...) {
        calls <<- calls + 1L
        stopifnot(identical(trait.names, manifest$traitname), identical(traits, paths),
                  identical(sample.prev, c(rep(NA_real_, 4), .2)),
                  identical(population.prev, c(rep(NA_real_, 4), .1)))
        n <- length(trait.names)
        S <- diag(.2, n); dimnames(S) <- list(trait.names, trait.names)
        I <- diag(n); dimnames(I) <- dimnames(S)
        m <- n*(n+1L)/2L
        list(S = S, I = I, S_Stand = I, V = diag(.001, m), V_Stand = diag(.001, m))
      }
      withCallingHandlers(run_genomicsem_ldsc(path, out, "ld", "ld", "", "existing", 1,
        .9, .01, 22, 200, 80, "error", ldsc_fun = estimator),
        warning = function(w) stop(conditionMessage(w)))
      selected <- read_gpca_manifest(file.path(out, "Selected_Traits.csv"))
      stopifnot(calls == 2L, identical(selected$traitname, traits),
                identical(selected$N, c(rep("NA", 4), "1000")))
      result <- new.env(); load(file.path(out, "genomicsem_LDSC.RData"), result)
      stopifnot(identical(colnames(result$LDSCoutput$S), traits))
      for (invalid in c("", "   ")) {
        bad <- copy(manifest); bad[1, traitname := invalid]
        fwrite(bad, path, quote = quoted)
        fail(run_genomicsem_ldsc(path, out, "ld", "ld", "", "existing", 1,
          .9, .01, 22, 200, 80, "error", ldsc_fun = estimator), "traitname.*empty")
        stopifnot(calls == 2L)
      }
    }
  }
  cat(group, "literal trait ID checks passed\n")
}
main()
