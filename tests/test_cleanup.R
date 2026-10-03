root <- commandArgs(TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
local({
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  traits <- c("002", "001", "1e3", "A")
  rows <- CJ(p1 = traits, p2 = traits)
  rows[, `:=`(rg = ifelse(p1 == p2, 1, .2), se = .1,
    h2_obs = .2, h2_obs_se = .02, h2_int = 1, h2_int_se = .01,
    gcov_int = ifelse(p1 == p2, 1, .01), gcov_int_se = .01)]
  rows[, `:=`(z = rg/se, p = 2*pnorm(-abs(rg/se)))]
  original <- copy(rows)
  baseline <- canonicalize_and_validate_ldsc(rows, traits)
  stopifnot(identical(rows, original))
  for (extension in c("csv", "csv.gz")) {
    path <- file.path(folder, paste0("input.", extension))
    connection <- if (extension == "csv.gz") gzfile(path, "wt") else file(path, "wt")
    write.csv(rows, connection, row.names = FALSE); close(connection)
    for (chunk_size in c(1L, 5L, 100L)) {
      selected <- read_python_ldsc_selected(path, traits, chunk_size)
      stopifnot(nrow(selected) == nrow(rows), identical(selected$p1, rows$p1),
                identical(selected$p2, rows$p2))
      resolved <- resolve_incomplete_ldsc_traits(selected, traits)
      checked <- canonicalize_and_validate_ldsc(resolved$ldsc_rows, traits)
      stopifnot(identical(checked$S_Stand, baseline$S_Stand), identical(checked$I, baseline$I),
                identical(evaluate_self_pair_qc(selected, traits), resolved$self_pair_qc))
    }
  }
  bad <- copy(rows); set(bad, j = "se", value = rep("invalid", nrow(bad)))
  fail(evaluate_self_pair_qc(bad, traits), "non-numeric")
  fail(canonicalize_and_validate_ldsc(bad, traits), "non-numeric")
  fail(evaluate_self_pair_qc(rows[, !"z"], traits), "missing required columns")

  traits <- c("A", "B")
  for (correlation in c(.3, -.3, 1.2)) for (method in c("tutorial", "as_computed")) {
    matrix <- matrix(c(1, correlation, correlation, 1), 2, dimnames = list(traits, traits))
    expected <- eigen(matrix, symmetric = TRUE)
    loading <- as.vector(expected$vectors %*% sqrt(diag(pmax(expected$values, 0)))[, 1])
    loading <- orient_pc1_loadings(loading, method)$loadings
    pc1 <- suppressWarnings(compute_pc1(matrix, traits, pc1_orientation = method, report_dir = folder))
    stopifnot(identical(pc1$loadings, loading))
    qc <- read.csv(file.path(folder, "GenomicPCA_PC1_QC.csv"))
    proteins <- read.csv(file.path(folder, "GenomicPCA_PC1_Protein_Contributions.csv"))
    stopifnot(qc$PC1_QC_Status == "PASS", qc$PC1_Sign_Multiplier == pc1$orientation$sign_multiplier,
      isTRUE(all.equal(proteins$PC1_Loading_Used, pc1$loadings, tolerance = 1e-14)))
  }
  zero <- matrix(0, 2, 2, dimnames = list(traits, traits))
  fail(compute_pc1(zero, traits, report_dir = folder), "finite positive")
  stopifnot(read.csv(file.path(folder, "GenomicPCA_PC1_QC.csv"))$PC1_QC_Status == "FAIL")
})
cat("Cleanup identifier, boundary-QC and PC1/report checks passed\n")
