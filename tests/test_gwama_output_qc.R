root <- commandArgs(TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
  error <- tryCatch({suppressWarnings(force(expr)); NULL}, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error)))
    stop("Expected ", pattern, "; got ", if (is.null(error)) "success" else conditionMessage(error))
}
out <- tempfile(); dir.create(out)
traits <- c("A", "B", "C")
# Actual bundled GWAMA; compare to an independent signed-weight equation.
inputs <- file.path(out, "inputs"); dir.create(inputs)
mk <- function(snp, z, n, ea = "A", oa = "G", af = .2) data.frame(
  SNPID = snp, CHR = 1L, BP = match(snp, c("shared", "partial", "swap", "underflow", "C_only")),
  EA = ea, OA = oa, EAF = af, N = n, Z = z, P = pchisq(z^2, 1, lower.tail = FALSE))
dat <- list(A = mk(c("shared", "partial", "swap", "underflow"), c(2.1, 2.1, 2.1, 1000.1), 1000),
            B = rbind(mk(c("shared", "underflow"), c(3.1, 1000.1), 2000),
                      mk("swap", 3.1, 2000, "G", "A", .8)),
            C = mk(c("shared", "partial", "underflow"), c(-1.1, -1.1, 1000.1), 500))
write_inputs <- function() for (trait in traits) fwrite(dat[[trait]],
  file.path(inputs, paste0(trait, "_GenomicPCA_inputs.tsv")), sep = "\t")
write_inputs()
CTI <- matrix(.1, 3, 3, dimnames = list(traits, traits)); diag(CTI) <- 1
correlation <- diag(3); dimnames(correlation) <- list(traits, traits)
loadings <- setNames(c(.8, -.6, .2), traits)
gwama <- load_modified_gwama(bundled_gwama_path)
folder <- file.path(out, "valid"); dir.create(folder)
name <- run_genomic_pca_gwama(trait_order = traits, CTI = CTI, pca_matrix = correlation,
  loadings = loadings, gpca_input_folder = inputs, outdir = folder,
  splitby_chr = "nosplit", my_GWAMA = gwama)
path <- file.path(folder, paste0(name, ".N_weighted_GWAMA.results.txt.gz"))
result <- fread(path)
for (snp in result$SNPID) {
  available <- which(vapply(dat, function(x) snp %in% x$SNPID, logical(1)))
  n <- vapply(dat[available], function(x) x$N[match(snp, x$SNPID)], numeric(1))
  z <- vapply(dat[available], function(x) {
    row <- x[match(snp, x$SNPID), ]
    row$Z * if (row$EA == "A") 1 else -1
  }, numeric(1))
  w <- sqrt(n)*loadings[available]
  expected <- sum(w*z)/sqrt(as.numeric(t(w) %*% CTI[available, available, drop = FALSE] %*% w))
  actual <- result[SNPID == snp]
  expected_n <- as.numeric(t(sqrt(n)) %*% solve(CTI[available, available, drop = FALSE]) %*% sqrt(n))
  stopifnot(abs(actual$N_eff - expected_n) < 1e-9)
  stopifnot(abs(actual$Z - expected) < 1e-10,
    abs(actual$BETA - expected/sqrt(sum(n)*2*.2*.8)) < 1e-10,
    abs(actual$SE - 1/sqrt(sum(n)*2*.2*.8)) < 1e-12,
    abs(actual$PVAL - pchisq(expected^2, 1, lower.tail = FALSE)) < 1e-12)
}
stopifnot(result[SNPID == "underflow", PVAL] == 0,
          result[SNPID == "partial", Direction] == "+?-",
          result[SNPID == "swap", Direction] == "+-?")
# Opposite signed weights with sum zero are valid; zero available weight is
# not inferred from their sum or from a zero association Z.
cancel <- setNames(c(1, -1/sqrt(2), 0), traits)
validate_gwama_output(path, dat, cancel, file.path(folder, "cancel"))

# A C-only SNP under a valid block-diagonal correlation has no available weight.
correlation[1, 2] <- correlation[2, 1] <- .5
loadings <- compute_pc1(correlation, traits)$loadings
stopifnot(loadings[3] == 0, all(abs(abs(loadings[1:2])-sqrt(.75)) < 1e-12))
dat$C <- rbind(dat$C, mk("C_only", 2.1, 1000)); write_inputs()
for (backend in c("python_ldsc", "genomicsem")) {
  env <- new.env(parent = globalenv())
  load_gpca_modules(file.path(root, "src/ldsc_gpca/r"), backend, env)
  folder <- file.path(out, backend); dir.create(folder)
  worker <- function(chr) env$safe_gwama_worker(chr, trait_order = traits, CTI = CTI,
    pca_matrix = correlation, loadings = loadings, gpca_input_folder = inputs,
    outdir = folder, splitby_chr = "nosplit", my_GWAMA = gwama)
  status <- make_run_status(run_parallel_workers(list(NULL), worker, 1L, folder))
  fwrite(status, file.path(folder, "GWAMA_Run_Status.csv"))
  stopifnot(!status$Success, status$Attempts == 2L, grepl("QC failed", status$Error))
  raw <- fread(file.path(folder, "GenomicPCA_PC1.N_weighted_GWAMA.results.txt.gz"))
  issues <- fread(file.path(folder, "GenomicPCA_PC1.GWAMA_QC_Issues.csv"))
  stopifnot(raw[SNPID == "C_only", is.na(Z) & is.na(PVAL) & N_eff == 1000],
    all(issues$SNPID == "C_only"), "zero_available_weight" %in% issues$Reason,
    all(c("Z_not_finite", "PVAL_not_finite_or_outside_0_1") %in% issues$Reason))
}
# Every exported numerical field is checked; P=0 remains valid.
for (field in c("BETA", "Z", "SE", "N_eff", "PVAL")) {
  values <- if (field == "PVAL") c(NA, Inf, -.1, 1.1) else if (field %in% c("SE", "N_eff"))
    c(NA, Inf, 0, -1) else c(NA, Inf, -Inf)
  for (value in values) {
    bad <- copy(result); set(bad, i = 1L, j = field, value = value)
    file <- file.path(out, "invalid.results.gz"); fwrite(bad, file, sep = "\t", compress = "gzip")
    fail(validate_gwama_output(file, dat, setNames(c(.8, -.6, .2), traits),
      file.path(out, "invalid")), "GWAMA output QC failed")
    stopifnot(file.exists(file))
  }
}
unlink(out, recursive = TRUE)
cat("Bundled GWAMA signed formula, partial availability, swaps, underflow, zero-weight failure and audits passed\n")
