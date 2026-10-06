# Shared PC1/CTI assessment. Selection is an explicit exploratory policy.
assess_cti <- function(CTI) {
  CTI <- validate_symmetric_matrix(CTI, "CTI")
  eig <- eigen(CTI, symmetric = TRUE)
  values <- eig$values
  tolerance <- max(1e-8, 100 * nrow(CTI) * .Machine$double.eps * max(abs(values)))
  condition <- if (min(values) > 0) max(values) / min(values) else NA_real_
  cholesky <- tryCatch({chol(CTI); TRUE}, error = function(e) FALSE)
  pass <- cholesky && min(values) > tolerance && is.finite(condition) && condition <= 1e8
  status <- if (pass) "PASS" else if (min(values) < -tolerance) "INDEFINITE" else "SINGULAR_OR_ILL_CONDITIONED"
  list(eigen = eig, pass = pass, metrics = data.frame(
    CTI_Status = status, CTI_Minimum_Eigenvalue = min(values), CTI_Maximum_Eigenvalue = max(values),
    CTI_Negative_Eigenvalues = sum(values < 0), CTI_Condition_Number = condition,
    CTI_Cholesky_Pass = cholesky, CTI_Eigenvalue_Threshold = tolerance, CTI_Condition_Limit = 1e8))
}

empty_cti_exclusions <- function() data.frame(Step = integer(), Trait = character(),
  Score = double(), Minimum_Eigenvalue_Before = double(), Minimum_Eigenvalue_After = double(), Reason = character())

assess_pc1 <- function(matrix, traits, args, outdir) {
  result <- tryCatch(compute_pc1(matrix, traits,
    negative_eigen_action = args$negative_eigen_action,
    matrix_eigen_tolerance = args$matrix_eigen_tolerance,
    pc1_orientation = args$pc1_orientation, pca_matrix_type = args$pca_matrix,
    report_dir = outdir, defer_negative_policy = TRUE), error = function(e) e)
  if (inherits(result, "error")) return(list(pass = FALSE, result = NULL, error = conditionMessage(result)))
  policy_pass <- args$negative_eigen_action != "error" ||
    !any(result$eigenvalue_results$Substantive_Negative_Eigenvalue)
  list(pass = TRUE, policy_pass = policy_pass, result = result, error = "")
}

write_matrix_assessment <- function(outdir, cti, pc1, traits, step, pc1_only, previous) {
  eig <- if (pc1$pass) pc1$result$eigenvalues else c(NA_real_, NA_real_)
  row <- cbind(data.frame(Step = step, Traits = length(traits), PC1_Status = if (pc1$pass) "PASS" else "FAIL",
    PC1_Error = pc1$error, PC1_Eigenvalue = eig[1], PC2_Eigenvalue = eig[2],
    PC1_PC2_Gap = eig[1] - eig[2], PC1_Only = pc1_only,
    PCA_Matrix_Policy_Pass = isTRUE(pc1$policy_pass),
    GWAMA_Matrix_Eligible = pc1$pass && isTRUE(pc1$policy_pass) && cti$pass), cti$metrics)
  history <- rbind(previous, row)
  write.csv(history, file.path(outdir, "GenomicPCA_Matrix_Validation_History.csv"), row.names = FALSE)
  write.csv(row, file.path(outdir, "GenomicPCA_Matrix_Validation.csv"), row.names = FALSE)
  history
}

write_selection_sensitivity <- function(baseline, final, traits, outdir) {
  initial <- setNames(baseline$loadings, baseline$loading_results$Trait)[traits]
  aligned <- final$loadings
  if (sum(initial * aligned) < 0) aligned <- -aligned
  delta <- aligned - initial
  write.csv(data.frame(Trait = traits, Original_PC1_Loading = initial,
    Retained_PC1_Loading_Sign_Aligned = aligned, Change = delta),
    file.path(outdir, "GenomicPCA_PC1_Selection_Sensitivity.csv"), row.names = FALSE)
  write.csv(data.frame(Loading_RMSE = sqrt(mean(delta^2)), Maximum_Absolute_Loading_Change = max(abs(delta)),
    Original_PC1_Eigenvalue = baseline$eigenvalues[1], Retained_PC1_Eigenvalue = final$eigenvalues[1],
    Interpretation = "Descriptive sensitivity, not sampling uncertainty or proof of scientific validity"),
    file.path(outdir, "GenomicPCA_PC1_Selection_Sensitivity_Summary.csv"), row.names = FALSE)
}

validate_analysis_matrices <- function(pca_matrix, CTI, traits, args) {
  args <- validate_matrix_options(args)
  saveRDS(args, file.path(args$outdir, "GenomicPCA_Validation_Settings.rds"))
  original_traits <- traits
  exclusions <- empty_cti_exclusions()
  history <- NULL
  baseline <- NULL
  limit <- floor(length(traits) * args$max_cti_drop_fraction)
  repeat {
    # Assess CTI and PC1 separately before enforcing the requested policy.
    cti <- assess_cti(CTI)
    pc1 <- assess_pc1(pca_matrix, traits, args, args$outdir)
    history <- write_matrix_assessment(args$outdir, cti, pc1, traits, nrow(exclusions), args$pc1_only, history)
    if (nrow(exclusions)) exclusions$Minimum_Eigenvalue_After[nrow(exclusions)] <- cti$metrics$CTI_Minimum_Eigenvalue
    write.csv(exclusions, file.path(args$outdir, "GenomicPCA_CTI_Excluded_Traits.csv"), row.names = FALSE)
    write.csv(data.frame(traitname = traits), file.path(args$outdir, "GenomicPCA_Assessed_Traits.csv"), row.names = FALSE)
    write.csv(CTI, file.path(args$outdir, "GenomicPCA_CTI_Assessed.csv"))
    if (!pc1$pass) stop(pc1$error, call. = FALSE)
    if (!pc1$policy_pass) stop("PCA matrix has substantive negative eigenvalues: --negative_eigen_action error blocks this run even though PC1 is computable. No traits were removed for this warning.", call. = FALSE)
    if (is.null(baseline)) baseline <- pc1$result
    write_selection_sensitivity(baseline, pc1$result, traits, args$outdir)
    if (cti$pass || args$pc1_only) break
    if (args$cti_action == "error") stop(
      "The LDSC intercept matrix I is not positive definite or numerically stable for GWAMA. PC1 and CTI assessments were saved separately. Use --pc1_only for PC1 output without GWAMA, or explicitly configure exploratory CTI selection.", call. = FALSE)
    if (nrow(exclusions) >= limit || length(traits) <= 2L) stop(
      "CTI still fails within the permitted exclusion limit; no usable subset was published. See matrix history and exclusions.", call. = FALSE)
    target <- max(cti$metrics$CTI_Eigenvalue_Threshold, max(cti$eigen$values) / 1e8)
    bad <- cti$eigen$values <= target
    score <- as.vector((cti$eigen$vectors[, bad, drop = FALSE]^2) %*%
      pmax(target - cti$eigen$values[bad], .Machine$double.eps))
    drop <- which.max(score) # Deterministic ties follow the manifest order.
    exclusions <- rbind(exclusions, data.frame(Step = nrow(exclusions) + 1L, Trait = traits[drop],
      Score = score[drop], Minimum_Eigenvalue_Before = min(cti$eigen$values), Minimum_Eigenvalue_After = NA_real_,
      Reason = "Exploratory CTI eigenvector participation; not evidence of an erroneous trait"))
    traits <- traits[-drop]
    pca_matrix <- pca_matrix[traits, traits, drop = FALSE]
    CTI <- CTI[traits, traits, drop = FALSE]
  }
  if (!cti$pass) warning("CTI failed; --pc1_only permits PC1 output but GWAMA was not run.", call. = FALSE)
  write.csv(data.frame(traitname = traits), file.path(args$outdir, "GenomicPCA_Retained_Traits.csv"), row.names = FALSE)
  list(pc1 = pc1$result, traits = traits, CTI = CTI, pca_matrix = pca_matrix,
    excluded = exclusions, original_traits = original_traits, cti = cti)
}
