# shared/reporting.R: function bodies preserved from the original workflow.

write_pc1_reports <- function(eigenvalues, eigenvectors, trait_order, outdir,
                              pca_matrix_type, pc1_orientation,
                              matrix_eigen_tolerance = 1e-8,
                              negative_eigen_action = "warn") {
  raw_total <- sum(eigenvalues)
  positive <- pmax(eigenvalues, 0)
  positive_total <- sum(positive)
  raw_pct <- if (is.finite(raw_total) && raw_total > 0)
    100 * eigenvalues / raw_total else rep(NA_real_, length(eigenvalues))
  positive_pct <- if (is.finite(positive_total) && positive_total > 0)
    100 * positive / positive_total else rep(NA_real_, length(eigenvalues))
  negative <- eigenvalues < 0
  threshold <- matrix_eigen_tolerance * max(1, max(abs(eigenvalues)))
  all_pcs <- data.frame(PC = seq_along(eigenvalues), PCA_Matrix_Type = pca_matrix_type,
    Eigenvalue_Raw = eigenvalues,
    Variance_Explained_Percent_Raw = raw_pct,
    Cumulative_Variance_Explained_Percent_Raw = cumsum(raw_pct),
    Positive_Eigenvalue_Normalized_Percent = positive_pct,
    Cumulative_Positive_Eigenvalue_Normalized_Percent = cumsum(positive_pct),
    Negative_Eigenvalue = negative,
    Substantive_Negative_Eigenvalue = eigenvalues < -threshold,
    Conventional_Variance_Partition = !any(negative) && raw_total > 0,
    Interpretation = if (any(negative))
      "Indefinite estimated matrix: raw percentages are not a conventional variance partition; positive-normalized percentages are descriptive only" else
      if (raw_total <= 0) "Non-positive total: variance percentages undefined" else
      "Raw percentages describe the selected PCA matrix",
    stringsAsFactors = FALSE)

  eigenvalue_ok <- is.finite(eigenvalues[1L]) && eigenvalues[1L] > 0
  coefficient <- eigenvectors[,1L]
  length_ok <- length(coefficient) == length(trait_order)
  coefficient_ok <- all(is.finite(coefficient))
  loading <- if (eigenvalue_ok) coefficient * sqrt(eigenvalues[1L]) else rep(NA_real_, length(coefficient))
  finite_ok <- all(is.finite(loading))
  nonzero_ok <- finite_ok && any(abs(loading) > sqrt(.Machine$double.eps))
  multiplier <- if (finite_ok && pc1_orientation == "tutorial" && median(loading) < 0) -1 else 1
  coefficient <- coefficient * multiplier
  loading <- loading * multiplier
  reasons <- c(if (!eigenvalue_ok) "PC1 eigenvalue is not finite and positive",
    if (!length_ok) "PC1 coefficient count differs from trait count",
    if (!coefficient_ok) "PC1 coefficients are non-finite",
    if (!finite_ok) "PC1 loadings are non-finite",
    if (finite_ok && !nonzero_ok) "PC1 loadings are effectively zero")
  passed <- !length(reasons)
  qc <- data.frame(PCA_Matrix_Type = pca_matrix_type, Trait_Count = length(trait_order),
    PC1_Eigenvalue = eigenvalues[1L], PC1_Eigenvalue_Finite_Positive = eigenvalue_ok,
    Loading_Count_Matches_Traits = length_ok, Loadings_Finite = finite_ok,
    Loadings_Not_All_Zero = nonzero_ok, PC1_Sign_Method = pc1_orientation,
    PC1_Sign_Multiplier = multiplier,
    PC1_Loading_Median_Used = if (finite_ok) median(loading) else NA_real_,
    PC1_QC_Status = if (passed) "PASS" else "FAIL",
    PC1_QC_Reasons = if (passed) "PC1 numerical checks passed" else paste(reasons, collapse = "; "),
    Matrix_Negative_Eigenvalue_Count = sum(negative),
    Matrix_Substantive_Negative_Eigenvalue_Count = sum(eigenvalues < -threshold),
    Matrix_Policy_Allows_Analysis = !(negative_eigen_action == "error" && any(eigenvalues < -threshold)),
    Scope = "PC1 numerical QC only; not a substitute for input, CTI or GWAMA QC",
    stringsAsFactors = FALSE)
  contribution <- if (passed) 100 * coefficient^2 else rep(NA_real_, length(coefficient))
  proteins <- data.frame(Manifest_Order = seq_along(trait_order), Trait = trait_order,
    PCA_Matrix_Type = pca_matrix_type, PC1_Eigenvector_Coefficient = coefficient,
    PC1_Loading_Used = loading, PC1_Contribution_Percent = contribution,
    Contribution_Rank = if (passed) rank(-contribution, ties.method = "min") else rep(NA_integer_, length(contribution)),
    PC1_QC_Status = if (passed) "PASS" else "FAIL", stringsAsFactors = FALSE)
  reports <- list(GenomicPCA_All_PCs_Variance_Explained = all_pcs,
    GenomicPCA_PC1_QC = qc, GenomicPCA_PC1_Protein_Contributions = proteins)
  for (name in names(reports)) write.csv(reports[[name]],
    file.path(outdir, paste0(name, ".csv")), row.names = FALSE)
  invisible(reports)
}
