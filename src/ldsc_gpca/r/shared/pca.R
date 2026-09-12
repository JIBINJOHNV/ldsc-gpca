# shared/pca.R: function bodies preserved from the original workflow.

orient_pc1_loadings <- function(loadings,
                                method = c("tutorial", "as_computed")) {
  method <- match.arg(method)
  loadings <- as.numeric(loadings)
  if (length(loadings) == 0L || any(!is.finite(loadings))) {
    stop("PC1 loadings must be a non-empty finite numeric vector.", call. = FALSE)
  }

  median_before_orientation <- median(loadings)
  sign_multiplier <- if (
    method == "tutorial" && median_before_orientation < 0
  ) -1 else 1
  oriented_loadings <- loadings * sign_multiplier

  list(
    method = method,
    loadings = oriented_loadings,
    loadings_before_orientation = loadings,
    median_before_orientation = median_before_orientation,
    median_after_orientation = median(oriented_loadings),
    sign_multiplier = sign_multiplier,
    flip_applied = sign_multiplier == -1
  )
}

# Reports describe the selected PCA matrix. Negative later eigenvalues are
# matrix diagnostics, not evidence against any particular protein or PC1.

compute_pc1 <- function(pca_matrix, trait_order,
                        negative_eigen_action = c("warn", "error"),
                        matrix_eigen_tolerance = 1e-8,
                        pc1_orientation = c("tutorial", "as_computed"),
                        pca_matrix_type = c("correlation", "covariance"),
                        report_dir = NULL) {
  negative_eigen_action <- match.arg(negative_eigen_action)
  pc1_orientation <- match.arg(pc1_orientation)
  pca_matrix_type <- match.arg(pca_matrix_type)
  if (!is.finite(matrix_eigen_tolerance) || matrix_eigen_tolerance <= 0) {
    stop(
      "matrix_eigen_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }

  pca_matrix <- validate_symmetric_matrix(pca_matrix, "PCA input matrix")

  if (!identical(rownames(pca_matrix), trait_order) ||
      !identical(colnames(pca_matrix), trait_order)) {
    stop("The phenotype order in the PCA matrix differs from the manifest.", call. = FALSE)
  }

  eigen_selected <- eigen(pca_matrix, symmetric = TRUE)
  eigenvectors <- eigen_selected$vectors
  eigenvalues <- eigen_selected$values

  # Write failure status too when decomposition succeeds but PC1 is invalid.
  if (!is.null(report_dir)) write_pc1_reports(eigenvalues, eigenvectors, trait_order,
    report_dir, pca_matrix_type, pc1_orientation, matrix_eigen_tolerance,
    negative_eigen_action)

  if (!is.finite(eigenvalues[1L]) || eigenvalues[1L] <= 0) {
    stop(
      glue(
        "PC1 requires a finite positive first eigenvalue; observed ",
        "{eigenvalues[1L]}."
      ),
      call. = FALSE
    )
  }

  eigenvalues_after_pmax <- pmax(eigenvalues, 0)
  loadings_before_orientation <- as.vector(
    eigenvectors %*% sqrt(diag(eigenvalues_after_pmax))[, 1]
  )

  if (length(loadings_before_orientation) != length(trait_order) ||
      any(!is.finite(loadings_before_orientation))) {
    stop("PC1 loadings are invalid.", call. = FALSE)
  }
  if (all(abs(loadings_before_orientation) <= sqrt(.Machine$double.eps))) {
    stop("All PC1 loadings are effectively zero.", call. = FALSE)
  }

  orientation <- orient_pc1_loadings(
    loadings_before_orientation,
    method = pc1_orientation
  )
  loadings <- orientation$loadings
  eigenvectors_before_orientation <- eigenvectors
  eigenvectors[, 1L] <- eigenvectors[, 1L] * orientation$sign_multiplier

  negative <- eigenvalues < 0
  eigenvalue_scale <- max(1, max(abs(eigenvalues)))
  negative_eigenvalue_threshold <- matrix_eigen_tolerance * eigenvalue_scale
  substantive_negative <- eigenvalues < -negative_eigenvalue_threshold
  if (any(negative)) {
    eigen_message <- glue(
      "The {pca_matrix_type} PCA matrix has {sum(negative)} negative ",
      "eigenvalue(s), including {sum(substantive_negative)} below the ",
      "numerical threshold -{format(negative_eigenvalue_threshold, scientific = TRUE)}. ",
      "pmax(value, 0) applies only to loading calculation; the matrix itself ",
      "was not altered."
    )
    if (negative_eigen_action == "error" && any(substantive_negative)) {
      stop(eigen_message, call. = FALSE)
    }
    warning(eigen_message, call. = FALSE)
  }

  eigenvalue_results <- data.frame(
    PC = seq_along(eigenvalues),
    PCA_Matrix_Type = pca_matrix_type,
    Eigenvalue_Raw = eigenvalues,
    Eigenvalue_After_pmax = eigenvalues_after_pmax,
    Pmax_Guard_Applied = negative,
    Substantive_Negative_Eigenvalue = substantive_negative,
    Negative_Eigenvalue_Threshold = negative_eigenvalue_threshold,
    stringsAsFactors = FALSE
  )

  raw_total <- sum(eigenvalues)
  guarded_total <- sum(eigenvalues_after_pmax)
  loading_results <- data.frame(
    Trait = trait_order,
    PCA_Matrix_Type = pca_matrix_type,
    PC1_Eigenvector = eigenvectors[, 1L],
    PC1_Loading_Used = loadings,
    PC1_Standardised_Loading = if (pca_matrix_type == "correlation") {
      loadings
    } else {
      NA_real_
    },
    PC1_Correlation_Matrix_Loading = if (pca_matrix_type == "correlation") {
      loadings
    } else {
      NA_real_
    },
    PC1_Covariance_Matrix_Loading = if (pca_matrix_type == "covariance") {
      loadings
    } else {
      NA_real_
    },
    PC1_Eigenvector_Before_Orientation =
      eigenvectors_before_orientation[, 1L],
    PC1_Loading_Before_Orientation =
      orientation$loadings_before_orientation,
    PC1_Loading_Median_Before_Orientation =
      orientation$median_before_orientation,
    PC1_Loading_Median_After_Orientation =
      orientation$median_after_orientation,
    PC1_Sign_Method = orientation$method,
    PC1_Sign_Multiplier = orientation$sign_multiplier,
    PC1_Median_Flip_Applied = orientation$flip_applied,
    PC1_Eigenvalue_Raw = eigenvalues[1L],
    PC1_Eigenvalue_After_pmax = eigenvalues_after_pmax[1L],
    PC1_Variance_Explained = eigenvalues[1L] / raw_total,
    PC1_Variance_Explained_After_pmax =
      eigenvalues_after_pmax[1L] / guarded_total,
    Any_Nonpositive_Eigenvalue_Guarded = any(negative),
    Substantive_Negative_Eigenvalue_Count = sum(substantive_negative),
    Negative_Eigen_Action = negative_eigen_action,
    Matrix_Eigen_Tolerance = matrix_eigen_tolerance,
    stringsAsFactors = FALSE
  )

  list(
    pca_matrix_type = pca_matrix_type,
    loadings = loadings,
    eigenvectors = eigenvectors,
    eigenvalues = eigenvalues,
    eigenvalues_after_pmax = eigenvalues_after_pmax,
    orientation = orientation,
    loading_results = loading_results,
    eigenvalue_results = eigenvalue_results
  )
}
