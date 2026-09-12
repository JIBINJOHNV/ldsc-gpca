# shared/matrix_validation.R: function bodies preserved from the original workflow.

values_agree <- function(minimum, maximum, tolerance,
                         comparison_epsilon = 1e-12) {
  (maximum - minimum) <= tolerance + comparison_epsilon
}

validate_symmetric_matrix <- function(matrix_object, matrix_name,
                                      tolerance = 1e-8) {
  matrix_object <- as.matrix(matrix_object)

  if (nrow(matrix_object) != ncol(matrix_object)) {
    stop(glue("{matrix_name} must be square."), call. = FALSE)
  }
  if (any(!is.finite(matrix_object))) {
    stop(glue("{matrix_name} contains non-finite values."), call. = FALSE)
  }

  maximum_asymmetry <- max(abs(matrix_object - t(matrix_object)))
  if (maximum_asymmetry > tolerance) {
    stop(
      glue(
        "{matrix_name} is asymmetric; maximum absolute difference is ",
        "{format(maximum_asymmetry, scientific = TRUE)}."
      ),
      call. = FALSE
    )
  }

  (matrix_object + t(matrix_object)) / 2
}

assert_analysis_order <- function(dat, CTI, pca_matrix, loadings, trait_order) {
  if (!identical(names(dat), trait_order)) {
    stop("The phenotype order in dat differs from the manifest.", call. = FALSE)
  }
  if (!identical(rownames(CTI), trait_order) ||
      !identical(colnames(CTI), trait_order)) {
    stop("The phenotype order in CTI differs from names(dat).", call. = FALSE)
  }
  if (!identical(rownames(pca_matrix), trait_order) ||
      !identical(colnames(pca_matrix), trait_order)) {
    stop("The phenotype order in the PCA matrix differs from names(dat).", call. = FALSE)
  }
  if (length(loadings) != length(trait_order)) {
    stop("The PC1 loading order/length differs from names(dat).", call. = FALSE)
  }

  invisible(TRUE)
}
