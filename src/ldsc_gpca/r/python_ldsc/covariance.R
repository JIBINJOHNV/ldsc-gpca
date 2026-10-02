# python_ldsc/covariance.R: function bodies preserved from the original workflow.

select_correlation_matrix <- function(ldsc, rows, trait_order, normalization = "pair",
                                      tolerance = 1e-3, epsilon = 1e-12,
                                      out_of_range_action = "warn") {
  if (normalization == "pair") return(ldsc$S_Stand[trait_order, trait_order, drop = FALSE])
  if (normalization != "trait_wide") stop("Unknown rg normalization.", call. = FALSE)
  rows <- rows[p1 %chin% trait_order & p2 %chin% trait_order]
  if (!all(c("rg_trait_wide", "normalization_status") %in% names(rows)))
    stop("trait_wide requires rg_trait_wide and normalization_status columns. Enrich the CSV with python -m ldsc_gpca.normalization first.", call. = FALSE)
  values <- suppressWarnings(as.numeric(rows$rg_trait_wide))
  bad <- !is.finite(values) | is.na(rows$normalization_status) | rows$normalization_status != "calculated"
  if (any(bad)) stop("Trait-wide normalization unavailable for: ",
    paste(head(paste(rows$p1[bad], rows$p2[bad], sep = " <-> "), 20L), collapse = ", "), call. = FALSE)
  i <- match(rows$p1, trait_order); j <- match(rows$p2, trait_order)
  pairs <- data.table(i = pmin(i, j), j = pmax(i, j), value = values)[,
    .(value = mean(value), spread = max(value) - min(value)), by = .(i, j)]
  if (any(pairs$spread > tolerance + epsilon))
    stop("Conflicting duplicate rg_trait_wide values.", call. = FALSE)
  if (any(abs(pairs$value[pairs$i == pairs$j] - 1) > epsilon))
    stop("Self-pair rg_trait_wide must equal 1.", call. = FALSE)
  matrix <- ldsc$S_Stand[trait_order, trait_order, drop = FALSE]
  matrix[cbind(pairs$i, pairs$j)] <- matrix[cbind(pairs$j, pairs$i)] <- pairs$value
  diag(matrix) <- 1
  matrix <- validate_symmetric_matrix(matrix, "trait-wide genetic correlation matrix")
  if (any(abs(matrix[upper.tri(matrix)]) > 1)) {
    message <- "Trait-wide rg outside [-1,1]; values were not clamped."
    if (out_of_range_action == "error") stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }
  matrix
}

construct_genetic_covariance_matrix <- function(correlation_matrix,
                                                heritability,
                                                trait_order,
                                                heritability_scale = c(
                                                  "liability", "observed", "mixed"
                                                ),
                                                tolerance = 1e-10,
                                                trait_scales = NULL) {
  heritability_scale <- match.arg(heritability_scale)
  if (is.null(trait_scales)) {
    if (heritability_scale == "mixed") stop("Mixed covariance requires explicit per-trait scale metadata.", call. = FALSE)
    trait_scales <- rep(heritability_scale, length(trait_order))
  }
  if (length(trait_scales) != length(trait_order) || anyNA(trait_scales) ||
      any(!trait_scales %in% c("observed", "liability")) ||
      (heritability_scale != "mixed" && any(trait_scales != heritability_scale)))
    stop("Invalid or incompatible covariance scale metadata.", call. = FALSE)
  correlation_matrix <- validate_symmetric_matrix(
    correlation_matrix,
    "genetic correlation matrix"
  )
  if (!identical(rownames(correlation_matrix), trait_order) ||
      !identical(colnames(correlation_matrix), trait_order)) {
    stop(
      "The phenotype order in the genetic correlation matrix differs from the manifest.",
      call. = FALSE
    )
  }
  if (!is.finite(tolerance) || tolerance <= 0) {
    stop("Covariance reconstruction tolerance must be positive.", call. = FALSE)
  }

  if (!is.null(names(heritability)) &&
      !identical(names(heritability), trait_order)) {
    stop(
      "The named heritability vector order differs from the manifest.",
      call. = FALSE
    )
  }
  heritability <- as.numeric(heritability)
  if (length(heritability) != length(trait_order)) {
    stop("The heritability vector length differs from the manifest.", call. = FALSE)
  }
  if (any(!is.finite(heritability)) || any(heritability <= 0)) {
    stop(
      glue(
        "Covariance PCA requires finite positive self-pair ",
        "heritability values on the explicitly selected scale(s)."
      ),
      call. = FALSE
    )
  }

  h2_scale <- sqrt(heritability)
  covariance_matrix <- correlation_matrix * tcrossprod(h2_scale)
  diag(covariance_matrix) <- heritability
  dimnames(covariance_matrix) <- list(trait_order, trait_order)
  covariance_matrix <- validate_symmetric_matrix(
    covariance_matrix,
    "derived genetic covariance matrix"
  )

  reconstructed_correlation <- covariance_matrix / tcrossprod(h2_scale)
  diag(reconstructed_correlation) <- 1
  maximum_reconstruction_error <- max(
    abs(reconstructed_correlation - correlation_matrix)
  )
  if (maximum_reconstruction_error > tolerance) {
    stop(
      glue(
        "Derived genetic covariance matrix failed its reconstruction check; ",
        "maximum correlation error was ",
        "{format(maximum_reconstruction_error, scientific = TRUE)}."
      ),
      call. = FALSE
    )
  }

  h2_column <- if (heritability_scale == "mixed") "trait-specific h2" else heritability_column_sets[[heritability_scale]][1L]
  attr(covariance_matrix, "source") <- paste0(
    "rg * sqrt(self ", h2_column, " trait 1 * self ",
    h2_column, " trait 2)"
  )
  attr(covariance_matrix, "heritability_scale") <- heritability_scale
  attr(covariance_matrix, "trait_scales") <- trait_scales
  attr(covariance_matrix, "maximum_reconstruction_error") <-
    maximum_reconstruction_error
  covariance_matrix
}

make_genetic_covariance_results <- function(correlation_results,
                                            covariance_matrix,
                                            trait_order,
                                            heritability_scale = c(
                                              "liability", "observed", "mixed"
                                            )) {
  heritability_scale <- match.arg(heritability_scale)
  correlation_results <- as.data.frame(correlation_results)
  trait_index_1 <- match(correlation_results$Trait_1, trait_order)
  trait_index_2 <- match(correlation_results$Trait_2, trait_order)
  if (anyNA(trait_index_1) || anyNA(trait_index_2)) {
    stop(
      "Genetic-correlation result order cannot be matched to covariance traits.",
      call. = FALSE
    )
  }

  covariance_results <- correlation_results[, c(
    "Trait_1", "Trait_2", "rg", "SE", "Z", "P"
  ), drop = FALSE]
  for (column in intersect(c("rg_used", "RG_Normalization"), names(correlation_results)))
    covariance_results[[column]] <- correlation_results[[column]]
  covariance_results$Heritability_Scale <- heritability_scale
  scales <- attr(covariance_matrix, "trait_scales", exact = TRUE)
  if (is.null(scales)) {
    if (heritability_scale == "mixed") stop("Missing per-trait covariance scales.", call. = FALSE)
    scales <- rep(heritability_scale, length(trait_order))
  }
  covariance_results$Heritability_Scale_1 <- scales[trait_index_1]
  covariance_results$Heritability_Scale_2 <- scales[trait_index_2]
  covariance_results$Heritability_1 <-
    diag(covariance_matrix)[trait_index_1]
  covariance_results$Heritability_2 <-
    diag(covariance_matrix)[trait_index_2]
  covariance_results$h2_liab_1 <- ifelse(covariance_results$Heritability_Scale_1 == "liability",
    covariance_results$Heritability_1, NA_real_)
  covariance_results$h2_liab_2 <- ifelse(covariance_results$Heritability_Scale_2 == "liability",
    covariance_results$Heritability_2, NA_real_)
  covariance_results$h2_obs_1 <- ifelse(covariance_results$Heritability_Scale_1 == "observed",
    covariance_results$Heritability_1, NA_real_)
  covariance_results$h2_obs_2 <- ifelse(covariance_results$Heritability_Scale_2 == "observed",
    covariance_results$Heritability_2, NA_real_)
  covariance_results$Genetic_Covariance_Derived <- covariance_matrix[
    cbind(trait_index_1, trait_index_2)
  ]
  h2_column <- if (heritability_scale == "mixed") "trait-specific h2" else heritability_column_sets[[heritability_scale]][1L]
  covariance_results$Derivation <- paste0(
    "selected rg * sqrt(self ", h2_column, "_1 * self ", h2_column, "_2)"
  )
  covariance_results
}

# Orient PC1 either with the Fuertjes tutorial convention (the default) or
# preserve the direction returned by eigen(). Both always operate on the whole
# loading vector; individual loadings are never selectively changed.
