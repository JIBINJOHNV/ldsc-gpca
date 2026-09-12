# genomicsem/reader.R: function bodies preserved from the original workflow.

validate_genomicsem_structure <- function(x) {
  required <- c("S", "I", "S_Stand", "V", "V_Stand")
  if (!is.list(x) || !all(required %in% names(x)))
    stop("LDSCoutput must contain S, I, S_Stand, V and V_Stand.", call. = FALSE)
  traits <- colnames(x$S)
  if (is.null(traits) || !length(traits) || anyNA(traits) ||
      any(!nzchar(trimws(traits))) || anyDuplicated(traits))
    stop("S requires unique, non-empty trait column names.", call. = FALSE)
  k <- length(traits)
  for (name in required) {
    m <- as.matrix(x[[name]])
    n <- if (name %in% c("V", "V_Stand")) as.integer(k * (k + 1) / 2) else k
    if (!is.numeric(m) || !identical(dim(m), c(n, n)))
      stop(paste(name, "has invalid type or dimensions."), call. = FALSE)
    if (name %in% c("S", "I", "S_Stand")) {
      for (labels in dimnames(m)) {
        if (!is.null(labels) && !identical(labels, traits))
          stop(paste(name, "has conflicting trait labels/order."), call. = FALSE)
      }
      dimnames(m) <- list(traits, traits)
    }
    # Missing entries may be resolved by explicit trait removal later. However,
    # conflicting finite orientations are a structural error, not a drop rule.
    both <- is.finite(m) & is.finite(t(m))
    if (any(abs(m[both] - t(m)[both]) > 1e-8))
      stop(paste(name, "is asymmetric."), call. = FALSE)
    x[[name]] <- m
  }
  x
}

genomicsem_se <- function(v, traits) {
  m <- matrix(NA_real_, length(traits), length(traits), dimnames = list(traits, traits))
  d <- diag(v)
  d[!is.finite(d) | d < 0] <- NA_real_
  m[lower.tri(m, diag = TRUE)] <- sqrt(d)
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

# GenomicSEM V uses column-major lower-triangle (including diagonal) order.
# Map each selected/reordered pair back to that order; preserve all genuine
# cross-estimate covariances instead of reconstructing a diagonal V.

subset_genomicsem <- function(x, traits) {
  old <- colnames(x$S)
  index <- match(traits, old)
  pair_index <- matrix(0L, length(old), length(old))
  pair_index[lower.tri(pair_index, diag = TRUE)] <- seq_len(nrow(x$V))
  pair_index[upper.tri(pair_index)] <- t(pair_index)[upper.tri(pair_index)]
  selected <- pair_index[index, index, drop = FALSE]
  v_index <- selected[lower.tri(selected, diag = TRUE)]
  for (name in c("S", "I", "S_Stand")) x[[name]] <- x[[name]][index, index, drop = FALSE]
  for (name in c("V", "V_Stand")) x[[name]] <- x[[name]][v_index, v_index, drop = FALSE]
  if (!is.null(x$N)) {
    if (length(x$N) != nrow(pair_index) * (nrow(pair_index) + 1) / 2)
      stop("Native N has an unexpected length.", call. = FALSE)
    x$N <- matrix(as.numeric(x$N)[v_index], nrow = 1L)
  }
  x
}
