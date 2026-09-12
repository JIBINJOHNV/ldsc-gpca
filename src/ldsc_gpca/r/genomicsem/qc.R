# genomicsem/qc.R: function bodies preserved from the original workflow.

qc_genomicsem <- function(x, traits, args, audit) {
  audit$traits <- traits
  x <- validate_genomicsem_structure(x)
  absent <- setdiff(traits, colnames(x$S))
  for (trait in absent) genomicsem_event(audit,
    if (args$allow_missing_traits) "warning" else "error", "Trait absent from LDSCoutput",
    trait, action = if (args$allow_missing_traits) "removed" else "stop")
  if (length(absent) && !args$allow_missing_traits) stop("Selected traits are absent from LDSCoutput.", call. = FALSE)
  selected <- traits[!traits %in% absent]
  hse <- genomicsem_se(x$V, colnames(x$S))
  rse <- genomicsem_se(x$V_Stand, colnames(x$S))
  invalid <- list()
  for (trait in selected) {
    reasons <- character()
    if (!is.finite(x$S[trait, trait]) || x$S[trait, trait] <= 0) reasons <- c(reasons, "Non-finite or non-positive heritability")
    if (!is.finite(hse[trait, trait]) || hse[trait, trait] <= 0) reasons <- c(reasons, "Non-finite or non-positive heritability SE")
    if (!is.finite(x$I[trait, trait]) || x$I[trait, trait] <= 0) reasons <- c(reasons, "Non-finite or non-positive self intercept")
    if (!is.finite(x$S_Stand[trait, trait])) reasons <- c(reasons, "Non-finite standardized diagonal")
    if (length(reasons)) invalid[[trait]] <- reasons
  }
  for (trait in names(invalid)) for (reason in invalid[[trait]]) genomicsem_event(audit,
    if (args$failed_ldsc_action == "error") "error" else "warning", reason, trait,
    action = if (args$failed_ldsc_action == "error") "stop" else "removed")
  if (length(invalid) && args$failed_ldsc_action == "error") stop("Selected traits failed heritability/intercept QC; see audit files.", call. = FALSE)
  selected <- selected[!selected %in% names(invalid)]
  # Standardized diagonal is a consistency check, not Python self-rg QC.
  bad_diag <- selected[abs(diag(x$S_Stand)[match(selected, colnames(x$S))] - 1) > 1e-8]
  for (trait in bad_diag) genomicsem_event(audit, "error", "S_Stand diagonal differs from 1 beyond 1e-8", trait, value = x$S_Stand[trait, trait], action = "stop")
  if (length(bad_diag)) stop("Inconsistent standardized diagonal.", call. = FALSE)
  # Identify invalid pairs without filtering by significance or rg magnitude.
  pairs <- if (length(selected) >= 2L) t(combn(selected, 2L)) else matrix(character(), 0, 2)
  bad_pairs <- list()
  for (i in seq_len(nrow(pairs))) {
    a <- pairs[i, 1]; b <- pairs[i, 2]
    reason <- character()
    for (name in c("S", "I", "S_Stand")) if (any(!is.finite(x[[name]][cbind(c(a,b), c(b,a))])))
      reason <- c(reason, paste("Missing/non-finite", name, "pair estimate"))
    if (!is.finite(rse[a,b]) || rse[a,b] <= 0) reason <- c(reason, "Non-finite or non-positive rg SE")
    if (!is.finite(hse[a,b]) || hse[a,b] <= 0) reason <- c(reason, "Non-finite or non-positive genetic covariance SE")
    if (length(reason)) bad_pairs[[length(bad_pairs) + 1L]] <- c(a, b, paste(reason, collapse = "; "))
  }
  for (pair in bad_pairs) genomicsem_event(audit,
    if (args$failed_ldsc_action == "error") "error" else "warning", pair[3], pair[1], pair[2],
    action = if (args$failed_ldsc_action == "error") "stop" else "resolve_by_trait_removal")
  if (length(bad_pairs) && args$failed_ldsc_action == "error") stop("Selected pairs failed QC; see audit files.", call. = FALSE)
  while (length(bad_pairs)) {
    counts <- vapply(selected, function(t) sum(vapply(bad_pairs, function(p) t %in% p[1:2], logical(1))), integer(1))
    drop <- selected[which.max(counts)] # deterministic tie: first in manifest
    affected <- Filter(function(p) drop %in% p[1:2], bad_pairs)
    reason <- paste(vapply(affected, function(p) paste(p[1], p[2], p[3], sep = ": "), character(1)), collapse = "; ")
    genomicsem_event(audit, "warning", paste("Removed to resolve invalid pairs:", reason), drop, action = "removed")
    selected <- selected[selected != drop]
    bad_pairs <- Filter(function(p) !drop %in% p[1:2], bad_pairs)
  }
  audit$retained <- selected
  if (length(selected) < 2L) stop("Fewer than two traits remain after QC.", call. = FALSE)
  x <- subset_genomicsem(x, selected)
  for (name in c("S", "I", "S_Stand", "V", "V_Stand")) x[[name]] <- validate_symmetric_matrix(x[[name]], name)
  # Verify standardized matrix provenance using the selected native covariance.
  expected <- x$S / sqrt(outer(diag(x$S), diag(x$S)))
  if (max(abs(expected - x$S_Stand)) > 1e-8 * max(1, max(abs(expected))))
    stop("S_Stand is inconsistent with standardized S; check RData provenance.", call. = FALSE)
  tryCatch(chol(x$I), error = function(e) stop("Selected CTI is not positive definite; no automatic repair was applied.", call. = FALSE))
  for (trait in selected) {
    z <- x$S[trait, trait] / hse[trait, trait]
    if (args$h2_z_warn_threshold > 0 && z < args$h2_z_warn_threshold)
      genomicsem_event(audit, "warning", paste("h2/SE below", args$h2_z_warn_threshold), trait, value = z)
  }
  out_of_range <- which(upper.tri(x$S_Stand) & abs(x$S_Stand) > 1, arr.ind = TRUE)
  for (i in seq_len(nrow(out_of_range))) {
    ij <- out_of_range[i, ]
    genomicsem_event(audit, if (args$rg_out_of_range_action == "error") "error" else "warning",
      "Genetic correlation outside [-1,1]; not clamped", selected[ij[1]], selected[ij[2]],
      value = x$S_Stand[ij[1],ij[2]], action = if (args$rg_out_of_range_action == "error") "stop" else "retained")
  }
  if (nrow(out_of_range) && args$rg_out_of_range_action == "error") stop("Out-of-range rg; see audit files.", call. = FALSE)
  x
}
