# genomicsem/reporting.R: function bodies preserved from the original workflow.

new_genomicsem_audit <- function(outdir) {
  a <- new.env(parent = emptyenv())
  a$outdir <- outdir
  a$traits <- character()
  a$retained <- character()
  a$events <- data.frame(Severity = character(), Scope = character(),
    Trait = character(), Other_Trait = character(), Reason = character(),
    Value = character(), Action = character(), stringsAsFactors = FALSE)
  a
}

write_genomicsem_audit <- function(a) {
  write.csv(a$events, file.path(a$outdir, "GenomicSEM_QC_Events.csv"), row.names = FALSE)
  for (kind in c("Warnings", "Removed_Traits")) {
    rows <- if (kind == "Warnings") a$events$Severity == "warning" else a$events$Action == "removed"
    write.csv(a$events[rows, , drop = FALSE],
      file.path(a$outdir, paste0("GenomicSEM_QC_", kind, ".csv")), row.names = FALSE)
  }
  summarize <- function(trait, severity = NULL) {
    rows <- a$events$Trait == trait | a$events$Other_Trait == trait
    if (!is.null(severity)) rows <- rows & a$events$Severity == severity
    paste(unique(a$events$Reason[rows]), collapse = "; ")
  }
  summary <- data.frame(Manifest_Order = seq_along(a$traits), Trait = a$traits,
    Retained_For_Analysis = a$traits %in% a$retained,
    Removed = a$traits %in% a$events$Trait[a$events$Action == "removed"],
    Warning_Reasons = vapply(a$traits, summarize, character(1), severity = "warning"),
    All_QC_Reasons = vapply(a$traits, summarize, character(1)), stringsAsFactors = FALSE)
  write.csv(summary, file.path(a$outdir, "GenomicSEM_Trait_QC_Summary.csv"), row.names = FALSE)
}

genomicsem_event <- function(a, severity, reason, trait = "", other = "",
                            value = "", action = "retained", scope = NULL) {
  if (is.null(scope)) scope <- if (nzchar(other)) "pair" else if (nzchar(trait)) "trait" else "run"
  a$events <- rbind(a$events, data.frame(Severity = severity, Scope = scope,
    Trait = trait, Other_Trait = other, Reason = reason, Value = as.character(value),
    Action = action, stringsAsFactors = FALSE))
}

write_genomicsem_diagnostics <- function(x, pc1, manifest, args) {
  traits <- colnames(x$S)
  hse <- genomicsem_se(x$V, traits)
  rse <- genomicsem_se(x$V_Stand, traits)
  h <- data.frame(Trait = traits, h2 = diag(x$S), h2_SE = diag(hse))
  h$Z <- h$h2 / h$h2_SE; h$P <- 2 * pnorm(abs(h$Z), lower.tail = FALSE)
  ij <- which(lower.tri(x$S), arr.ind = TRUE)
  rg <- data.frame(Trait_1 = traits[ij[,2]], Trait_2 = traits[ij[,1]], rg = x$S_Stand[ij], SE = rse[ij])
  rg$Z <- rg$rg / rg$SE; rg$P <- 2 * pnorm(abs(rg$Z), lower.tail = FALSE)
  tables <- list(Global_Heritability_Estimates = h, Global_Genetic_Correlations = rg,
    GenomicPCA_PC1_Weights_Used = pc1$loading_results,
    GenomicPCA_Selected_Traits_Eigenvalues = pc1$eigenvalue_results,
    GenomicSEM_Retained_Manifest = manifest[match(traits, manifest$traitname), , drop = FALSE])
  for (name in names(tables)) write.csv(tables[[name]], file.path(args$outdir, paste0(name, ".csv")), row.names = FALSE)
  matrices <- list(GenomicPCA_CTI_Used = x$I, GenomicPCA_Correlation_Matrix_Used = x$S_Stand,
    GenomicPCA_Covariance_Matrix_Used = x$S, GenomicPCA_PCA_Matrix_Used = if (args$pca_matrix == "correlation") x$S_Stand else x$S,
    GenomicSEM_RG_SE_Matrix = rse)
  for (name in names(matrices)) write.csv(matrices[[name]], file.path(args$outdir, paste0(name, ".csv")), row.names = TRUE)
  global <- eigen(x$S_Stand, symmetric = TRUE)
  global_loadings <- global$vectors %*% diag(sqrt(pmax(global$values, 0)))
  dimnames(global_loadings) <- list(traits, paste0("PC", seq_along(traits)))
  write.csv(global_loadings, file.path(args$outdir, "Global_Standardised_Loadings.csv"))
  write.csv(data.frame(PC = seq_along(traits), Eigenvalue_Raw = global$values,
    Eigenvalue_After_pmax = pmax(global$values, 0)),
    file.path(args$outdir, "Global_PCA_Eigenvalues.csv"), row.names = FALSE)
  # Native matrices, including correctly subset/reordered full sampling covariance.
  LDSCoutput <- x
  save(LDSCoutput, file = file.path(args$outdir, "GenomicSEM_LDSC_Used.RData"))
  saveRDS(args, file.path(args$outdir, "GenomicSEM_Run_Settings.rds"))
  writeLines(capture.output(sessionInfo()), file.path(args$outdir, "GenomicSEM_SessionInfo.txt"))
}
