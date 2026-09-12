# python_ldsc/workflow.R: function bodies preserved from the original workflow.

gpsca_main <- function() {
  args <- validate_cli_paths(parse_command_line())
  manifest_trait_order <- read_trait_manifest(args$input)

  message(glue(
    "--- Reading Python LDSC estimates for {length(manifest_trait_order)} manifest traits ---"
  ))
  selected_ldsc <- read_python_ldsc_selected(
    args$python_ldsc,
    manifest_trait_order,
    args$ldsc_chunk_size,
    heritability_scale = args$heritability_scale,
    pca_matrix = args$pca_matrix
  )
  scale_audit <- data.frame(Trait = manifest_trait_order,
    Heritability_Scale = trait_heritability_scales(selected_ldsc, manifest_trait_order),
    Requested_Policy = args$heritability_scale, PCA_Matrix = args$pca_matrix)
  fwrite(scale_audit, file.path(args$outdir, "Python_LDSC_Heritability_Scales.csv"))
  source_rows_read <- attr(selected_ldsc, "source_rows_read")
  selected_heritability_scale <- attr(
    selected_ldsc,
    "heritability_scale",
    exact = TRUE
  )
  availability <- resolve_ldsc_trait_order(
    selected_ldsc,
    manifest_trait_order,
    args$allow_missing_traits
  )
  incomplete_resolution <- resolve_incomplete_ldsc_traits(
    selected_ldsc,
    availability$trait_order,
    args$failed_ldsc_action,
    self_rg_tolerance = args$self_rg_tolerance,
    comparison_epsilon = args$comparison_epsilon,
    h2_z_warn_threshold = args$h2_z_warn_threshold,
    heritability_scale = selected_heritability_scale
  )
  selected_ldsc <- incomplete_resolution$ldsc_rows
  trait_order <- incomplete_resolution$trait_order

  message("--- Validating pairwise LDSC coverage and constructing matrices ---")
  ldsc <- canonicalize_and_validate_ldsc(
    ldsc_rows = selected_ldsc,
    trait_order = trait_order,
    tolerance = args$duplicate_tolerance,
    source_rows_read = source_rows_read,
    duplicate_z_tolerance = args$duplicate_z_tolerance,
    self_rg_tolerance = args$self_rg_tolerance,
    comparison_epsilon = args$comparison_epsilon,
    h2_z_warn_threshold = args$h2_z_warn_threshold,
    z_consistency_tolerance = args$z_consistency_tolerance,
    z_consistency_action = args$z_consistency_action,
    rg_out_of_range_action = args$rg_out_of_range_action,
    heritability_scale = selected_heritability_scale
  )
  ldsc$validation_summary$Original_Manifest_Trait_Count <-
    length(manifest_trait_order)
  ldsc$validation_summary$Dropped_Missing_Trait_Count <-
    nrow(availability$missing_traits)
  ldsc$validation_summary$Dropped_Failed_LDSC_Trait_Count <-
    nrow(incomplete_resolution$excluded_traits)
  ldsc$validation_summary$Initial_Missing_Or_Invalid_Pair_Count <-
    incomplete_resolution$initial_failed_pair_count
  ldsc$validation_summary$Initial_Additional_Self_QC_Failure_Count <-
    incomplete_resolution$initial_self_qc_failure_count
  ldsc$validation_summary$Allow_Missing_Traits <- args$allow_missing_traits
  ldsc$validation_summary$Failed_LDSC_Action <- args$failed_ldsc_action
  ldsc$validation_summary$Heritability_Scale_Requested <-
    args$heritability_scale
  ldsc$validation_summary$Negative_Eigen_Action <- args$negative_eigen_action
  ldsc$validation_summary$Matrix_Eigen_Tolerance <-
    args$matrix_eigen_tolerance
  ldsc$validation_summary$Validate_Only <- args$validate_only

  # Manifest order is authoritative for every analysis object.
  order <- trait_order
  correlation_matrix <- ldsc$S_Stand[order, order, drop = FALSE]
  CTI <- ldsc$I[order, order, drop = FALSE]
  covariance_matrix <- NULL
  genetic_covariance_results <- NULL
  covariance_reconstruction_error <- NA_real_
  if (args$pca_matrix == "covariance") {
    heritability <- setNames(
      ldsc$heritability_results$Heritability,
      ldsc$heritability_results$Trait
    )[order]
    covariance_matrix <- construct_genetic_covariance_matrix(
      correlation_matrix,
      heritability,
      order,
      heritability_scale = ldsc$heritability_scale,
      trait_scales = ldsc$heritability_results$Heritability_Scale
    )
    covariance_reconstruction_error <- attr(
      covariance_matrix,
      "maximum_reconstruction_error",
      exact = TRUE
    )
    genetic_covariance_results <- make_genetic_covariance_results(
      ldsc$genetic_correlation_results,
      covariance_matrix,
      order,
      heritability_scale = ldsc$heritability_scale
    )
    pca_matrix <- covariance_matrix
  } else {
    pca_matrix <- correlation_matrix
  }

  pc1 <- compute_pc1(
    pca_matrix,
    order,
    negative_eigen_action = args$negative_eigen_action,
    matrix_eigen_tolerance = args$matrix_eigen_tolerance,
    pc1_orientation = args$pc1_orientation,
    pca_matrix_type = args$pca_matrix,
    report_dir = args$outdir
  )
  pc1$loading_results$Heritability_Scale_Used <-
    ldsc$heritability_scale
  pc1$eigenvalue_results$Heritability_Scale_Used <-
    ldsc$heritability_scale
  ldsc$validation_summary$PCA_Matrix_Type <- args$pca_matrix
  ldsc$validation_summary$Covariance_Matrix_Source <- if (
    args$pca_matrix == "covariance"
  ) {
    attr(covariance_matrix, "source", exact = TRUE)
  } else {
    NA_character_
  }
  ldsc$validation_summary$Covariance_Reconstruction_Maximum_Error <-
    covariance_reconstruction_error
  ldsc$validation_summary$PC1_Orientation_Method <-
    pc1$orientation$method
  ldsc$validation_summary$PC1_Loading_Median_Before_Orientation <-
    pc1$orientation$median_before_orientation
  ldsc$validation_summary$PC1_Loading_Median_After_Orientation <-
    pc1$orientation$median_after_orientation
  ldsc$validation_summary$PC1_Sign_Multiplier <-
    pc1$orientation$sign_multiplier
  ldsc$validation_summary$PC1_Median_Flip_Applied <-
    pc1$orientation$flip_applied

  self_pair_qc <- incomplete_resolution$self_pair_qc
  self_pair_qc$Retained_For_Analysis <-
    self_pair_qc$Trait %in% trait_order
  excluded_match <- match(
    self_pair_qc$Trait,
    incomplete_resolution$excluded_traits$Trait
  )
  self_pair_qc$Exclusion_Reason <-
    incomplete_resolution$excluded_traits$Reason[excluded_match]

  write_analysis_outputs(
    args$outdir,
    ldsc,
    pc1,
    correlation_matrix,
    pca_matrix,
    CTI,
    availability$missing_traits,
    incomplete_resolution$excluded_traits,
    self_pair_qc,
    covariance_matrix,
    genetic_covariance_results
  )

  if (isTRUE(args$validate_only)) {
    validation_status <- data.frame(
      Chromosome = "not_run_validate_only",
      PCA_Matrix_Type = args$pca_matrix,
      Heritability_Scale_Used = ldsc$heritability_scale,
      Success = NA,
      Output = NA_character_,
      Error = "GWAMA intentionally skipped by --validate_only",
      stringsAsFactors = FALSE
    )
    write.csv(
      validation_status,
      file.path(args$outdir, "GWAMA_Run_Status.csv"),
      row.names = FALSE
    )
    message(
      "--- Validation-only analysis complete. GWAMA was not run. Outputs: ",
      normalizePath(args$outdir, mustWork = FALSE),
      " ---"
    )
    return(invisible(validation_status))
  }

  message("--- Preflighting every GWAMA file before starting any GWAMA run ---")
  preflight_gwama_inputs(
    order,
    args$gpca_input_folder,
    args$splitby_chr
  )
  my_GWAMA <- load_modified_gwama(args$source_path)

  common_arguments <- list(
    trait_order = order,
    CTI = CTI,
    pca_matrix = pca_matrix,
    loadings = pc1$loadings,
    gpca_input_folder = args$gpca_input_folder,
    outdir = args$outdir,
    splitby_chr = args$splitby_chr,
    my_GWAMA = my_GWAMA,
    pca_matrix_type = args$pca_matrix
  )

  if (args$splitby_chr == "split") {
    physical_cores <- detectCores(logical = FALSE)
    if (is.na(physical_cores) || physical_cores < 1L) physical_cores <- 1L
    worker_count <- if (args$cores == 0L) {
      min(22L, max(1L, physical_cores - 1L))
    } else {
      min(22L, args$cores)
    }

    message(glue("--- Running chromosomes 1-22 with {worker_count} worker(s) ---"))
    worker_call <- function(chr) {
      do.call(safe_gwama_worker, c(list(chr = chr), common_arguments))
    }

    if (.Platform$OS.type == "windows" || worker_count == 1L) {
      gwama_results <- lapply(1:22, worker_call)
    } else {
      gwama_results <- mclapply(
        1:22,
        worker_call,
        mc.cores = worker_count,
        mc.preschedule = TRUE
      )
    }
  } else {
    message("--- Running whole-genome genomic PCA GWAMA ---")
    gwama_results <- list(do.call(
      safe_gwama_worker,
      c(list(chr = NULL), common_arguments)
    ))
  }

  run_status <- make_run_status(gwama_results)
  run_status$PCA_Matrix_Type <- args$pca_matrix
  run_status$Heritability_Scale_Used <- ldsc$heritability_scale
  write.csv(
    run_status,
    file.path(args$outdir, "GWAMA_Run_Status.csv"),
    row.names = FALSE
  )

  if (any(!run_status$Success)) {
    failed <- run_status[!run_status$Success, , drop = FALSE]
    stop(
      "One or more GWAMA runs failed:\n",
      paste0(failed$Chromosome, ": ", failed$Error, collapse = "\n"),
      call. = FALSE
    )
  }

  message(
    "--- Analysis complete. Audit files and GWAMA results are in: ",
    normalizePath(args$outdir, mustWork = FALSE),
    " ---"
  )
}
