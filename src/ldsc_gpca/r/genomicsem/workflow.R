# Native GenomicSEM orchestration and audited analysis execution.

genomicsem_main <- function(arguments = commandArgs(TRUE)) {
  parser <- genomicsem_parser()
  if (!length(arguments) || any(arguments %in% c("-h", "--help"))) {
    parser$print_help()
    return(invisible(NULL))
  }
  args <- tryCatch(parser$parse_args(arguments), error = function(e) {
    parser$print_help()
    stop(e)
  })
  for (name in c("h2_z_warn_threshold", "matrix_eigen_tolerance"))
    if (!is.finite(args[[name]]) || args[[name]] < 0 || (name == "matrix_eigen_tolerance" && args[[name]] == 0))
      stop(paste("Invalid", name), call. = FALSE)
  if (is.na(args$cores) || args$cores < 0L) stop("--cores must be non-negative.", call. = FALSE)
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(args$outdir)) stop("Cannot create output directory.", call. = FALSE)
  if (file.exists(file.path(args$outdir, "GenomicSEM_QC_Events.csv")))
    stop("Existing GenomicSEM run found; use a fresh --outdir to preserve audit files.", call. = FALSE)
  audit <- new_genomicsem_audit(args$outdir)
  on.exit(write_genomicsem_audit(audit), add = TRUE)
  withCallingHandlers(tryCatch({
    manifest <- fread(args$input, data.table = FALSE)
    if (!"traitname" %in% names(manifest)) stop("Manifest requires traitname.", call. = FALSE)
    traits <- as.character(manifest$traitname)
    if (length(traits) < 2L || anyNA(traits) || any(!nzchar(trimws(traits))) || anyDuplicated(traits))
      stop("Manifest requires at least two unique, non-empty trait names.", call. = FALSE)
    audit$traits <- traits
    # Isolated load: the RData cannot replace CLI arguments or helper functions.
    env <- new.env(parent = emptyenv())
    loaded <- load(args$ldsc_path, envir = env)
    if (!"LDSCoutput" %in% loaded) stop("RData must contain LDSCoutput.", call. = FALSE)
    x <- qc_genomicsem(env$LDSCoutput, traits, args, audit)
    traits <- colnames(x$S)
    matrix <- if (args$pca_matrix == "correlation") x$S_Stand else x$S
    pc1 <- compute_pc1(matrix, traits, negative_eigen_action = args$negative_eigen_action,
      matrix_eigen_tolerance = args$matrix_eigen_tolerance,
      pc1_orientation = args$pc1_orientation, pca_matrix_type = args$pca_matrix,
      report_dir = args$outdir)
    write_genomicsem_diagnostics(x, pc1, manifest, args)
    if (args$validate_only) {
      status <- data.frame(Chromosome = "not_run_validate_only", Success = NA,
        Output = NA_character_, Error = "GWAMA intentionally skipped by --validate_only")
    } else {
      if (is.null(args$gpca_input_folder)) stop("GWAMA requires --gpca_input_folder.", call. = FALSE)
      if (!file.exists(args$source_path)) stop("GWAMA function script not found: ", args$source_path, call. = FALSE)
      preflight_gwama_inputs(traits, args$gpca_input_folder, args$splitby_chr)
      fun <- load_modified_gwama(args$source_path)
      worker <- function(chr) {
        warnings <- character()
        result <- withCallingHandlers(safe_gwama_worker(chr, trait_order = traits, CTI = x$I,
          pca_matrix = matrix, loadings = pc1$loadings, gpca_input_folder = args$gpca_input_folder,
          outdir = args$outdir, splitby_chr = args$splitby_chr, my_GWAMA = fun, pca_matrix_type = args$pca_matrix),
          warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
        result$warnings <- warnings
        result
      }
      cores <- args$cores
      if (cores == 0L) {
        cores <- parallel::detectCores(logical = FALSE)
        cores <- if (is.na(cores)) 1L else max(1L, cores - 1L)
      }
      runs <- if (args$splitby_chr == "nosplit") list(worker(NULL)) else
        if (.Platform$OS.type == "windows" || cores == 1L) lapply(1:22, worker) else
          parallel::mclapply(1:22, worker, mc.cores = min(22L, cores))
      for (i in seq_along(runs)) {
        if (!is.list(runs[[i]])) runs[[i]] <- list(success = FALSE,
          chromosome = i, output = NA_character_, error = "Parallel worker did not return a valid result")
        label <- if (is.null(runs[[i]]$chromosome)) "whole_genome" else runs[[i]]$chromosome
        for (w in runs[[i]]$warnings) genomicsem_event(audit, "warning", paste("GWAMA", label, w), action = "reported")
        if (!isTRUE(runs[[i]]$success)) genomicsem_event(audit, "error",
          paste("GWAMA", label, runs[[i]]$error), action = "stop")
      }
      status <- make_run_status(runs)
    }
    write.csv(status, file.path(args$outdir, "GWAMA_Run_Status.csv"), row.names = FALSE)
    if (any(!status$Success, na.rm = TRUE)) stop("One or more GWAMA runs failed; see GWAMA_Run_Status.csv.", call. = FALSE)
    message("Completed. Trait warnings/removals and reasons: ", file.path(args$outdir, "GenomicSEM_QC_Events.csv"))
    invisible(status)
  }, error = function(e) {
    genomicsem_event(audit, "error", conditionMessage(e), action = "stop")
    stop(e)
  }), warning = function(w) {
    genomicsem_event(audit, "warning", conditionMessage(w), action = "reported")
  })
}
