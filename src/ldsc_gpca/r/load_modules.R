# Load only shared code and the selected backend into the caller's environment.
load_gpca_modules <- function(root, backend = c("python_ldsc", "genomicsem"),
                              envir = parent.frame()) {
  backend <- match.arg(backend)
  assign("bundled_gwama_path", file.path(normalizePath(root, mustWork = TRUE),
    "vendor", "N_weighted_GWAMA.function.1_2_6.R"), envir = envir)
  suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(glue)
    library(parallel)
  })
  shared <- c("input.R", "matrix_validation.R", "reporting.R", "pca.R", "gwama.R")
  specific <- if (backend == "python_ldsc") {
    c("constants.R", "cli.R", "reader.R", "qc.R", "covariance.R", "reporting.R", "workflow.R")
  } else {
    c("cli.R", "reader.R", "qc.R", "reporting.R", "workflow.R")
  }
  for (module in c(file.path("shared", shared), file.path(backend, specific))) {
    path <- file.path(root, module)
    if (!file.exists(path)) stop("Missing R module: ", path,
      ". Keep the complete package R directory together.", call. = FALSE)
    sys.source(path, envir = envir)
  }
  invisible(NULL)
}
