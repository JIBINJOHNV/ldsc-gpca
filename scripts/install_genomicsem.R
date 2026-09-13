#!/usr/bin/env Rscript
# Official source pinned to an immutable commit, not the moving master branch.
revision <- "6b65ca5db39fdade08b0d811477be1cdd57b5039"
prefix <- Sys.getenv("CONDA_PREFIX")
if (!nzchar(prefix)) stop("Activate the ldsc-gpca Conda/Mamba environment first.")
target <- file.path(prefix, "lib", "R", "library")
if (!dir.exists(target) || !startsWith(normalizePath(R.home()), paste0(normalizePath(prefix), "/"))) {
  stop("Rscript must belong to the active Conda environment.")
}
if (!requireNamespace("remotes", quietly = TRUE)) stop("Install r-remotes from environment.yml first.")
options(repos = c(CRAN = "https://cloud.r-project.org"))
remotes::install_github(
  paste0("GenomicSEM/GenomicSEM@", revision), lib = target,
  dependencies = NA, upgrade = "never", build_vignettes = FALSE
)
library(GenomicSEM, lib.loc = target)
installed <- packageDescription("GenomicSEM", lib.loc = target)
if (!identical(installed$RemoteSha, revision)) stop("Installed GenomicSEM revision does not match the requested commit.")
message("Verified GenomicSEM commit: ", installed$RemoteSha)
sessionInfo()
