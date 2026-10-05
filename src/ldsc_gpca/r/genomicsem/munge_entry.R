#!/usr/bin/env Rscript
# Internal single-trait runner; preparation owns scheduling and two-attempt retries.
args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("GenomicSEM", quietly = TRUE)) stop("Install the pinned GenomicSEM package first")
runtime <- packageDescription("GenomicSEM")
if (!identical(runtime$RemoteSha, "6b65ca5db39fdade08b0d811477be1cdd57b5039"))
  stop("GenomicSEM preparation requires the pinned package; run scripts/install_genomicsem.R")
if (identical(args, "--check")) quit(status = 0L)
if (length(args) != 5L) stop("Use ldsc-gpca prepare --help.")
# GenomicSEM munge strips spaces from trait.names; keep its temporary name path-free.
setwd(dirname(args[2]))
GenomicSEM::munge(files = args[1], hm3 = args[3], trait.names = "prepared", log.name = basename(args[2]), N = NA,
                 info.filter = as.numeric(args[4]), maf.filter = as.numeric(args[5]),
                 parallel = FALSE, cores = 1L, overwrite = TRUE)
if (!file.rename("prepared.sumstats.gz", paste0(basename(args[2]), ".sumstats.gz")))
  stop("Cannot publish GenomicSEM munged output")
write.table(data.frame(package_version = runtime$Version, source_commit = runtime$RemoteSha,
                       r_version = R.version.string), paste0(basename(args[2]), ".runtime.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
