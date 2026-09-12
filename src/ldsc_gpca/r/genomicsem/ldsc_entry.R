#!/usr/bin/env Rscript
# Internal runner; public CLI/help is ldsc-gpca genomicsem ldsc.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop("Use ldsc-gpca genomicsem ldsc --help.")
entry <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(entry), "reader.R"))
source(file.path(dirname(entry), "ldsc_pipeline.R"))
run_genomicsem_ldsc(args[1], args[2], args[3], args[4], args[5], args[6],
  as.integer(args[7]), as.numeric(args[8]), as.numeric(args[9]),
  as.integer(args[10]), as.integer(args[11]), as.numeric(args[12]), args[13])
