#!/usr/bin/env Rscript
# Thin genomicsem entry point; source() loads functions without running analysis.
.gpca_entry_path <- if (sys.nframe() > 0L) sys.frame(1)$ofile else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.gpca_module_root <- dirname(normalizePath(.gpca_entry_path, mustWork = TRUE))
sys.source(file.path(.gpca_module_root, "load_modules.R"), envir = environment())
load_gpca_modules(.gpca_module_root, "genomicsem", envir = environment())
rm(.gpca_entry_path, .gpca_module_root)

if (sys.nframe() == 0L) genomicsem_main()
