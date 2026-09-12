# Python-LDSC argument definitions and help formatting.

requested_cli_color_mode <- function(arguments) {
  inline_mode <- grep("^--color=", arguments, value = TRUE)
  if (length(inline_mode) > 0L) {
    return(sub("^--color=", "", tail(inline_mode, 1L)))
  }

  option_position <- which(arguments == "--color")
  if (length(option_position) > 0L) {
    last_position <- tail(option_position, 1L)
    if (last_position < length(arguments)) {
      return(arguments[last_position + 1L])
    }
  }

  "auto"
}

cli_color_enabled <- function(mode = "auto") {
  if (identical(mode, "always")) return(TRUE)
  if (identical(mode, "never")) return(FALSE)

  no_color <- Sys.getenv("NO_COLOR", unset = NA_character_)
  if (!is.na(no_color) && nzchar(no_color)) return(FALSE)
  if (identical(tolower(Sys.getenv("TERM", unset = "")), "dumb")) {
    return(FALSE)
  }

  color_force <- Sys.getenv("CLICOLOR_FORCE", unset = "")
  if (nzchar(color_force) && color_force != "0") return(TRUE)

  isTRUE(isatty(stdout()))
}

colorize_cli_help <- function(help_text, enabled = TRUE) {
  if (!isTRUE(enabled)) return(paste0(help_text, "\n"))

  reset <- "\033[0m"
  bold_cyan <- "\033[1;36m"
  green <- "\033[32m"
  yellow <- "\033[33m"
  bold_yellow <- "\033[1;33m"
  bold_blue <- "\033[1;34m"
  lines <- strsplit(help_text, "\n", fixed = TRUE)[[1L]]
  group_headings <- c(
    "options:",
    "Required inputs:",
    "GWAMA inputs (not required with --validate_only):",
    "Trait and failed-result handling:",
    "QC comparison thresholds:",
    "PCA settings:",
    "Diagnostic warning/error controls:",
    "Execution and performance:"
  )

  for (line_index in seq_along(lines)) {
    line <- lines[line_index]
    if (grepl("^usage:", line)) {
      line <- sub("^usage:", paste0(bold_cyan, "usage:", reset), line)
    } else if (line %in% group_headings) {
      line <- paste0(bold_cyan, line, reset)
    } else if (identical(line, "IMPORTANT DEFAULTS")) {
      line <- paste0(bold_yellow, line, reset)
    } else if (line %in% c(
      "REQUIRED MANIFEST COLUMN",
      "REQUIRED PYTHON LDSC COLUMNS",
      "MATRIX DEFINITIONS",
      "PC1 SIGN CONVENTION",
      "VALIDATION-ONLY EXAMPLE"
    )) {
      line <- paste0(bold_blue, line, reset)
    } else if (grepl("^[[:space:]]{2}-{1,2}", line)) {
      line <- sub(
        "^([[:space:]]{2})(-{1,2}[^[:space:],]+(?:, --[^[:space:]]+)?)",
        paste0("\\1", green, "\\2", reset),
        line,
        perl = TRUE
      )
    }

    line <- gsub(
      "Default:",
      paste0(yellow, "Default:", reset),
      line,
      fixed = TRUE
    )
    lines[line_index] <- line
  }

  paste0(paste(lines, collapse = "\n"), "\n")
}

parse_command_line <- function() {
  raw_arguments <- commandArgs(trailingOnly = TRUE)
  parser <- ArgumentParser(
    prog = "ldsc-gpca gpca",
    add_help = FALSE,
    usage = paste(
      "%(prog)s --input MANIFEST.csv --python_ldsc LDSC.csv[.gz]",
      "--outdir DIRECTORY [options]"
    ),
    description = paste0(
      "Build genomic-PC1 weights and run the modified Fuertjes GWAMA method\n",
      "using a complete pairwise results table produced by Python ldsc.py."
    ),
    formatter_class = "argparse.RawDescriptionHelpFormatter",
    epilog = paste0(
      "INPUT FILE FORMATS\n",
      "  Manifest: comma-separated CSV with a header row.\n",
      "  Python LDSC: CSV, tab-separated TSV or whitespace-delimited text;\n",
      "    separator is detected from the header; gzip (.gz) is accepted.\n",
      "  GWAMA inputs: tab-separated TSV, exact ordered columns:\n",
      "    SNPID,CHR,BP,EA,OA,EAF,N,Z,P\n",
      "    A1/A2/p aliases are accepted for EA/OA/P.\n",
      "  Split files: {traitname}_chr{CHR}_GenomicPCA_inputs.tsv\n",
      "  Whole-genome files: {traitname}_GenomicPCA_inputs.tsv\n",
      "  --source_path: R source defining the already-modified GWAMA function.\n\n",
      "REQUIRED MANIFEST COLUMN\n",
      "  traitname    Unique, non-empty trait identifier. Each row selects one\n",
      "               trait, and row order defines all analysis matrices.\n",
      "  Direct R needs no other manifest columns. The package CLI additionally\n",
      "  requires vcf_files if --gpca_input_folder is omitted for preparation.\n\n",
      "REQUIRED PYTHON LDSC COLUMNS\n",
      "  Always required:\n",
      "    p1, p2, rg, se, z, p, h2_int, h2_int_se, gcov_int, gcov_int_se\n",
      "  Plus one complete heritability pair:\n",
      "    h2_liab, h2_liab_se  OR  h2_obs, h2_obs_se\n",
      "  Auto ignores empty pairs; ambiguous populated scales require a choice.\n",
      "  Mixed-scale covariance requires --heritability_scale mixed.\n\n",
      "MATRIX DEFINITIONS\n",
      "  Correlation matrix: diagonal = 1; off-diagonal = rg\n",
      "  Covariance matrix:  diagonal = selected-scale self h2; off-diagonal =\n",
      "                      rg * sqrt(h2_1 * h2_2)\n",
      "  GWAMA CTI matrix:  diagonal = h2_int; off-diagonal = gcov_int\n\n",
      "PC1 SIGN CONVENTION\n",
      "  Default --pc1_orientation tutorial: multiply the entire PC1 loading\n",
      "  vector by -1 when its median is negative. Use 'as_computed' only to\n",
      "  retain R's arbitrary eigenvector direction for legacy comparison.\n\n",
      "IMPORTANT DEFAULTS\n",
      "  Extra LDSC traits are ignored. Missing or invalid selected values stop\n",
      "  the run. No LDSC value is imputed or silently clamped. Use\n",
      "  --failed_ldsc_action drop_traits only when automatic removal is wanted.\n\n",
      "VALIDATION-ONLY EXAMPLE\n",
      "  Rscript gpsca_gwama_python_ldsc.r \\\n",
      "    --input selected_traits.csv \\\n",
      "    --python_ldsc all_pairwise_ldsc.csv.gz \\\n",
      "    --outdir validation_results \\\n",
      "    --failed_ldsc_action drop_traits \\\n",
      "    --validate_only\n\n",
      "Run with GWAMA by additionally providing --gpca_input_folder and\n",
      "--source_path, and omit --validate_only."
    )
  )

  required_inputs <- parser$add_argument_group("Required inputs")
  gwama_inputs <- parser$add_argument_group(
    "GWAMA inputs (not required with --validate_only)"
  )
  trait_handling <- parser$add_argument_group("Trait and failed-result handling")
  qc_thresholds <- parser$add_argument_group("QC comparison thresholds")
  pca_settings <- parser$add_argument_group("PCA settings")
  diagnostic_actions <- parser$add_argument_group("Diagnostic warning/error controls")
  execution <- parser$add_argument_group("Execution and performance")

  parser$add_argument(
    "-h", "--help",
    action = "store_true",
    help = "Show this help message and exit."
  )

  required_inputs$add_argument(
    "--input",
    required = TRUE,
    metavar = "MANIFEST.csv",
    help = paste(
      "CSV trait manifest. Required column: 'traitname'. No other manifest",
      "columns are required; extra columns are ignored. Trait names must be",
      "unique/non-empty, and row order defines matrix and GWAMA order."
    )
  )
  required_inputs$add_argument(
    "--python_ldsc",
    required = TRUE,
    metavar = "LDSC.csv[.gz]",
    help = paste(
      "Complete pairwise Python LDSC table in CSV/TSV format; gzip is",
      "accepted. Required columns are listed below."
    )
  )
  required_inputs$add_argument(
    "--outdir",
    required = TRUE,
    metavar = "DIRECTORY",
    help = "Directory for QC audit files, matrices, weights, and GWAMA results."
  )

  gwama_inputs$add_argument(
    "--gpca_input_folder",
    default = NULL,
    metavar = "DIRECTORY",
    help = "Directory containing the nine-column per-trait GWAMA input files."
  )
  gwama_inputs$add_argument(
    "--source_path",
    default = NULL,
    metavar = "GWAMA_FUNCTION.R",
    help = "R file defining the already-modified my_GWAMA() function."
  )
  gwama_inputs$add_argument(
    "--splitby_chr",
    default = "split",
    choices = c("split", "nosplit"),
    help = paste(
      "GWAMA input layout: 'split' uses chromosomes 1-22; 'nosplit' uses",
      "one whole-genome file per trait. Default: split."
    )
  )

  trait_handling$add_argument(
    "--allow_missing_traits",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Remove manifest traits absent from the LDSC table. Without this flag,",
      "an absent selected trait stops the run. Extra LDSC traits are ignored."
    )
  )
  trait_handling$add_argument(
    "--failed_ldsc_action",
    choices = c("error", "drop_traits"),
    default = "error",
    help = paste(
      "For missing/non-finite/invalid selected pairs: 'error' stops;",
      "'drop_traits' removes traits until a complete valid subset remains.",
      "No estimates are imputed. Default: error."
    )
  )

  qc_thresholds$add_argument(
    "--duplicate_tolerance", "--tolerance",
    dest = "duplicate_tolerance",
    type = "double",
    metavar = "FLOAT",
    default = 1e-3,
    help = paste(
      "Maximum absolute difference between duplicate orientations for rg,",
      "SE, p, and intercept estimates. Default: 0.001."
    )
  )
  qc_thresholds$add_argument(
    "--duplicate_z_tolerance",
    type = "double",
    metavar = "FLOAT",
    default = 1e-2,
    help = "Maximum absolute duplicate-orientation z difference. Default: 0.01."
  )
  qc_thresholds$add_argument(
    "--self_rg_tolerance",
    type = "double",
    metavar = "FLOAT",
    default = 1e-2,
    help = paste(
      "Maximum |self rg - 1|. With the default 0.01, values from 0.99",
      "through 1.01 pass."
    )
  )
  qc_thresholds$add_argument(
    "--comparison_epsilon",
    type = "double",
    metavar = "FLOAT",
    default = 1e-12,
    help = paste(
      "Floating-point allowance added to boundary comparisons. Default:",
      "1e-12."
    )
  )

  pca_settings$add_argument(
    "--heritability_scale",
    choices = c("auto", "liability", "observed", "mixed"),
    default = "auto",
    help = paste(
      "Heritability columns used for self-pair QC and covariance PCA. 'auto'",
      "uses populated values for selected traits, ignoring empty pairs.",
      "Correlation: per-trait scale; covariance: one complete common scale.",
      "'mixed' explicitly permits trait-specific covariance scales, only when",
      "each trait has one populated scale. Ambiguity stops; no conversion.",
      "Default: auto."
    )
  )
  pca_settings$add_argument(
    "--pca_matrix",
    choices = c("correlation", "covariance"),
    default = "correlation",
    help = paste(
      "Matrix decomposed for PC1: 'correlation' uses S_Stand; 'covariance'",
      "derives genetic covariance from rg and selected-scale self h2, following the",
      "tutorial alternative procedure. CTI is unchanged. Default: correlation."
    )
  )
  pca_settings$add_argument(
    "--pc1_orientation",
    choices = c("tutorial", "as_computed"),
    default = "tutorial",
    help = paste(
      "PC1 sign convention. 'tutorial' flips the entire loading vector when",
      "its median is negative; 'as_computed' keeps R's eigenvector sign.",
      "Association strength and p-values are unchanged. Default: tutorial."
    )
  )

  diagnostic_actions$add_argument(
    "--h2_z_warn_threshold",
    type = "double",
    metavar = "FLOAT",
    default = 2,
    help = paste(
      "Warn when retained self-pair h2/SE is below this value. Diagnostic",
      "only; never removes traits. Use 0 to disable. Default: 2."
    )
  )
  diagnostic_actions$add_argument(
    "--z_consistency_tolerance",
    type = "double",
    metavar = "FLOAT",
    default = 1e-2,
    help = paste(
      "Maximum relative difference between reported z and rg/SE after",
      "duplicate collapse. Default: 0.01."
    )
  )
  diagnostic_actions$add_argument(
    "--z_consistency_action",
    choices = c("warn", "error"),
    default = "warn",
    help = "Warn or stop when z is inconsistent with rg/SE. Default: warn."
  )
  diagnostic_actions$add_argument(
    "--rg_out_of_range_action",
    choices = c("warn", "error"),
    default = "warn",
    help = paste(
      "Warn or stop for finite off-diagonal rg outside [-1,1]. Values are",
      "never clamped. Default: warn."
    )
  )
  diagnostic_actions$add_argument(
    "--negative_eigen_action",
    choices = c("warn", "error"),
    default = "warn",
    help = paste(
      "Warn or stop for substantive negative selected-PCA-matrix",
      "eigenvalues. Raw values are reported. Default: warn."
    )
  )
  diagnostic_actions$add_argument(
    "--matrix_eigen_tolerance",
    type = "double",
    metavar = "FLOAT",
    default = 1e-8,
    help = paste(
      "Relative cutoff separating substantive negative eigenvalues from",
      "floating-point noise. Default: 1e-8."
    )
  )

  execution$add_argument(
    "--validate_only",
    action = "store_true",
    default = FALSE,
    help = "Run LDSC QC and PCA and write audit files, but do not run GWAMA."
  )
  execution$add_argument(
    "--cores",
    type = "integer",
    metavar = "N",
    default = 0L,
    help = paste(
      "Workers for chromosome-split GWAMA. 0 selects up to 22 based on",
      "available physical cores. Default: 0."
    )
  )
  execution$add_argument(
    "--ldsc_chunk_size",
    type = "integer",
    metavar = "N",
    default = 250000L,
    help = "Number of LDSC rows read per streaming chunk. Default: 250000."
  )

  execution$add_argument(
    "--color",
    choices = c("auto", "always", "never"),
    default = "auto",
    help = paste(
      "Help-message colour: auto for terminals, always to force it, or never",
      "to disable it. NO_COLOR is respected in auto mode. Default: auto."
    )
  )

  if (!length(raw_arguments) || any(raw_arguments %in% c("-h", "--help"))) {
    color_mode <- requested_cli_color_mode(raw_arguments)
    writeLines(
      colorize_cli_help(
        parser$format_help(),
        enabled = cli_color_enabled(color_mode)
      )
    )
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  tryCatch(parser$parse_args(raw_arguments), error = function(e) {
    parser$print_help()
    stop(e)
  })
}

validate_cli_paths <- function(args) {
  args$splitby_chr <- tolower(trimws(args$splitby_chr))

  if (!file.exists(args$input)) {
    stop(glue("Trait manifest not found: {args$input}"), call. = FALSE)
  }
  if (!file.exists(args$python_ldsc)) {
    stop(glue("Python LDSC results not found: {args$python_ldsc}"), call. = FALSE)
  }
  if (!isTRUE(args$validate_only)) {
    if (is.null(args$gpca_input_folder) || !nzchar(args$gpca_input_folder)) {
      stop("--gpca_input_folder is required unless --validate_only is used.", call. = FALSE)
    }
    if (!dir.exists(args$gpca_input_folder)) {
      stop(
        glue("Genomic PCA input folder not found: {args$gpca_input_folder}"),
        call. = FALSE
      )
    }
    if (is.null(args$source_path) || !nzchar(args$source_path)) {
      stop("--source_path is required unless --validate_only is used.", call. = FALSE)
    }
    if (!file.exists(args$source_path)) {
      stop(glue("GWAMA function script not found: {args$source_path}"), call. = FALSE)
    }
  }
  if (!is.finite(args$duplicate_tolerance) || args$duplicate_tolerance <= 0) {
    stop(
      "--duplicate_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (!is.finite(args$duplicate_z_tolerance) ||
      args$duplicate_z_tolerance <= 0) {
    stop(
      "--duplicate_z_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (!is.finite(args$self_rg_tolerance) || args$self_rg_tolerance <= 0) {
    stop(
      "--self_rg_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (!is.finite(args$comparison_epsilon) || args$comparison_epsilon < 0) {
    stop(
      "--comparison_epsilon must be finite and non-negative.",
      call. = FALSE
    )
  }
  if (!is.finite(args$h2_z_warn_threshold) ||
      args$h2_z_warn_threshold < 0) {
    stop(
      "--h2_z_warn_threshold must be finite and non-negative.",
      call. = FALSE
    )
  }
  if (!is.finite(args$z_consistency_tolerance) ||
      args$z_consistency_tolerance <= 0) {
    stop(
      "--z_consistency_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (!is.finite(args$matrix_eigen_tolerance) ||
      args$matrix_eigen_tolerance <= 0) {
    stop(
      "--matrix_eigen_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (is.na(args$ldsc_chunk_size) || args$ldsc_chunk_size < 1L) {
    stop("--ldsc_chunk_size must be a positive integer.", call. = FALSE)
  }
  if (is.na(args$cores) || args$cores < 0L) {
    stop("--cores must be zero or a positive integer.", call. = FALSE)
  }

  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(args$outdir)) {
    stop(glue("Unable to create output directory: {args$outdir}"), call. = FALSE)
  }

  args
}
