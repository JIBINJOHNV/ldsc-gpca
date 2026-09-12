#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# GENOMIC PCA & GWAMA FROM PAIRWISE PYTHON LDSC OUTPUT
# -----------------------------------------------------------------------------
# Implements both PCA procedures in Anna Fuertjes' tutorial:
# https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html
#
# The Python LDSC table is used to construct only the objects needed by the
# selected tutorial procedure:
#   * correlation PCA (default): diagonal 1, off-diagonal pairwise rg
#   * covariance PCA: diagonal selected-scale self h2, off-diagonal
#                     rg * sqrt(h2_1 * h2_2)
#   * CTI:       diagonal univariate LDSC intercept (h2_int), off-diagonal
#                bivariate LDSC intercept (gcov_int)
#
# Python ldsc.py's marginal standard errors do not reproduce GenomicSEM's full
# V or V_Stand matrices. This script therefore does not fabricate those objects
# and does not run paLDSC.
#
# The manifest is general-purpose: it may select any number of unique traits
# supported by available memory. Only its traitname column is interpreted;
# other columns (for example labels or paths) do not define analysis groups.
# Extra traits in the LDSC file are ignored. Strict failure is the default;
# --allow_missing_traits and --failed_ldsc_action=drop_traits are explicit,
# audited opt-ins that reduce the selected set without imputing LDSC values.
#
# --source_path must point to the already-modified Fuertjes GWAMA function. Its
# weight calculation must be equivalent to:
#
#   W <- t(t(sqrt(N)) * h2)
#   sqrt_W <- W
#
# In that function call, h2 contains the selected PC1 loading vector, not SNP
# heritability; the argument name is inherited from the supplied GWAMA code.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(glue)
  library(parallel)
})

required_python_ldsc_base_columns <- c(
  "p1", "p2", "rg", "se", "z", "p",
  "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"
)

heritability_column_sets <- list(
  liability = c("h2_liab", "h2_liab_se"),
  observed = c("h2_obs", "h2_obs_se")
)

internal_heritability_columns <- c("h2", "h2_se")

# All processing after input normalization uses scale-neutral h2/h2_se names.
required_python_ldsc_columns <- c(
  required_python_ldsc_base_columns,
  internal_heritability_columns
)

python_ldsc_numeric_columns <- setdiff(
  required_python_ldsc_columns,
  c("p1", "p2")
)

duplicate_comparison_columns <- c(
  "rg", "se", "z", "p", "gcov_int", "gcov_int_se"
)

gwama_required_columns <- c(
  "SNPID", "CHR", "BP", "EA", "OA", "EAF", "N", "Z", "P"
)

# Colour is applied after argparse has formatted and wrapped the text. This
# keeps ANSI escape sequences from disturbing help alignment. Automatic mode
# follows terminal capability and the widely used NO_COLOR/CLICOLOR_FORCE
# environment conventions.
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
      "REQUIRED MANIFEST COLUMN\n",
      "  traitname    Unique, non-empty trait identifier. Each row selects one\n",
      "               trait, and row order defines all analysis matrices.\n",
      "  No other manifest columns are required. Extra columns are ignored.\n\n",
      "REQUIRED PYTHON LDSC COLUMNS\n",
      "  Always required:\n",
      "    p1, p2, rg, se, z, p, h2_int, h2_int_se, gcov_int, gcov_int_se\n",
      "  Plus one complete heritability pair:\n",
      "    h2_liab, h2_liab_se  OR  h2_obs, h2_obs_se\n",
      "  If both pairs exist, select one with --heritability_scale.\n\n",
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
    choices = c("auto", "liability", "observed"),
    default = "auto",
    help = paste(
      "Heritability columns used for self-pair QC and covariance PCA. 'auto'",
      "selects the only complete pair; if both pairs exist, choose",
      "'liability' or 'observed' explicitly. Scales are never mixed.",
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

  if (any(raw_arguments %in% c("-h", "--help"))) {
    color_mode <- requested_cli_color_mode(raw_arguments)
    writeLines(
      colorize_cli_help(
        parser$format_help(),
        enabled = cli_color_enabled(color_mode)
      )
    )
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  parser$parse_args(raw_arguments)
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

read_trait_manifest <- function(path) {
  manifest <- fread(path, data.table = FALSE, check.names = FALSE)

  if (!"traitname" %in% names(manifest)) {
    stop("The input CSV must contain a column named 'traitname'.", call. = FALSE)
  }

  traits <- trimws(as.character(manifest$traitname))

  if (anyNA(traits) || any(traits == "")) {
    stop("The traitname column contains missing or empty names.", call. = FALSE)
  }
  if (anyDuplicated(traits)) {
    duplicates <- unique(traits[duplicated(traits)])
    stop(
      "The traitname column must be unique. Duplicates:\n",
      paste(duplicates, collapse = "\n"),
      call. = FALSE
    )
  }
  if (length(traits) < 2L) {
    stop("Genomic PCA/GWAMA requires at least two traits.", call. = FALSE)
  }

  traits
}

open_ldsc_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

detect_delimiter <- function(header_line) {
  comma_count <- lengths(regmatches(header_line, gregexpr(",", header_line)))
  tab_count <- lengths(regmatches(header_line, gregexpr("\\t", header_line)))

  if (tab_count > comma_count && tab_count > 0L) {
    "\t"
  } else if (comma_count > 0L) {
    ","
  } else {
    "auto"
  }
}

# Python LDSC writes either an observed-scale or liability-scale h2 pair,
# depending on whether prevalence conversion was requested. Resolve that pair
# once for the whole file; never combine scales across rows or traits.
resolve_heritability_columns <- function(column_names,
                                         requested_scale = c(
                                           "auto", "liability", "observed"
                                         )) {
  requested_scale <- match.arg(requested_scale)
  column_names <- as.character(column_names)

  completeness <- vapply(
    heritability_column_sets,
    function(columns) all(columns %in% column_names),
    logical(1)
  )
  partial <- vapply(
    heritability_column_sets,
    function(columns) {
      any(columns %in% column_names) && !all(columns %in% column_names)
    },
    logical(1)
  )

  if (any(partial)) {
    details <- vapply(
      names(partial)[partial],
      function(scale) {
        missing <- setdiff(heritability_column_sets[[scale]], column_names)
        paste0(scale, " scale is missing: ", paste(missing, collapse = ", "))
      },
      character(1)
    )
    stop(
      "Python LDSC input contains an incomplete heritability column pair: ",
      paste(details, collapse = "; "),
      ". Supply both columns for a scale.",
      call. = FALSE
    )
  }

  complete_scales <- names(completeness)[completeness]
  if (requested_scale == "auto") {
    if (length(complete_scales) == 0L) {
      stop(
        "Python LDSC input must contain one complete heritability pair: ",
        "h2_liab with h2_liab_se, or h2_obs with h2_obs_se.",
        call. = FALSE
      )
    }
    if (length(complete_scales) > 1L) {
      stop(
        "Python LDSC input contains both liability- and observed-scale ",
        "heritability pairs. Select one explicitly with ",
        "--heritability_scale liability or --heritability_scale observed.",
        call. = FALSE
      )
    }
    selected_scale <- complete_scales
  } else {
    selected_scale <- requested_scale
    if (!completeness[[selected_scale]]) {
      required <- heritability_column_sets[[selected_scale]]
      stop(
        glue(
          "--heritability_scale {selected_scale} requires columns: ",
          "{paste(required, collapse = ', ')}."
        ),
        call. = FALSE
      )
    }
  }

  selected_columns <- heritability_column_sets[[selected_scale]]
  list(
    scale = selected_scale,
    h2_column = selected_columns[1L],
    h2_se_column = selected_columns[2L]
  )
}

normalize_heritability_columns <- function(ldsc_rows,
                                           requested_scale = c(
                                             "auto", "liability", "observed"
                                           )) {
  requested_scale <- match.arg(requested_scale)
  inherited_scale <- attr(ldsc_rows, "heritability_scale", exact = TRUE)
  internal_present <- internal_heritability_columns %in% names(ldsc_rows)

  if (any(internal_present) && !all(internal_present)) {
    stop(
      "Internal LDSC heritability columns are incomplete; both h2 and h2_se are required.",
      call. = FALSE
    )
  }

  if (all(internal_present)) {
    if (is.null(inherited_scale) || !inherited_scale %in% c("liability", "observed")) {
      if (requested_scale == "auto") {
        stop(
          "Normalized h2/h2_se columns lack their heritability-scale metadata.",
          call. = FALSE
        )
      }
      inherited_scale <- requested_scale
    }
    if (requested_scale != "auto" && requested_scale != inherited_scale) {
      stop(
        glue(
          "Requested heritability scale '{requested_scale}' conflicts with ",
          "the normalized input scale '{inherited_scale}'."
        ),
        call. = FALSE
      )
    }
    normalized <- as.data.table(copy(ldsc_rows))
    selected_columns <- heritability_column_sets[[inherited_scale]]
    attr(normalized, "heritability_scale") <- inherited_scale
    attr(normalized, "heritability_value_column") <- selected_columns[1L]
    attr(normalized, "heritability_se_column") <- selected_columns[2L]
    return(normalized)
  }

  resolution <- resolve_heritability_columns(
    names(ldsc_rows),
    requested_scale
  )
  normalized <- as.data.table(copy(ldsc_rows))
  setnames(
    normalized,
    c(resolution$h2_column, resolution$h2_se_column),
    internal_heritability_columns
  )

  unused_heritability_columns <- setdiff(
    intersect(unlist(heritability_column_sets), names(normalized)),
    internal_heritability_columns
  )
  if (length(unused_heritability_columns) > 0L) {
    normalized[, (unused_heritability_columns) := NULL]
  }

  attr(normalized, "heritability_scale") <- resolution$scale
  attr(normalized, "heritability_value_column") <- resolution$h2_column
  attr(normalized, "heritability_se_column") <- resolution$h2_se_column
  normalized
}

# Stream the potentially multi-million-row table and retain only selected
# trait-by-trait rows. This avoids holding a larger all-trait table in memory.
read_python_ldsc_selected <- function(path, trait_order, chunk_size = 250000L,
                                      heritability_scale = c(
                                        "auto", "liability", "observed"
                                      )) {
  heritability_scale <- match.arg(heritability_scale)
  connection <- open_ldsc_connection(path)
  on.exit(close(connection), add = TRUE)

  header <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header) != 1L) {
    stop("The Python LDSC file is empty.", call. = FALSE)
  }

  separator <- detect_delimiter(header)
  header_table <- fread(
    text = header,
    nrows = 0L,
    sep = separator,
    check.names = FALSE,
    showProgress = FALSE
  )
  missing_columns <- setdiff(
    required_python_ldsc_base_columns,
    names(header_table)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Python LDSC input is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  heritability_resolution <- resolve_heritability_columns(
    names(header_table),
    heritability_scale
  )
  source_columns <- c(
    required_python_ldsc_base_columns,
    heritability_resolution$h2_column,
    heritability_resolution$h2_se_column
  )

  retained_chunks <- list()
  retained_index <- 0L
  rows_read <- 0L
  selected_traits_seen <- rep(FALSE, length(trait_order))

  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    rows_read <- rows_read + length(lines)
    chunk <- fread(
      text = paste(c(header, lines), collapse = "\n"),
      sep = separator,
      select = source_columns,
      check.names = FALSE,
      showProgress = FALSE,
      na.strings = c("NA", "NaN", "nan", "")
    )
    setnames(
      chunk,
      c(
        heritability_resolution$h2_column,
        heritability_resolution$h2_se_column
      ),
      internal_heritability_columns
    )

    chunk[, p1 := trimws(as.character(p1))]
    chunk[, p2 := trimws(as.character(p2))]
    newly_seen <- unique(c(
      chunk$p1[chunk$p1 %chin% trait_order],
      chunk$p2[chunk$p2 %chin% trait_order]
    ))
    selected_traits_seen[match(newly_seen, trait_order)] <- TRUE
    chunk <- chunk[p1 %chin% trait_order & p2 %chin% trait_order]

    if (nrow(chunk) > 0L) {
      retained_index <- retained_index + 1L
      retained_chunks[[retained_index]] <- chunk
    }

    if (rows_read %% (4L * chunk_size) == 0L) {
      message(glue("  parsed {format(rows_read, big.mark = ',')} LDSC rows"))
    }
  }

  if (length(retained_chunks) == 0L) {
    selected <- as.data.table(setNames(
      replicate(length(required_python_ldsc_columns), logical(0), simplify = FALSE),
      required_python_ldsc_columns
    ))
  } else {
    selected <- rbindlist(retained_chunks, use.names = TRUE)
  }

  attr(selected, "source_rows_read") <- rows_read
  attr(selected, "selected_traits_seen") <- trait_order[selected_traits_seen]
  attr(selected, "heritability_scale") <- heritability_resolution$scale
  attr(selected, "heritability_value_column") <-
    heritability_resolution$h2_column
  attr(selected, "heritability_se_column") <-
    heritability_resolution$h2_se_column
  selected
}

resolve_ldsc_trait_order <- function(ldsc_rows, manifest_traits,
                                     allow_missing_traits = FALSE) {
  traits_seen <- attr(ldsc_rows, "selected_traits_seen", exact = TRUE)
  if (is.null(traits_seen)) {
    traits_seen <- unique(c(
      as.character(ldsc_rows$p1),
      as.character(ldsc_rows$p2)
    ))
  }

  present <- manifest_traits %chin% traits_seen
  missing_table <- data.frame(
    Manifest_Order = which(!present),
    Trait = manifest_traits[!present],
    Reason = rep("Absent from Python LDSC p1/p2 columns", sum(!present)),
    stringsAsFactors = FALSE
  )

  if (nrow(missing_table) > 0L && !isTRUE(allow_missing_traits)) {
    stop(
      "Selected traits absent from the Python LDSC file:\n",
      paste(missing_table$Trait, collapse = "\n"),
      paste0(
        "\nRe-run with --allow_missing_traits to drop these traits and ",
        "continue with those present."
      ),
      call. = FALSE
    )
  }

  retained_traits <- manifest_traits[present]
  if (length(retained_traits) < 2L) {
    stop(
      glue(
        "Only {length(retained_traits)} manifest trait(s) are present in ",
        "Python LDSC; at least two are required for genomic PCA/GWAMA."
      ),
      call. = FALSE
    )
  }

  if (nrow(missing_table) > 0L) {
    warning(
      glue(
        "Dropping {nrow(missing_table)} manifest trait(s) absent from Python ",
        "LDSC because --allow_missing_traits was supplied: ",
        "{paste(missing_table$Trait, collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  list(
    trait_order = retained_traits,
    missing_traits = missing_table
  )
}

coerce_python_ldsc_numeric <- function(ldsc_rows) {
  for (column_name in python_ldsc_numeric_columns) {
    original <- ldsc_rows[[column_name]]
    converted <- suppressWarnings(as.numeric(original))
    invalid_conversion <- is.na(converted) & !is.na(original)

    if (any(invalid_conversion)) {
      example <- as.character(original[which(invalid_conversion)[1L]])
      stop(
        glue(
          "Python LDSC column '{column_name}' contains a non-numeric value: ",
          "'{example}'."
        ),
        call. = FALSE
      )
    }

    set(ldsc_rows, j = column_name, value = converted)
  }

  ldsc_rows
}

python_ldsc_valid_row <- function(ldsc_rows) {
  finite_numeric <- Reduce(
    `&`,
    lapply(ldsc_rows[, ..python_ldsc_numeric_columns], is.finite)
  )
  positive_se <- Reduce(
    `&`,
    lapply(
      ldsc_rows[, .(se, h2_se, h2_int_se, gcov_int_se)],
      function(value) is.finite(value) & value > 0
    )
  )
  valid_p <- is.finite(ldsc_rows$p) & ldsc_rows$p >= 0 & ldsc_rows$p <= 1

  finite_numeric & positive_se & valid_p
}

# Evaluate every selected trait's LDSC self-pair in one place. This function is
# deliberately independent of the pair-removal algorithm so that strict/error
# mode, drop-traits mode, final validation, and the audit file use identical
# rules. A low h2/SE is diagnostic only; it is never part of Self_QC_Pass.
evaluate_self_pair_qc <- function(ldsc_rows, trait_order,
                                  self_rg_tolerance = 1e-2,
                                  comparison_epsilon = 1e-12,
                                  h2_z_warn_threshold = 2,
                                  heritability_scale = c(
                                    "auto", "liability", "observed"
                                  )) {
  heritability_scale <- match.arg(heritability_scale)
  if (!is.finite(self_rg_tolerance) || self_rg_tolerance <= 0) {
    stop("self_rg_tolerance must be finite and greater than zero.", call. = FALSE)
  }
  if (!is.finite(comparison_epsilon) || comparison_epsilon < 0) {
    stop("comparison_epsilon must be finite and non-negative.", call. = FALSE)
  }
  if (!is.finite(h2_z_warn_threshold) || h2_z_warn_threshold < 0) {
    stop("h2_z_warn_threshold must be finite and non-negative.", call. = FALSE)
  }

  ldsc_rows <- normalize_heritability_columns(
    ldsc_rows,
    heritability_scale
  )
  selected_scale <- attr(ldsc_rows, "heritability_scale", exact = TRUE)
  h2_source_column <- attr(
    ldsc_rows,
    "heritability_value_column",
    exact = TRUE
  )
  h2_se_source_column <- attr(
    ldsc_rows,
    "heritability_se_column",
    exact = TRUE
  )

  missing_columns <- setdiff(required_python_ldsc_columns, names(ldsc_rows))
  if (length(missing_columns) > 0L) {
    stop(
      "Python LDSC input is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ldsc_rows[, p1 := trimws(as.character(p1))]
  ldsc_rows[, p2 := trimws(as.character(p2))]
  ldsc_rows <- ldsc_rows[p1 %chin% trait_order & p2 %chin% trait_order]
  ldsc_rows <- coerce_python_ldsc_numeric(ldsc_rows)

  trait_index_1 <- match(ldsc_rows$p1, trait_order)
  trait_index_2 <- match(ldsc_rows$p2, trait_order)
  ldsc_rows[, pair_i := pmin(trait_index_1, trait_index_2)]
  ldsc_rows[, pair_j := pmax(trait_index_1, trait_index_2)]
  self_rows <- ldsc_rows[pair_i == pair_j]

  finite_mean <- function(value) {
    if (length(value) > 0L && all(is.finite(value))) mean(value) else NA_real_
  }
  summarize_self <- function(rows) {
    data.table(
      Self_Source_Rows = nrow(rows),
      Required_Numeric_Finite = nrow(rows) > 0L && all(
        vapply(
          rows[, ..python_ldsc_numeric_columns],
          function(value) all(is.finite(value)),
          logical(1)
        )
      ),
      Required_SE_Positive = nrow(rows) > 0L && all(
        vapply(
          rows[, .(se, h2_se, h2_int_se, gcov_int_se)],
          function(value) all(is.finite(value) & value > 0),
          logical(1)
        )
      ),
      P_Valid = nrow(rows) > 0L && all(
        is.finite(rows$p) & rows$p >= 0 & rows$p <= 1
      ),
      Self_RG = finite_mean(rows$rg),
      Self_RG_SE = finite_mean(rows$se),
      Self_Z = finite_mean(rows$z),
      Self_P = finite_mean(rows$p),
      Self_H2 = finite_mean(rows$h2),
      Self_H2_SE = finite_mean(rows$h2_se),
      Self_H2_Intercept = finite_mean(rows$h2_int),
      Self_H2_Intercept_SE = finite_mean(rows$h2_int_se),
      Self_Gcov_Intercept = finite_mean(rows$gcov_int),
      Self_Gcov_Intercept_SE = finite_mean(rows$gcov_int_se)
    )
  }

  summaries <- lapply(seq_along(trait_order), function(trait_index) {
    summarize_self(self_rows[pair_i == trait_index])
  })
  result <- rbindlist(summaries)
  result[, `:=`(
    Manifest_Order = seq_along(trait_order),
    Trait = trait_order,
    Self_Pair_Found = Self_Source_Rows > 0L
  )]
  setcolorder(
    result,
    c(
      "Manifest_Order", "Trait", "Self_Pair_Found", "Self_Source_Rows",
      setdiff(names(result), c("Manifest_Order", "Trait", "Self_Pair_Found", "Self_Source_Rows"))
    )
  )

  result[, Self_RG_Deviation := abs(Self_RG - 1)]
  result[, `:=`(
    Self_RG_Lower_Bound = 1 - self_rg_tolerance,
    Self_RG_Upper_Bound = 1 + self_rg_tolerance,
    Self_RG_Within_Tolerance = is.finite(Self_RG) &
      Self_RG_Deviation <= self_rg_tolerance + comparison_epsilon,
    Self_H2_Positive = is.finite(Self_H2) & Self_H2 > 0,
    Self_H2_SE_Positive = is.finite(Self_H2_SE) & Self_H2_SE > 0,
    H2_Z = Self_H2 / Self_H2_SE
  )]
  result[, H2_Z_Below_Warning_Threshold := if (h2_z_warn_threshold > 0) {
    is.finite(H2_Z) & H2_Z < h2_z_warn_threshold
  } else {
    rep(FALSE, .N)
  }]

  result[, Self_QC_Pass :=
    Self_Pair_Found &
      Required_Numeric_Finite &
      Required_SE_Positive &
      P_Valid &
      Self_H2_Positive &
      Self_RG_Within_Tolerance]
  result[, Self_QC_Failure_Reason := fifelse(
    !Self_Pair_Found,
    "Missing LDSC self-pair",
    fifelse(
      !Required_Numeric_Finite,
      "Non-finite required LDSC self-pair value",
      fifelse(
        !Required_SE_Positive,
        "Non-positive required LDSC self-pair standard error",
        fifelse(
          !P_Valid,
          "Self-pair p-value outside [0,1]",
          fifelse(
            !Self_H2_Positive,
            paste0("Non-positive self-pair ", h2_source_column),
            fifelse(
              !Self_RG_Within_Tolerance,
              "Self-pair rg differs from 1 beyond tolerance",
              ""
            )
          )
        )
      )
    )
  )]
  result[, `:=`(
    Heritability_Scale = selected_scale,
    Heritability_Source_Column = h2_source_column,
    Heritability_SE_Source_Column = h2_se_source_column,
    Self_RG_Tolerance = self_rg_tolerance,
    Comparison_Epsilon = comparison_epsilon,
    H2_Z_Warn_Threshold = h2_z_warn_threshold
  )]

  as.data.frame(result)
}

empty_failed_trait_table <- function() {
  data.frame(
    Exclusion_Step = integer(),
    Manifest_Order = integer(),
    Trait = character(),
    Reason = character(),
    Failed_Pairs_At_Exclusion = integer(),
    h2_liab = double(),
    h2_liab_se = double(),
    h2_Z = double(),
    Heritability_Scale = character(),
    Heritability = double(),
    Heritability_SE = double(),
    h2_obs = double(),
    h2_obs_se = double(),
    stringsAsFactors = FALSE
  )
}

# Missing/non-finite pair estimates cannot be inserted into a genomic PCA
# matrix. In opt-in drop mode, find a deterministic complete finite subset:
#   1. remove traits whose self-pair is missing/invalid or whose self h2 <= 0;
#   2. repeatedly remove the trait incident to the most remaining failed pairs;
#   3. break ties by lower self-h2 Z, then later manifest position.
# This is a scalable greedy vertex-cover heuristic. It does not claim to find
# the mathematically largest possible subset (that problem is NP-hard), and it
# never imputes an LDSC estimate.
resolve_incomplete_ldsc_traits <- function(ldsc_rows, trait_order,
                                           action = c("error", "drop_traits"),
                                           self_rg_tolerance = 1e-2,
                                           comparison_epsilon = 1e-12,
                                           h2_z_warn_threshold = 2,
                                           heritability_scale = c(
                                             "auto", "liability", "observed"
                                           )) {
  action <- match.arg(action)
  heritability_scale <- match.arg(heritability_scale)
  source_rows_read <- attr(ldsc_rows, "source_rows_read", exact = TRUE)
  ldsc_rows <- normalize_heritability_columns(
    ldsc_rows,
    heritability_scale
  )
  selected_scale <- attr(ldsc_rows, "heritability_scale", exact = TRUE)
  h2_source_column <- attr(
    ldsc_rows,
    "heritability_value_column",
    exact = TRUE
  )
  h2_se_source_column <- attr(
    ldsc_rows,
    "heritability_se_column",
    exact = TRUE
  )
  ldsc_rows <- ldsc_rows[p1 %chin% trait_order & p2 %chin% trait_order]
  ldsc_rows <- coerce_python_ldsc_numeric(ldsc_rows)

  number_traits <- length(trait_order)
  index_1 <- match(ldsc_rows$p1, trait_order)
  index_2 <- match(ldsc_rows$p2, trait_order)
  ldsc_rows[, pair_i := pmin(index_1, index_2)]
  ldsc_rows[, pair_j := pmax(index_1, index_2)]
  ldsc_rows[, valid_ldsc_row := python_ldsc_valid_row(.SD)]

  pair_status <- ldsc_rows[, .(
    Source_Rows = .N,
    Invalid_Source_Rows = sum(!valid_ldsc_row),
    Pair_Valid = all(valid_ldsc_row)
  ), by = .(pair_i, pair_j)]

  expected_pairs <- CJ(
    pair_i = seq_len(number_traits),
    pair_j = seq_len(number_traits)
  )[pair_i <= pair_j]
  missing_pairs <- expected_pairs[!pair_status, on = .(pair_i, pair_j)]
  missing_pairs[, `:=`(
    Source_Rows = 0L,
    Invalid_Source_Rows = 0L,
    Pair_Valid = FALSE
  )]
  failed_pairs <- rbindlist(
    list(pair_status[Pair_Valid == FALSE], missing_pairs),
    use.names = TRUE
  )
  setorder(failed_pairs, pair_i, pair_j)

  self_pair_qc <- evaluate_self_pair_qc(
    ldsc_rows,
    trait_order,
    self_rg_tolerance = self_rg_tolerance,
    comparison_epsilon = comparison_epsilon,
    h2_z_warn_threshold = h2_z_warn_threshold,
    heritability_scale = selected_scale
  )
  self_pair_qc <- as.data.table(self_pair_qc)
  self_lookup <- data.table(
    pair_i = seq_len(number_traits),
    rg = self_pair_qc$Self_RG,
    h2 = self_pair_qc$Self_H2,
    h2_se = self_pair_qc$Self_H2_SE,
    h2_Z = self_pair_qc$H2_Z
  )

  invalid_self_indices <- self_pair_qc[Self_QC_Pass == FALSE, Manifest_Order]

  if (nrow(failed_pairs) == 0L && length(invalid_self_indices) == 0L) {
    if (!is.null(source_rows_read)) {
      attr(ldsc_rows, "source_rows_read") <- source_rows_read
    }
    attr(ldsc_rows, "heritability_scale") <- selected_scale
    attr(ldsc_rows, "heritability_value_column") <- h2_source_column
    attr(ldsc_rows, "heritability_se_column") <- h2_se_source_column
    return(list(
      ldsc_rows = ldsc_rows,
      trait_order = trait_order,
      excluded_traits = empty_failed_trait_table(),
      initial_failed_pair_count = 0L,
      initial_self_qc_failure_count = 0L,
      self_pair_qc = as.data.frame(self_pair_qc),
      heritability_scale = selected_scale
    ))
  }

  pair_descriptions <- paste0(
    trait_order[failed_pairs$pair_i], " <-> ",
    trait_order[failed_pairs$pair_j],
    ifelse(
      failed_pairs$Source_Rows == 0L,
      " [missing]",
      " [non-finite or invalid numeric result]"
    )
  )
  self_qc_only <- setdiff(
    invalid_self_indices,
    failed_pairs[pair_i == pair_j, pair_i]
  )
  self_descriptions <- paste0(
    trait_order[self_qc_only],
    " [",
    self_pair_qc$Self_QC_Failure_Reason[self_qc_only],
    "]"
  )
  all_descriptions <- c(pair_descriptions, self_descriptions)
  if (action == "error") {
    stop(
      glue(
        "The selected LDSC set contains {nrow(failed_pairs)} missing/invalid ",
        "unique pair estimate(s) and {length(self_qc_only)} additional ",
        "self-pair QC failure(s):\n"
      ),
      paste(head(all_descriptions, 50L), collapse = "\n"),
      if (length(all_descriptions) > 50L) {
        glue("\n... and {length(all_descriptions) - 50L} more")
      } else {
        ""
      },
      paste0(
        "\nRe-run with --failed_ldsc_action drop_traits to construct a ",
        "complete finite subset without imputing pair estimates."
      ),
      call. = FALSE
    )
  }

  initial_failed_pair_count <- nrow(failed_pairs)
  active <- rep(TRUE, number_traits)
  exclusions <- list()
  exclusion_step <- 0L

  add_exclusion <- function(trait_index, reason, failed_degree) {
    exclusion_step <<- exclusion_step + 1L
    exclusions[[exclusion_step]] <<- data.frame(
      Exclusion_Step = exclusion_step,
      Manifest_Order = trait_index,
      Trait = trait_order[trait_index],
      Reason = reason,
      Failed_Pairs_At_Exclusion = as.integer(failed_degree),
      h2_liab = if (selected_scale == "liability") {
        self_lookup$h2[trait_index]
      } else {
        NA_real_
      },
      h2_liab_se = if (selected_scale == "liability") {
        self_lookup$h2_se[trait_index]
      } else {
        NA_real_
      },
      h2_Z = self_lookup$h2_Z[trait_index],
      Heritability_Scale = selected_scale,
      Heritability = self_lookup$h2[trait_index],
      Heritability_SE = self_lookup$h2_se[trait_index],
      h2_obs = if (selected_scale == "observed") {
        self_lookup$h2[trait_index]
      } else {
        NA_real_
      },
      h2_obs_se = if (selected_scale == "observed") {
        self_lookup$h2_se[trait_index]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    active[trait_index] <<- FALSE
  }

  initial_degree <- tabulate(
    c(failed_pairs$pair_i, failed_pairs$pair_j),
    nbins = number_traits
  ) - tabulate(
    failed_pairs$pair_i[failed_pairs$pair_i == failed_pairs$pair_j],
    nbins = number_traits
  )

  for (trait_index in sort(invalid_self_indices)) {
    reason <- self_pair_qc$Self_QC_Failure_Reason[trait_index]
    add_exclusion(trait_index, reason, initial_degree[trait_index])
  }

  remaining_failed <- failed_pairs[
    pair_i != pair_j & active[pair_i] & active[pair_j]
  ]
  while (nrow(remaining_failed) > 0L) {
    degree <- tabulate(
      c(remaining_failed$pair_i, remaining_failed$pair_j),
      nbins = number_traits
    )
    candidates <- which(active & degree == max(degree[active]))
    candidate_h2_z <- self_lookup$h2_Z[candidates]
    candidate_h2_z[!is.finite(candidate_h2_z)] <- -Inf
    candidates <- candidates[candidate_h2_z == min(candidate_h2_z)]
    selected_index <- max(candidates)

    add_exclusion(
      selected_index,
      "Greedy removal to resolve missing/non-finite LDSC pairs",
      degree[selected_index]
    )
    remaining_failed <- remaining_failed[
      pair_i != selected_index & pair_j != selected_index
    ]
  }

  retained_traits <- trait_order[active]
  if (length(retained_traits) < 2L) {
    stop(
      glue(
        "Resolving failed LDSC estimates would leave only ",
        "{length(retained_traits)} trait(s); at least two are required."
      ),
      call. = FALSE
    )
  }

  excluded_traits <- rbindlist(exclusions, use.names = TRUE, fill = TRUE)
  warning(
    glue(
      "--failed_ldsc_action=drop_traits excluded {nrow(excluded_traits)} ",
      "trait(s) to obtain a complete finite LDSC subset of ",
      "{length(retained_traits)} traits. No pair estimate was imputed."
    ),
    call. = FALSE
  )

  if (!is.null(source_rows_read)) {
    attr(ldsc_rows, "source_rows_read") <- source_rows_read
  }
  attr(ldsc_rows, "heritability_scale") <- selected_scale
  attr(ldsc_rows, "heritability_value_column") <- h2_source_column
  attr(ldsc_rows, "heritability_se_column") <- h2_se_source_column
  list(
    ldsc_rows = ldsc_rows,
    trait_order = retained_traits,
    excluded_traits = as.data.frame(excluded_traits),
    initial_failed_pair_count = initial_failed_pair_count,
    initial_self_qc_failure_count = length(self_qc_only),
    self_pair_qc = as.data.frame(self_pair_qc),
    heritability_scale = selected_scale
  )
}

coerce_and_validate_numeric <- function(ldsc_rows) {
  ldsc_rows <- coerce_python_ldsc_numeric(ldsc_rows)

  non_finite <- vapply(
    python_ldsc_numeric_columns,
    function(column_name) any(!is.finite(ldsc_rows[[column_name]])),
    logical(1)
  )
  if (any(non_finite)) {
    stop(
      "Selected Python LDSC rows contain missing or non-finite values in: ",
      paste(names(non_finite)[non_finite], collapse = ", "),
      call. = FALSE
    )
  }

  se_columns <- c("se", "h2_se", "h2_int_se", "gcov_int_se")
  non_positive_se <- vapply(
    se_columns,
    function(column_name) any(ldsc_rows[[column_name]] <= 0),
    logical(1)
  )
  if (any(non_positive_se)) {
    stop(
      "Selected Python LDSC rows contain non-positive standard errors in: ",
      paste(names(non_positive_se)[non_positive_se], collapse = ", "),
      call. = FALSE
    )
  }

  if (any(ldsc_rows$p < 0 | ldsc_rows$p > 1)) {
    stop("Selected Python LDSC p-values must lie in [0, 1].", call. = FALSE)
  }

  ldsc_rows
}

values_agree <- function(minimum, maximum, tolerance,
                         comparison_epsilon = 1e-12) {
  (maximum - minimum) <= tolerance + comparison_epsilon
}

validate_symmetric_matrix <- function(matrix_object, matrix_name,
                                      tolerance = 1e-8) {
  matrix_object <- as.matrix(matrix_object)

  if (nrow(matrix_object) != ncol(matrix_object)) {
    stop(glue("{matrix_name} must be square."), call. = FALSE)
  }
  if (any(!is.finite(matrix_object))) {
    stop(glue("{matrix_name} contains non-finite values."), call. = FALSE)
  }

  maximum_asymmetry <- max(abs(matrix_object - t(matrix_object)))
  if (maximum_asymmetry > tolerance) {
    stop(
      glue(
        "{matrix_name} is asymmetric; maximum absolute difference is ",
        "{format(maximum_asymmetry, scientific = TRUE)}."
      ),
      call. = FALSE
    )
  }

  (matrix_object + t(matrix_object)) / 2
}

canonicalize_and_validate_ldsc <- function(ldsc_rows, trait_order,
                                           tolerance = 1e-3,
                                           source_rows_read = nrow(ldsc_rows),
                                           duplicate_z_tolerance = 1e-2,
                                           self_rg_tolerance = 1e-2,
                                           comparison_epsilon = 1e-12,
                                           h2_z_warn_threshold = 2,
                                           z_consistency_tolerance = 1e-2,
                                           z_consistency_action = c("warn", "error"),
                                           rg_out_of_range_action = c("warn", "error"),
                                           heritability_scale = c(
                                             "auto", "liability", "observed"
                                           )) {
  z_consistency_action <- match.arg(z_consistency_action)
  rg_out_of_range_action <- match.arg(rg_out_of_range_action)
  heritability_scale <- match.arg(heritability_scale)
  numeric_tolerances <- c(
    tolerance = tolerance,
    duplicate_z_tolerance = duplicate_z_tolerance,
    self_rg_tolerance = self_rg_tolerance,
    z_consistency_tolerance = z_consistency_tolerance
  )
  if (any(!is.finite(numeric_tolerances) | numeric_tolerances <= 0)) {
    stop(
      "Duplicate, self-rg, and z-consistency tolerances must be finite and greater than zero.",
      call. = FALSE
    )
  }
  if (!is.finite(comparison_epsilon) || comparison_epsilon < 0) {
    stop("comparison_epsilon must be finite and non-negative.", call. = FALSE)
  }
  if (!is.finite(h2_z_warn_threshold) || h2_z_warn_threshold < 0) {
    stop("h2_z_warn_threshold must be finite and non-negative.", call. = FALSE)
  }

  ldsc_rows <- normalize_heritability_columns(
    ldsc_rows,
    heritability_scale
  )
  selected_scale <- attr(ldsc_rows, "heritability_scale", exact = TRUE)
  h2_source_column <- attr(
    ldsc_rows,
    "heritability_value_column",
    exact = TRUE
  )
  h2_se_source_column <- attr(
    ldsc_rows,
    "heritability_se_column",
    exact = TRUE
  )

  missing_columns <- setdiff(required_python_ldsc_columns, names(ldsc_rows))
  if (length(missing_columns) > 0L) {
    stop(
      "Python LDSC input is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ldsc_rows[, p1 := trimws(as.character(p1))]
  ldsc_rows[, p2 := trimws(as.character(p2))]
  ldsc_rows <- ldsc_rows[p1 %chin% trait_order & p2 %chin% trait_order]

  present_traits <- unique(c(ldsc_rows$p1, ldsc_rows$p2))
  missing_traits <- trait_order[!trait_order %chin% present_traits]
  if (length(missing_traits) > 0L) {
    stop(
      "Selected traits absent from the Python LDSC file:\n",
      paste(missing_traits, collapse = "\n"),
      call. = FALSE
    )
  }

  ldsc_rows <- coerce_and_validate_numeric(ldsc_rows)

  trait_index_1 <- match(ldsc_rows$p1, trait_order)
  trait_index_2 <- match(ldsc_rows$p2, trait_order)
  ldsc_rows[, pair_i := pmin(trait_index_1, trait_index_2)]
  ldsc_rows[, pair_j := pmax(trait_index_1, trait_index_2)]

  duplicate_stats <- ldsc_rows[, c(
    setNames(lapply(.SD, min), paste0("min_", names(.SD))),
    setNames(lapply(.SD, max), paste0("max_", names(.SD))),
    list(source_row_count = .N)
  ), by = .(pair_i, pair_j), .SDcols = duplicate_comparison_columns]

  conflict <- rep(FALSE, nrow(duplicate_stats))
  conflict_fields <- rep("", nrow(duplicate_stats))
  for (column_name in duplicate_comparison_columns) {
    current_tolerance <- if (column_name == "z") {
      duplicate_z_tolerance
    } else {
      tolerance
    }
    agrees <- values_agree(
      duplicate_stats[[paste0("min_", column_name)]],
      duplicate_stats[[paste0("max_", column_name)]],
      current_tolerance,
      comparison_epsilon
    )
    newly_conflicting <- !agrees
    conflict[newly_conflicting] <- TRUE
    conflict_fields[newly_conflicting] <- ifelse(
      conflict_fields[newly_conflicting] == "",
      column_name,
      paste0(conflict_fields[newly_conflicting], ",", column_name)
    )
  }

  if (any(conflict)) {
    bad <- duplicate_stats[conflict]
    descriptions <- paste0(
      trait_order[bad$pair_i], " <-> ", trait_order[bad$pair_j],
      " [", conflict_fields[conflict], "]"
    )
    stop(
      "Conflicting duplicate Python LDSC estimates were found:\n",
      paste(head(descriptions, 25L), collapse = "\n"),
      if (length(descriptions) > 25L) {
        glue("\n... and {length(descriptions) - 25L} more")
      } else {
        ""
      },
      call. = FALSE
    )
  }

  collapsed <- ldsc_rows[, lapply(.SD, mean),
    by = .(pair_i, pair_j),
    .SDcols = duplicate_comparison_columns
  ]
  collapsed <- merge(
    collapsed,
    duplicate_stats[, .(pair_i, pair_j, source_row_count)],
    by = c("pair_i", "pair_j"),
    sort = FALSE
  )
  setorder(collapsed, pair_i, pair_j)

  number_traits <- length(trait_order)
  expected_pairs <- number_traits * (number_traits + 1L) / 2L
  expected_grid <- CJ(
    pair_i = seq_len(number_traits),
    pair_j = seq_len(number_traits)
  )[pair_i <= pair_j]
  missing_pairs <- expected_grid[
    !collapsed,
    on = .(pair_i, pair_j)
  ]

  if (nrow(missing_pairs) > 0L) {
    descriptions <- paste0(
      trait_order[missing_pairs$pair_i], " <-> ",
      trait_order[missing_pairs$pair_j]
    )
    stop(
      glue(
        "Python LDSC input is incomplete: expected {expected_pairs} unique ",
        "self/pairwise combinations but found {nrow(collapsed)}. Missing:\n"
      ),
      paste(head(descriptions, 50L), collapse = "\n"),
      if (length(descriptions) > 50L) {
        glue("\n... and {length(descriptions) - 50L} more")
      } else {
        ""
      },
      call. = FALSE
    )
  }

  self_pair_qc <- evaluate_self_pair_qc(
    ldsc_rows,
    trait_order,
    self_rg_tolerance = self_rg_tolerance,
    comparison_epsilon = comparison_epsilon,
    h2_z_warn_threshold = h2_z_warn_threshold,
    heritability_scale = selected_scale
  )
  failed_self_qc <- self_pair_qc[!self_pair_qc$Self_QC_Pass, , drop = FALSE]
  if (nrow(failed_self_qc) > 0L) {
    descriptions <- paste0(
      failed_self_qc$Trait,
      " [",
      failed_self_qc$Self_QC_Failure_Reason,
      "]"
    )
    stop(
      "Selected traits failed LDSC self-pair QC:\n",
      paste(descriptions, collapse = "\n"),
      call. = FALSE
    )
  }

  if (h2_z_warn_threshold > 0 &&
      any(self_pair_qc$H2_Z_Below_Warning_Threshold)) {
    warning(
      glue(
        "{sum(self_pair_qc$H2_Z_Below_Warning_Threshold)} retained trait(s) ",
        "have self-pair h2/SE below {h2_z_warn_threshold}. This is a ",
        "diagnostic warning only; no trait was removed for low h2/SE."
      ),
      call. = FALSE
    )
  }

  self_summary <- data.table(
    pair_i = self_pair_qc$Manifest_Order,
    rg = self_pair_qc$Self_RG,
    se = self_pair_qc$Self_RG_SE,
    z = self_pair_qc$Self_Z,
    p = self_pair_qc$Self_P,
    h2 = self_pair_qc$Self_H2,
    h2_se = self_pair_qc$Self_H2_SE,
    h2_int = self_pair_qc$Self_H2_Intercept,
    h2_int_se = self_pair_qc$Self_H2_Intercept_SE,
    source_row_count = self_pair_qc$Self_Source_Rows
  )

  S_Stand <- diag(1, number_traits)
  I <- matrix(NA_real_, number_traits, number_traits)
  rg_se_matrix <- matrix(NA_real_, number_traits, number_traits)
  intercept_se_matrix <- matrix(NA_real_, number_traits, number_traits)
  matrix_names <- list(trait_order, trait_order)
  dimnames(S_Stand) <- matrix_names
  dimnames(I) <- matrix_names
  dimnames(rg_se_matrix) <- matrix_names
  dimnames(intercept_se_matrix) <- matrix_names

  diag(I) <- self_summary$h2_int[match(seq_len(number_traits), self_summary$pair_i)]
  diag(rg_se_matrix) <- self_summary$se[
    match(seq_len(number_traits), self_summary$pair_i)
  ]
  diag(intercept_se_matrix) <- self_summary$h2_int_se[
    match(seq_len(number_traits), self_summary$pair_i)
  ]

  off_diagonal <- collapsed[pair_i < pair_j]
  upper_indices <- cbind(off_diagonal$pair_i, off_diagonal$pair_j)
  lower_indices <- cbind(off_diagonal$pair_j, off_diagonal$pair_i)

  S_Stand[upper_indices] <- S_Stand[lower_indices] <- off_diagonal$rg
  I[upper_indices] <- I[lower_indices] <- off_diagonal$gcov_int
  rg_se_matrix[upper_indices] <- rg_se_matrix[lower_indices] <- off_diagonal$se
  intercept_se_matrix[upper_indices] <-
    intercept_se_matrix[lower_indices] <- off_diagonal$gcov_int_se

  S_Stand <- validate_symmetric_matrix(S_Stand, "S_Stand")
  I <- validate_symmetric_matrix(I, "I")
  rg_se_matrix <- validate_symmetric_matrix(
    rg_se_matrix,
    "Python LDSC rg SE matrix"
  )
  intercept_se_matrix <- validate_symmetric_matrix(
    intercept_se_matrix,
    "Python LDSC intercept SE matrix"
  )

  intercept_eigenvalues <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  chol_error <- tryCatch(
    {
      chol(I)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(chol_error)) {
    stop(
      glue(
        "The LDSC intercept matrix I is not positive definite ",
        "(minimum eigenvalue {format(min(intercept_eigenvalues), scientific = TRUE)}): ",
        "{chol_error}"
      ),
      call. = FALSE
    )
  }

  out_of_range <- off_diagonal$rg < -1 | off_diagonal$rg > 1
  if (any(out_of_range)) {
    out_of_range_message <- glue(
      "{sum(out_of_range)} off-diagonal genetic correlation(s) lie outside ",
      "[-1, 1]. Values were not clamped and are flagged in the diagnostics."
    )
    if (rg_out_of_range_action == "error") {
      stop(out_of_range_message, call. = FALSE)
    }
    warning(out_of_range_message, call. = FALSE)
  }

  expected_z <- off_diagonal$rg / off_diagonal$se
  z_relative_difference <- abs(off_diagonal$z - expected_z) /
    pmax(1, abs(off_diagonal$z), abs(expected_z))
  z_inconsistent <- z_relative_difference >
    z_consistency_tolerance + comparison_epsilon
  if (any(z_inconsistent)) {
    z_message <- glue(
      "{sum(z_inconsistent)} off-diagonal z value(s) differ from rg/SE ",
      "beyond the relative tolerance {z_consistency_tolerance}. This is a ",
      "diagnostic consistency check; PCA values were not filtered or reweighted."
    )
    if (z_consistency_action == "error") {
      stop(z_message, call. = FALSE)
    }
    warning(z_message, call. = FALSE)
  }

  heritability_results <- data.frame(
    Trait = trait_order,
    h2_liab = if (selected_scale == "liability") {
      self_summary$h2
    } else {
      NA_real_
    },
    h2_liab_se = if (selected_scale == "liability") {
      self_summary$h2_se
    } else {
      NA_real_
    },
    Z = self_summary$h2 / self_summary$h2_se,
    P = 2 * pnorm(
      abs(self_summary$h2 / self_summary$h2_se),
      lower.tail = FALSE
    ),
    h2_int = self_summary$h2_int,
    h2_int_se = self_summary$h2_int_se,
    Heritability_Scale = selected_scale,
    Heritability = self_summary$h2,
    Heritability_SE = self_summary$h2_se,
    h2_obs = if (selected_scale == "observed") {
      self_summary$h2
    } else {
      NA_real_
    },
    h2_obs_se = if (selected_scale == "observed") {
      self_summary$h2_se
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  genetic_correlation_results <- data.frame(
    Trait_1 = trait_order[off_diagonal$pair_i],
    Trait_2 = trait_order[off_diagonal$pair_j],
    rg = off_diagonal$rg,
    SE = off_diagonal$se,
    Z = off_diagonal$z,
    Z_From_RG_SE = expected_z,
    Z_Relative_Difference = z_relative_difference,
    Z_Inconsistent_With_RG_SE = z_inconsistent,
    P = off_diagonal$p,
    gcov_int = off_diagonal$gcov_int,
    gcov_int_se = off_diagonal$gcov_int_se,
    Source_Rows_Collapsed = off_diagonal$source_row_count,
    RG_Outside_Unit_Interval = out_of_range,
    stringsAsFactors = FALSE
  )

  pairs_per_trait <- tabulate(
    c(collapsed$pair_i, collapsed$pair_j),
    nbins = number_traits
  ) - tabulate(
    collapsed$pair_i[collapsed$pair_i == collapsed$pair_j],
    nbins = number_traits
  )
  validation_summary <- data.frame(
    Trait_Order = seq_len(number_traits),
    Trait = trait_order,
    Present_In_Python_LDSC = TRUE,
    Self_Pair_Found = self_pair_qc$Self_Pair_Found,
    Self_QC_Pass = self_pair_qc$Self_QC_Pass,
    Self_RG = self_pair_qc$Self_RG,
    Self_RG_Deviation = self_pair_qc$Self_RG_Deviation,
    Self_H2 = self_pair_qc$Self_H2,
    Self_H2_SE = self_pair_qc$Self_H2_SE,
    Heritability_Scale_Used = selected_scale,
    Heritability_Source_Column = h2_source_column,
    Heritability_SE_Source_Column = h2_se_source_column,
    Self_H2_Z = self_pair_qc$H2_Z,
    Low_H2_Z_Diagnostic = self_pair_qc$H2_Z_Below_Warning_Threshold,
    Unique_Pairs_With_Selected_Traits = pairs_per_trait,
    Expected_Pairs_With_Selected_Traits = number_traits,
    Total_Source_Rows_Read = source_rows_read,
    Selected_Source_Rows = nrow(ldsc_rows),
    Expected_Unique_Self_And_Pairwise_Rows = expected_pairs,
    Observed_Unique_Self_And_Pairwise_Rows = nrow(collapsed),
    Duplicate_Rows_Collapsed = nrow(ldsc_rows) - nrow(collapsed),
    Duplicate_Tolerance = tolerance,
    Duplicate_Z_Tolerance = duplicate_z_tolerance,
    Self_RG_Tolerance = self_rg_tolerance,
    Comparison_Epsilon = comparison_epsilon,
    H2_Z_Warn_Threshold = h2_z_warn_threshold,
    Low_H2_Z_Diagnostic_Count = sum(
      self_pair_qc$H2_Z_Below_Warning_Threshold
    ),
    Z_Consistency_Tolerance = z_consistency_tolerance,
    Z_Consistency_Action = z_consistency_action,
    Z_Inconsistent_With_RG_SE_Count = sum(z_inconsistent),
    RG_Out_Of_Range_Action = rg_out_of_range_action,
    Off_Diagonal_RG_Outside_Unit_Interval = sum(out_of_range),
    Minimum_Intercept_Matrix_Eigenvalue = min(intercept_eigenvalues),
    stringsAsFactors = FALSE
  )

  list(
    heritability_scale = selected_scale,
    heritability_source_column = h2_source_column,
    heritability_se_source_column = h2_se_source_column,
    S_Stand = S_Stand,
    I = I,
    rg_se_matrix = rg_se_matrix,
    intercept_se_matrix = intercept_se_matrix,
    collapsed_pairs = collapsed,
    heritability_results = heritability_results,
    genetic_correlation_results = genetic_correlation_results,
    validation_summary = validation_summary,
    self_pair_qc = self_pair_qc
  )
}

# Reconstruct the unstandardized genetic covariance matrix required by the
# tutorial's alternative procedure. Python LDSC supplies rg and univariate h2
# here, but not GenomicSEM's full sampling-covariance V matrices; none are
# fabricated. Trait names and order must agree exactly.
construct_genetic_covariance_matrix <- function(correlation_matrix,
                                                heritability,
                                                trait_order,
                                                heritability_scale = c(
                                                  "liability", "observed"
                                                ),
                                                tolerance = 1e-10) {
  heritability_scale <- match.arg(heritability_scale)
  correlation_matrix <- validate_symmetric_matrix(
    correlation_matrix,
    "genetic correlation matrix"
  )
  if (!identical(rownames(correlation_matrix), trait_order) ||
      !identical(colnames(correlation_matrix), trait_order)) {
    stop(
      "The phenotype order in the genetic correlation matrix differs from the manifest.",
      call. = FALSE
    )
  }
  if (!is.finite(tolerance) || tolerance <= 0) {
    stop("Covariance reconstruction tolerance must be positive.", call. = FALSE)
  }

  if (!is.null(names(heritability)) &&
      !identical(names(heritability), trait_order)) {
    stop(
      "The named heritability vector order differs from the manifest.",
      call. = FALSE
    )
  }
  heritability <- as.numeric(heritability)
  if (length(heritability) != length(trait_order)) {
    stop("The heritability vector length differs from the manifest.", call. = FALSE)
  }
  if (any(!is.finite(heritability)) || any(heritability <= 0)) {
    stop(
      glue(
        "Covariance PCA requires finite positive self-pair ",
        "{heritability_column_sets[[heritability_scale]][1L]} values."
      ),
      call. = FALSE
    )
  }

  h2_scale <- sqrt(heritability)
  covariance_matrix <- correlation_matrix * tcrossprod(h2_scale)
  diag(covariance_matrix) <- heritability
  dimnames(covariance_matrix) <- list(trait_order, trait_order)
  covariance_matrix <- validate_symmetric_matrix(
    covariance_matrix,
    "derived genetic covariance matrix"
  )

  reconstructed_correlation <- covariance_matrix / tcrossprod(h2_scale)
  diag(reconstructed_correlation) <- 1
  maximum_reconstruction_error <- max(
    abs(reconstructed_correlation - correlation_matrix)
  )
  if (maximum_reconstruction_error > tolerance) {
    stop(
      glue(
        "Derived genetic covariance matrix failed its reconstruction check; ",
        "maximum correlation error was ",
        "{format(maximum_reconstruction_error, scientific = TRUE)}."
      ),
      call. = FALSE
    )
  }

  h2_column <- heritability_column_sets[[heritability_scale]][1L]
  attr(covariance_matrix, "source") <- paste0(
    "rg * sqrt(self ", h2_column, " trait 1 * self ",
    h2_column, " trait 2)"
  )
  attr(covariance_matrix, "heritability_scale") <- heritability_scale
  attr(covariance_matrix, "maximum_reconstruction_error") <-
    maximum_reconstruction_error
  covariance_matrix
}

make_genetic_covariance_results <- function(correlation_results,
                                            covariance_matrix,
                                            trait_order,
                                            heritability_scale = c(
                                              "liability", "observed"
                                            )) {
  heritability_scale <- match.arg(heritability_scale)
  correlation_results <- as.data.frame(correlation_results)
  trait_index_1 <- match(correlation_results$Trait_1, trait_order)
  trait_index_2 <- match(correlation_results$Trait_2, trait_order)
  if (anyNA(trait_index_1) || anyNA(trait_index_2)) {
    stop(
      "Genetic-correlation result order cannot be matched to covariance traits.",
      call. = FALSE
    )
  }

  covariance_results <- correlation_results[, c(
    "Trait_1", "Trait_2", "rg", "SE", "Z", "P"
  ), drop = FALSE]
  covariance_results$Heritability_Scale <- heritability_scale
  covariance_results$Heritability_1 <-
    diag(covariance_matrix)[trait_index_1]
  covariance_results$Heritability_2 <-
    diag(covariance_matrix)[trait_index_2]
  covariance_results$h2_liab_1 <- if (heritability_scale == "liability") {
    covariance_results$Heritability_1
  } else {
    NA_real_
  }
  covariance_results$h2_liab_2 <- if (heritability_scale == "liability") {
    covariance_results$Heritability_2
  } else {
    NA_real_
  }
  covariance_results$h2_obs_1 <- if (heritability_scale == "observed") {
    covariance_results$Heritability_1
  } else {
    NA_real_
  }
  covariance_results$h2_obs_2 <- if (heritability_scale == "observed") {
    covariance_results$Heritability_2
  } else {
    NA_real_
  }
  covariance_results$Genetic_Covariance_Derived <- covariance_matrix[
    cbind(trait_index_1, trait_index_2)
  ]
  h2_column <- heritability_column_sets[[heritability_scale]][1L]
  covariance_results$Derivation <- paste0(
    "rg * sqrt(self ", h2_column, "_1 * self ", h2_column, "_2)"
  )
  covariance_results
}

# Orient PC1 either with the Fuertjes tutorial convention (the default) or
# preserve the direction returned by eigen(). Both always operate on the whole
# loading vector; individual loadings are never selectively changed.
orient_pc1_loadings <- function(loadings,
                                method = c("tutorial", "as_computed")) {
  method <- match.arg(method)
  loadings <- as.numeric(loadings)
  if (length(loadings) == 0L || any(!is.finite(loadings))) {
    stop("PC1 loadings must be a non-empty finite numeric vector.", call. = FALSE)
  }

  median_before_orientation <- median(loadings)
  sign_multiplier <- if (
    method == "tutorial" && median_before_orientation < 0
  ) -1 else 1
  oriented_loadings <- loadings * sign_multiplier

  list(
    method = method,
    loadings = oriented_loadings,
    loadings_before_orientation = loadings,
    median_before_orientation = median_before_orientation,
    median_after_orientation = median(oriented_loadings),
    sign_multiplier = sign_multiplier,
    flip_applied = sign_multiplier == -1
  )
}

compute_pc1 <- function(pca_matrix, trait_order,
                        negative_eigen_action = c("warn", "error"),
                        matrix_eigen_tolerance = 1e-8,
                        pc1_orientation = c("tutorial", "as_computed"),
                        pca_matrix_type = c("correlation", "covariance")) {
  negative_eigen_action <- match.arg(negative_eigen_action)
  pc1_orientation <- match.arg(pc1_orientation)
  pca_matrix_type <- match.arg(pca_matrix_type)
  if (!is.finite(matrix_eigen_tolerance) || matrix_eigen_tolerance <= 0) {
    stop(
      "matrix_eigen_tolerance must be finite and greater than zero.",
      call. = FALSE
    )
  }

  pca_matrix <- validate_symmetric_matrix(pca_matrix, "PCA input matrix")

  if (!identical(rownames(pca_matrix), trait_order) ||
      !identical(colnames(pca_matrix), trait_order)) {
    stop("The phenotype order in the PCA matrix differs from the manifest.", call. = FALSE)
  }

  eigen_selected <- eigen(pca_matrix, symmetric = TRUE)
  eigenvectors <- eigen_selected$vectors
  eigenvalues <- eigen_selected$values

  if (!is.finite(eigenvalues[1L]) || eigenvalues[1L] <= 0) {
    stop(
      glue(
        "PC1 requires a finite positive first eigenvalue; observed ",
        "{eigenvalues[1L]}."
      ),
      call. = FALSE
    )
  }

  eigenvalues_after_pmax <- pmax(eigenvalues, 0)
  loadings_before_orientation <- as.vector(
    eigenvectors %*% sqrt(diag(eigenvalues_after_pmax))[, 1]
  )

  if (length(loadings_before_orientation) != length(trait_order) ||
      any(!is.finite(loadings_before_orientation))) {
    stop("PC1 loadings are invalid.", call. = FALSE)
  }
  if (all(abs(loadings_before_orientation) <= sqrt(.Machine$double.eps))) {
    stop("All PC1 loadings are effectively zero.", call. = FALSE)
  }

  orientation <- orient_pc1_loadings(
    loadings_before_orientation,
    method = pc1_orientation
  )
  loadings <- orientation$loadings
  eigenvectors_before_orientation <- eigenvectors
  eigenvectors[, 1L] <- eigenvectors[, 1L] * orientation$sign_multiplier

  negative <- eigenvalues < 0
  eigenvalue_scale <- max(1, max(abs(eigenvalues)))
  negative_eigenvalue_threshold <- matrix_eigen_tolerance * eigenvalue_scale
  substantive_negative <- eigenvalues < -negative_eigenvalue_threshold
  if (any(negative)) {
    eigen_message <- glue(
      "The {pca_matrix_type} PCA matrix has {sum(negative)} negative ",
      "eigenvalue(s), including {sum(substantive_negative)} below the ",
      "numerical threshold -{format(negative_eigenvalue_threshold, scientific = TRUE)}. ",
      "pmax(value, 0) applies only to loading calculation; the matrix itself ",
      "was not altered."
    )
    if (negative_eigen_action == "error" && any(substantive_negative)) {
      stop(eigen_message, call. = FALSE)
    }
    warning(eigen_message, call. = FALSE)
  }

  eigenvalue_results <- data.frame(
    PC = seq_along(eigenvalues),
    PCA_Matrix_Type = pca_matrix_type,
    Eigenvalue_Raw = eigenvalues,
    Eigenvalue_After_pmax = eigenvalues_after_pmax,
    Pmax_Guard_Applied = negative,
    Substantive_Negative_Eigenvalue = substantive_negative,
    Negative_Eigenvalue_Threshold = negative_eigenvalue_threshold,
    stringsAsFactors = FALSE
  )

  raw_total <- sum(eigenvalues)
  guarded_total <- sum(eigenvalues_after_pmax)
  loading_results <- data.frame(
    Trait = trait_order,
    PCA_Matrix_Type = pca_matrix_type,
    PC1_Eigenvector = eigenvectors[, 1L],
    PC1_Loading_Used = loadings,
    PC1_Standardised_Loading = if (pca_matrix_type == "correlation") {
      loadings
    } else {
      NA_real_
    },
    PC1_Correlation_Matrix_Loading = if (pca_matrix_type == "correlation") {
      loadings
    } else {
      NA_real_
    },
    PC1_Covariance_Matrix_Loading = if (pca_matrix_type == "covariance") {
      loadings
    } else {
      NA_real_
    },
    PC1_Eigenvector_Before_Orientation =
      eigenvectors_before_orientation[, 1L],
    PC1_Loading_Before_Orientation =
      orientation$loadings_before_orientation,
    PC1_Loading_Median_Before_Orientation =
      orientation$median_before_orientation,
    PC1_Loading_Median_After_Orientation =
      orientation$median_after_orientation,
    PC1_Sign_Method = orientation$method,
    PC1_Sign_Multiplier = orientation$sign_multiplier,
    PC1_Median_Flip_Applied = orientation$flip_applied,
    PC1_Eigenvalue_Raw = eigenvalues[1L],
    PC1_Eigenvalue_After_pmax = eigenvalues_after_pmax[1L],
    PC1_Variance_Explained = eigenvalues[1L] / raw_total,
    PC1_Variance_Explained_After_pmax =
      eigenvalues_after_pmax[1L] / guarded_total,
    Any_Nonpositive_Eigenvalue_Guarded = any(negative),
    Substantive_Negative_Eigenvalue_Count = sum(substantive_negative),
    Negative_Eigen_Action = negative_eigen_action,
    Matrix_Eigen_Tolerance = matrix_eigen_tolerance,
    stringsAsFactors = FALSE
  )

  list(
    pca_matrix_type = pca_matrix_type,
    loadings = loadings,
    eigenvectors = eigenvectors,
    eigenvalues = eigenvalues,
    eigenvalues_after_pmax = eigenvalues_after_pmax,
    orientation = orientation,
    loading_results = loading_results,
    eigenvalue_results = eigenvalue_results
  )
}

canonicalize_gwama_column_names <- function(column_names) {
  column_names[column_names == "A1"] <- "EA"
  column_names[column_names == "A2"] <- "OA"
  column_names[column_names == "p"] <- "P"
  column_names
}

validate_gwama_columns <- function(column_names, file_path, trait_name) {
  canonical_names <- canonicalize_gwama_column_names(column_names)

  if (length(canonical_names) != length(gwama_required_columns) ||
      !identical(canonical_names, gwama_required_columns)) {
    stop(
      glue(
        "GWAMA input for '{trait_name}' must contain exactly these columns ",
        "in order: {paste(gwama_required_columns, collapse = ', ')}. ",
        "Observed after allowed A1/A2/p renaming: ",
        "{paste(canonical_names, collapse = ', ')}. File: {file_path}"
      ),
      call. = FALSE
    )
  }

  canonical_names
}

build_gwama_input_paths <- function(trait_order, gpca_input_folder,
                                    splitby_chr, chr = NULL) {
  if (splitby_chr == "split") {
    if (is.null(chr) || length(chr) != 1L || is.na(chr)) {
      stop("A chromosome is required when --splitby_chr=split.", call. = FALSE)
    }
    paths <- file.path(
      gpca_input_folder,
      paste0(trait_order, "_chr", chr, "_GenomicPCA_inputs.tsv")
    )
  } else if (splitby_chr == "nosplit") {
    paths <- file.path(
      gpca_input_folder,
      paste0(trait_order, "_GenomicPCA_inputs.tsv")
    )
  } else {
    stop("splitby_chr must be 'split' or 'nosplit'.", call. = FALSE)
  }

  setNames(paths, trait_order)
}

preflight_gwama_inputs <- function(trait_order, gpca_input_folder,
                                   splitby_chr) {
  chromosomes <- if (splitby_chr == "split") as.list(1:22) else list(NULL)
  all_paths <- unlist(lapply(
    chromosomes,
    function(chr) build_gwama_input_paths(
      trait_order,
      gpca_input_folder,
      splitby_chr,
      chr
    )
  ), use.names = FALSE)

  missing_paths <- all_paths[!file.exists(all_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing required GWAMA input files:\n",
      paste(head(missing_paths, 100L), collapse = "\n"),
      if (length(missing_paths) > 100L) {
        glue("\n... and {length(missing_paths) - 100L} more")
      } else {
        ""
      },
      call. = FALSE
    )
  }

  for (chr in chromosomes) {
    named_paths <- build_gwama_input_paths(
      trait_order,
      gpca_input_folder,
      splitby_chr,
      chr
    )
    for (trait_name in trait_order) {
      file_path <- named_paths[[trait_name]]
      header <- fread(
        file_path,
        nrows = 0L,
        data.table = FALSE,
        check.names = FALSE,
        showProgress = FALSE
      )
      validate_gwama_columns(names(header), file_path, trait_name)
    }
  }

  invisible(TRUE)
}

load_modified_gwama <- function(source_path) {
  gwama_environment <- new.env(parent = globalenv())
  sys.source(source_path, envir = gwama_environment)

  if (exists("my_GWAMA", envir = gwama_environment, mode = "function",
             inherits = FALSE)) {
    get("my_GWAMA", envir = gwama_environment, inherits = FALSE)
  } else if (exists("multivariate_GWAMA", envir = gwama_environment,
                    mode = "function", inherits = FALSE)) {
    get("multivariate_GWAMA", envir = gwama_environment, inherits = FALSE)
  } else {
    stop(
      paste0(
        "--source_path must define my_GWAMA() or multivariate_GWAMA(). ",
        "Supply the already-modified Fuertjes function."
      ),
      call. = FALSE
    )
  }
}

read_gwama_data <- function(trait_order, gpca_input_folder,
                            splitby_chr, chr = NULL) {
  input_paths <- build_gwama_input_paths(
    trait_order,
    gpca_input_folder,
    splitby_chr,
    chr
  )
  dat <- vector("list", length(trait_order))
  names(dat) <- trait_order

  for (trait_name in trait_order) {
    file_path <- input_paths[[trait_name]]
    trait_data <- fread(file_path, data.table = FALSE, check.names = FALSE)
    canonical_names <- validate_gwama_columns(
      names(trait_data),
      file_path,
      trait_name
    )
    names(trait_data) <- canonical_names
    dat[[trait_name]] <- trait_data
  }

  dat
}

assert_analysis_order <- function(dat, CTI, pca_matrix, loadings, trait_order) {
  if (!identical(names(dat), trait_order)) {
    stop("The phenotype order in dat differs from the manifest.", call. = FALSE)
  }
  if (!identical(rownames(CTI), trait_order) ||
      !identical(colnames(CTI), trait_order)) {
    stop("The phenotype order in CTI differs from names(dat).", call. = FALSE)
  }
  if (!identical(rownames(pca_matrix), trait_order) ||
      !identical(colnames(pca_matrix), trait_order)) {
    stop("The phenotype order in the PCA matrix differs from names(dat).", call. = FALSE)
  }
  if (length(loadings) != length(trait_order)) {
    stop("The PC1 loading order/length differs from names(dat).", call. = FALSE)
  }

  invisible(TRUE)
}

run_genomic_pca_gwama <- function(chr = NULL, trait_order, CTI, pca_matrix,
                                  loadings, gpca_input_folder, outdir,
                                  splitby_chr, my_GWAMA,
                                  pca_matrix_type = c("correlation", "covariance")) {
  pca_matrix_type <- match.arg(pca_matrix_type)
  dat <- read_gwama_data(
    trait_order,
    gpca_input_folder,
    splitby_chr,
    chr
  )
  assert_analysis_order(dat, CTI, pca_matrix, loadings, trait_order)

  analysis_name <- if (pca_matrix_type == "correlation") {
    "GenomicPCA_PC1"
  } else {
    "GenomicPCA_Covariance_PC1"
  }
  output_name <- if (splitby_chr == "split") {
    glue("Chr{chr}_{analysis_name}")
  } else {
    analysis_name
  }

  # Fuertjes tutorial Step C.3. The externally supplied function must already
  # implement W = sqrt(N) * the selected tutorial PC loading. The parameter is
  # still named
  # h2 only because that is the original GWAMA function's API.
  my_GWAMA(
    x = dat,
    cov_Z = CTI,
    h2 = loadings,
    out = outdir,
    name = output_name,
    output_gz = TRUE,
    check_columns = FALSE
  )

  invisible(output_name)
}

safe_gwama_worker <- function(chr = NULL, ...) {
  tryCatch(
    {
      output_name <- run_genomic_pca_gwama(chr = chr, ...)
      list(
        success = TRUE,
        chromosome = chr,
        output = output_name,
        error = NA_character_
      )
    },
    error = function(e) {
      list(
        success = FALSE,
        chromosome = chr,
        output = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
}

write_analysis_outputs <- function(outdir, ldsc, pc1, correlation_matrix,
                                   pca_matrix, CTI,
                                   missing_traits, failed_traits,
                                   self_pair_qc,
                                   covariance_matrix = NULL,
                                   genetic_covariance_results = NULL) {
  write.csv(
    ldsc$validation_summary,
    file.path(outdir, "Python_LDSC_Input_Validation_Summary.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$heritability_results,
    file.path(outdir, "Global_Heritability_Estimates.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$genetic_correlation_results,
    file.path(outdir, "Global_Genetic_Correlations.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$rg_se_matrix,
    file.path(outdir, "Python_LDSC_RG_SE_Matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    ldsc$intercept_se_matrix,
    file.path(outdir, "Python_LDSC_Intercept_SE_Matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    CTI,
    file.path(outdir, "GenomicPCA_CTI_Used.csv"),
    row.names = TRUE
  )
  write.csv(
    correlation_matrix,
    file.path(outdir, "GenomicPCA_Correlation_Matrix_Used.csv"),
    row.names = TRUE
  )
  write.csv(
    pca_matrix,
    file.path(outdir, "GenomicPCA_PCA_Matrix_Used.csv"),
    row.names = TRUE
  )
  if (!is.null(covariance_matrix)) {
    write.csv(
      covariance_matrix,
      file.path(outdir, "GenomicPCA_Covariance_Matrix_Used.csv"),
      row.names = TRUE
    )
  }
  if (!is.null(genetic_covariance_results)) {
    write.csv(
      genetic_covariance_results,
      file.path(outdir, "Global_Genetic_Covariances_Derived.csv"),
      row.names = FALSE
    )
  }
  write.csv(
    pc1$loading_results,
    file.path(outdir, "GenomicPCA_PC1_Weights_Used.csv"),
    row.names = FALSE
  )
  write.csv(
    pc1$eigenvalue_results,
    file.path(outdir, "GenomicPCA_Selected_Traits_Eigenvalues.csv"),
    row.names = FALSE
  )
  write.csv(
    missing_traits,
    file.path(outdir, "Python_LDSC_Dropped_Missing_Traits.csv"),
    row.names = FALSE
  )
  write.csv(
    failed_traits,
    file.path(outdir, "Python_LDSC_Dropped_Failed_Traits.csv"),
    row.names = FALSE
  )
  write.csv(
    self_pair_qc,
    file.path(outdir, "Python_LDSC_Self_Pair_QC.csv"),
    row.names = FALSE
  )
}

make_run_status <- function(gwama_results) {
  data.frame(
    Chromosome = vapply(
      gwama_results,
      function(result) {
        if (is.null(result$chromosome)) "whole_genome" else as.character(result$chromosome)
      },
      character(1)
    ),
    Success = vapply(
      gwama_results,
      function(result) isTRUE(result$success),
      logical(1)
    ),
    Output = vapply(
      gwama_results,
      function(result) as.character(result$output),
      character(1)
    ),
    Error = vapply(
      gwama_results,
      function(result) as.character(result$error),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

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
    heritability_scale = args$heritability_scale
  )
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
      heritability_scale = ldsc$heritability_scale
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
    pca_matrix_type = args$pca_matrix
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

if (sys.nframe() == 0L) {
  gpsca_main()
}
