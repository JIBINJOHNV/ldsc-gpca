# Managed GPCA input validation; summary-statistic/GWAMA headers are not modified here.
check_gpca_cli_options <- function(arguments) {
  groups <- list("--input", "--outdir", "--ldsc_results", "--n_cores", "--duplicate_tolerance", "--rg_normalization")
  for (group in groups) {
    values <- character()
    flags <- character()
    for (i in seq_along(arguments)) {
      flag <- sub("=.*$", "", arguments[i])
      if (!flag %in% group) next
      if (grepl("=", arguments[i], fixed = TRUE)) {
        value <- sub("^[^=]*=", "", arguments[i])
      } else if (i < length(arguments)) value <- arguments[i + 1L] else next
      values <- c(values, value)
      flags <- c(flags, flag)
    }
    if (group[1] %in% c("--n_cores", "--duplicate_tolerance")) {
      numbers <- suppressWarnings(as.numeric(values))
      if (all(is.finite(numbers))) values <- as.character(numbers)
    }
    if (length(unique(values)) > 1L)
      stop("Conflicting values for ", paste(unique(flags), collapse = " and "), call. = FALSE)
  }
  arguments
}

read_gpca_manifest <- function(path) {
  # Missing-value tokens are valid text in trait names and annotations.
  manifest <- fread(path, data.table = FALSE, check.names = FALSE,
                    colClasses = "character", na.strings = NULL)
  if (anyDuplicated(names(manifest))) stop("Manifest headers must be unique.", call. = FALSE)
  # Required canonical columns are checked by the consuming workflow.
  # Extra columns are annotations, not aliases for missing required fields.
  manifest
}

# Shared trait-name validation.

read_trait_manifest <- function(path) {
  manifest <- read_gpca_manifest(path)

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

# Shared CLI contract for both matrix-validation workflows.
add_matrix_validation_options <- function(parser) {
  group <- parser$add_argument_group("PC1 and CTI validation")
  group$add_argument("--pc1_only", action = "store_true", default = FALSE,
    help = "Save PC1 and both validation reports without GWAMA. CTI failure is nonblocking; invalid PC1 still stops. Default: off.")
  group$add_argument("--cti_action", choices = c("error", "explore_drop"), default = "error",
    help = "CTI failure policy: error stops combined validation/GWAMA; explore_drop tries audited eigenvector-based trait exclusions and rechecks BOTH matrices. No matrix repair. Default: error.")
  group$add_argument("--max_cti_drop_fraction", type = "double", default = 0, metavar = "FLOAT",
    help = "Maximum fraction of traits entering CTI assessment to exclude (0 <= value < 1). explore_drop requires an explicit positive limit. Not a biological cutoff. Default: 0.")
}

validate_matrix_options <- function(args) {
  if (is.null(args$pc1_only)) args$pc1_only <- FALSE
  if (is.null(args$cti_action)) args$cti_action <- "error"
  if (is.null(args$max_cti_drop_fraction)) args$max_cti_drop_fraction <- 0
  f <- args$max_cti_drop_fraction
  if (!is.finite(f) || f < 0 || f >= 1) stop("--max_cti_drop_fraction must be finite and in [0,1).", call. = FALSE)
  if (isTRUE(args$pc1_only) && isTRUE(args$validate_only)) stop("Choose --pc1_only or --validate_only, not both.", call. = FALSE)
  if (args$cti_action == "explore_drop" && (f <= 0 || isTRUE(args$pc1_only)))
    stop("--cti_action explore_drop requires a positive --max_cti_drop_fraction and cannot be used with --pc1_only.", call. = FALSE)
  if (args$cti_action == "error" && f != 0) stop("--max_cti_drop_fraction requires --cti_action explore_drop.", call. = FALSE)
  args
}
