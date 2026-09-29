# Managed GPCA input validation; summary-statistic/GWAMA headers are not modified here.
check_gpca_cli_options <- function(arguments) {
  groups <- list("--input", "--outdir", "--ldsc_results", "--n_cores", "--duplicate_tolerance")
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
  manifest <- fread(path, data.table = FALSE, check.names = FALSE, colClasses = "character")
  if (anyDuplicated(names(manifest))) stop("Manifest headers must be unique.", call. = FALSE)
  removed <- intersect(names(manifest), c("gwas_name", "sampleprevalence", "pop_prevalence",
    "populationprevalence", "munge_inputs", "traits"))
  if (length(removed)) stop("Unsupported manifest headers: ", paste(removed, collapse = ", "),
    ". Use the column names listed in --help.", call. = FALSE)
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
