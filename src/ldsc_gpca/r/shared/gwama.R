gwama_required_columns <- c("SNPID", "CHR", "BP", "EA", "OA", "EAF", "N", "Z", "P")

# Read required fields in canonical order before calling the modified GWAMA function.

canonicalize_gwama_column_names <- function(column_names) {
  column_names[column_names == "A1"] <- "EA"
  column_names[column_names == "A2"] <- "OA"
  column_names[column_names == "p"] <- "P"
  column_names
}

validate_gwama_columns <- function(column_names, file_path, trait_name) {
  canonical_names <- canonicalize_gwama_column_names(column_names)

  # Prefer an explicit canonical name when both it and a legacy alias exist.
  # This allows A1/A2/p to remain harmless extra annotations beside EA/OA/P.
  for (canonical in c("EA", "OA", "P")) {
    alias <- c(EA = "A1", OA = "A2", P = "p")[[canonical]]
    if (canonical %in% column_names) canonical_names[column_names == alias] <- alias
  }
  missing <- setdiff(gwama_required_columns, canonical_names)
  ambiguous <- gwama_required_columns[vapply(gwama_required_columns,
    function(column) sum(canonical_names == column) > 1L, logical(1))]
  if (length(missing) || length(ambiguous)) {
    stop(
      glue(
        "GWAMA input for '{trait_name}' is missing required columns: ",
        "{paste(missing, collapse = ', ')}; ambiguous required columns: ",
        "{paste(ambiguous, collapse = ', ')}. ",
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
    dat[[trait_name]] <- trait_data[, gwama_required_columns, drop = FALSE]
  }

  dat
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
  output_file <- file.path(outdir, paste0(output_name, ".N_weighted_GWAMA.results.txt.gz"))
  previous <- file.info(output_file)[, c("size", "mtime", "ctime")]

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

  current <- file.info(output_file)[, c("size", "mtime", "ctime")]
  if (is.na(current$size) || current$size <= 0 || identical(previous, current))
    stop("GWAMA produced no fresh non-empty result: ", output_file, call. = FALSE)

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
