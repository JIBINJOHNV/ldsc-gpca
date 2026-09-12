# python_ldsc/reader.R: function bodies preserved from the original workflow.

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

# Python --rg reports the heritability of p2, not p1. Keep that mapping
# explicit, including when both orientations of a pair are supplied.

trait_heritability_scales <- function(rows, traits) {
  if ("H2_Scale" %in% names(rows)) {
    by_trait <- unique(rows[, .(p2, H2_Scale)])
    if (anyDuplicated(by_trait$p2)) stop("Conflicting scale metadata for a trait.", call. = FALSE)
    return(by_trait$H2_Scale[match(traits, by_trait$p2)])
  }
  rep(attr(rows, "heritability_scale", exact = TRUE), length(traits))
}

validate_covariance_heritability <- function(rows, requested_scale) {
  scales <- unique(rows$H2_Scale)
  if (requested_scale == "auto" && (length(scales) != 1L ||
      !scales %in% c("observed", "liability")))
    stop("Covariance auto requires a complete common scale; mixed scales require explicit choice.", call. = FALSE)
  self <- rows[p1 == p2]
  if (!nrow(self) || !all(unique(rows$p2) %in% self$p2) ||
      any(!is.finite(self$h2) | self$h2 <= 0 | !is.finite(self$h2_se) | self$h2_se <= 0))
    stop("Covariance PCA requires complete finite positive self heritability and SE for every selected trait; no fallback or conversion.", call. = FALSE)
  if (requested_scale == "mixed")
    warning("Explicit mixed-scale covariance: trait-specific scales affect PCA loadings; confirm phenotype/sample-size conventions. No conversion performed.", call. = FALSE)
  invisible(TRUE)
}

normalize_heritability_columns <- function(ldsc_rows,
    requested_scale = c("auto", "liability", "observed", "mixed"),
    pca_matrix = c("correlation", "covariance")) {
  requested_scale <- match.arg(requested_scale)
  pca_matrix <- match.arg(pca_matrix)
  normalized <- as.data.table(copy(ldsc_rows))
  inherited <- attr(ldsc_rows, "heritability_scale", exact = TRUE)
  internal <- internal_heritability_columns %in% names(normalized)
  if (any(internal) && !all(internal)) stop("Internal h2/h2_se columns are incomplete.", call. = FALSE)
  if (all(internal)) {
    if (is.null(inherited)) {
      if (!requested_scale %in% c("observed", "liability"))
        stop("Normalized h2/h2_se columns lack their heritability-scale metadata.", call. = FALSE)
      inherited <- requested_scale
    }
    if (requested_scale %in% c("observed", "liability") && inherited != requested_scale)
      stop("Requested scale conflicts with normalized input scale.", call. = FALSE)
    if (!"H2_Scale" %in% names(normalized)) {
      if (inherited == "mixed") stop("Mixed normalized input lacks per-trait scale metadata.", call. = FALSE)
      normalized[, H2_Scale := inherited]
    }
    used <- unique(normalized$H2_Scale[normalized$H2_Scale %in% c("observed", "liability")])
    attr(normalized, "heritability_scale") <- if (length(used) == 1L) used else "mixed"
    attr(normalized, "heritability_value_column") <- if (length(used) == 1L) heritability_column_sets[[used]][1L] else "trait-specific h2"
    attr(normalized, "heritability_se_column") <- if (length(used) == 1L) heritability_column_sets[[used]][2L] else "trait-specific h2_se"
    if (pca_matrix == "covariance") validate_covariance_heritability(normalized, requested_scale)
    return(normalized)
  }
  sources <- intersect(unlist(heritability_column_sets), names(normalized))
  if (!length(sources)) stop("Supply h2_obs/h2_obs_se or h2_liab/h2_liab_se.", call. = FALSE)
  populated <- function(x) !is.na(x) & !trimws(as.character(x)) %in% c("", "NA", "NaN", "nan")
  # An unused, entirely empty pair may be ignored, even if its other column is absent.
  for (scale in names(heritability_column_sets)) {
    cols <- heritability_column_sets[[scale]]
    present <- intersect(cols, names(normalized))
    if (length(present) == 1L && any(populated(normalized[[present]])))
      stop("Incomplete heritability column pair for ", scale, ".", call. = FALSE)
    for (col in setdiff(cols, names(normalized))) normalized[, (col) := NA_real_]
  }
  traits <- unique(as.character(normalized$p2))
  availability_table <- normalized[, lapply(heritability_column_sets, function(cols)
    any(populated(.SD[[cols[1]]]) | populated(.SD[[cols[2]]]))), by = p2]
  scale_names <- names(heritability_column_sets)
  availability <- lapply(seq_along(traits), function(i)
    setNames(as.logical(availability_table[i, ..scale_names]), scale_names))
  chosen <- requested_scale
  if (pca_matrix == "covariance" && requested_scale == "auto") {
    complete <- vapply(heritability_column_sets, function(cols) {
      self <- normalized[p1 == p2]
      nrow(self) > 0L && all(traits %in% self$p2) &&
        all(populated(self[[cols[1]]]) & populated(self[[cols[2]]]))
    }, logical(1))
    if (sum(complete) != 1L)
      stop("Covariance auto requires exactly one complete common scale. Choose observed/liability explicitly, or mixed for trait-specific scales; no conversion is performed.", call. = FALSE)
    chosen <- names(complete)[complete]
  }
  scales <- vapply(seq_along(traits), function(i) {
    if (chosen %in% c("observed", "liability")) return(chosen)
    available <- names(availability[[i]])[availability[[i]]]
    if (length(available) > 1L)
      stop("Ambiguous heritability scale for trait ", traits[i],
           ": both scales are populated; choose --heritability_scale observed or liability.", call. = FALSE)
    if (!length(available)) return("unavailable")
    available
  }, character(1))
  normalized[, H2_Scale := scales[match(p2, traits)]]
  normalized[, `:=`(h2 = NA_real_, h2_se = NA_real_)]
  for (scale in names(heritability_column_sets)) {
    cols <- heritability_column_sets[[scale]]
    index <- which(normalized$H2_Scale == scale)
    for (i in seq_along(cols))
      set(normalized, i = index, j = internal_heritability_columns[i],
          value = suppressWarnings(as.numeric(normalized[[cols[i]]][index])))
  }
  normalized[, (unlist(heritability_column_sets)) := NULL]
  used <- unique(scales[scales != "unavailable"])
  attr(normalized, "heritability_scale") <- if (length(used) == 1L) used else "mixed"
  attr(normalized, "heritability_value_column") <- if (length(used) == 1L) heritability_column_sets[[used]][1L] else "trait-specific h2"
  attr(normalized, "heritability_se_column") <- if (length(used) == 1L) heritability_column_sets[[used]][2L] else "trait-specific h2_se"
  if (pca_matrix == "covariance") validate_covariance_heritability(normalized, requested_scale)
  normalized
}

# Stream the potentially multi-million-row table and retain only selected
# trait-by-trait rows. This avoids holding a larger all-trait table in memory.

read_python_ldsc_selected <- function(path, trait_order, chunk_size = 250000L,
                                      heritability_scale = c(
                                        "auto", "liability", "observed", "mixed"
                                      ), pca_matrix = c("correlation", "covariance")) {
  heritability_scale <- match.arg(heritability_scale)
  pca_matrix <- match.arg(pca_matrix)
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
  source_columns <- c(required_python_ldsc_base_columns,
                      intersect(unlist(heritability_column_sets), names(header_table)))

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
      replicate(length(source_columns), logical(0), simplify = FALSE),
      source_columns
    ))
  } else {
    selected <- rbindlist(retained_chunks, use.names = TRUE)
  }

  selected <- normalize_heritability_columns(selected, heritability_scale, pca_matrix)
  attr(selected, "source_rows_read") <- rows_read
  attr(selected, "selected_traits_seen") <- trait_order[selected_traits_seen]
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
