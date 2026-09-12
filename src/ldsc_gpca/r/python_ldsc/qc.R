# python_ldsc/qc.R: function bodies preserved from the original workflow.

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
                                    "auto", "liability", "observed", "mixed"
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

  ldsc_rows <- as.data.table(copy(ldsc_rows))
  ldsc_rows <- ldsc_rows[trimws(as.character(p1)) %chin% trait_order & trimws(as.character(p2)) %chin% trait_order]
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
  trait_scales <- trait_heritability_scales(ldsc_rows, trait_order)
  result[, `:=`(
    Heritability_Scale = trait_scales,
    Heritability_Source_Column = ifelse(trait_scales == "observed", "h2_obs", ifelse(trait_scales == "liability", "h2_liab", NA_character_)),
    Heritability_SE_Source_Column = ifelse(trait_scales == "observed", "h2_obs_se", ifelse(trait_scales == "liability", "h2_liab_se", NA_character_)),
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
                                             "auto", "liability", "observed", "mixed"
                                           )) {
  action <- match.arg(action)
  heritability_scale <- match.arg(heritability_scale)
  source_rows_read <- attr(ldsc_rows, "source_rows_read", exact = TRUE)
  ldsc_rows <- as.data.table(copy(ldsc_rows))
  ldsc_rows <- ldsc_rows[trimws(as.character(p1)) %chin% trait_order & trimws(as.character(p2)) %chin% trait_order]
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
    selected_scale <- self_pair_qc$Heritability_Scale[trait_index]
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
                                             "auto", "liability", "observed", "mixed"
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

  ldsc_rows <- as.data.table(copy(ldsc_rows))
  ldsc_rows <- ldsc_rows[trimws(as.character(p1)) %chin% trait_order & trimws(as.character(p2)) %chin% trait_order]
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

  trait_scales <- self_pair_qc$Heritability_Scale
  heritability_results <- data.frame(
    Trait = trait_order,
    h2_liab = ifelse(trait_scales == "liability", self_summary$h2, NA_real_),
    h2_liab_se = ifelse(trait_scales == "liability", self_summary$h2_se, NA_real_),
    Z = self_summary$h2 / self_summary$h2_se,
    P = 2 * pnorm(
      abs(self_summary$h2 / self_summary$h2_se),
      lower.tail = FALSE
    ),
    h2_int = self_summary$h2_int,
    h2_int_se = self_summary$h2_int_se,
    Heritability_Scale = trait_scales,
    Heritability = self_summary$h2,
    Heritability_SE = self_summary$h2_se,
    h2_obs = ifelse(trait_scales == "observed", self_summary$h2, NA_real_),
    h2_obs_se = ifelse(trait_scales == "observed", self_summary$h2_se, NA_real_),
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
