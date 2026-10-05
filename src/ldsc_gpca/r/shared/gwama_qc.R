# Validate exported associations outside the unchanged GWAMA implementation.
validate_gwama_output <- function(path, dat, loadings, audit_prefix) {
  issues_path <- paste0(audit_prefix, ".GWAMA_QC_Issues.csv")
  summary_path <- paste0(audit_prefix, ".GWAMA_QC_Summary.csv")
  issues <- data.table(Row = integer(), SNPID = character(), Reason = character())
  publish <- function(total) {
    fwrite(issues, issues_path, na = "NA")
    fwrite(data.table(Status = if (nrow(issues)) "failed" else "passed",
      Total_Rows = total, Invalid_Rows = uniqueN(issues$Row),
      Issue_Count = nrow(issues), Policy = "error", Raw_Output = path), summary_path)
    if (nrow(issues)) stop("GWAMA output QC failed: ", issues$Reason[1L],
      "; raw output retained. See ", issues_path, call. = FALSE)
    invisible(TRUE)
  }
  structural_failure <- function(reason, total = 0L) {
    issues <<- data.table(Row = NA_integer_, SNPID = "", Reason = reason)
    publish(total)
  }
  header <- tryCatch(fread(path, nrows = 0L, check.names = FALSE,
                          showProgress = FALSE), error = identity)
  if (inherits(header, "error")) structural_failure(conditionMessage(header))
  required <- c("SNPID", "Direction", "BETA", "SE", "Z", "PVAL", "N_eff")
  missing <- setdiff(required, names(header))
  if (anyDuplicated(names(header)) || length(missing))
    structural_failure(paste("Missing/duplicate GWAMA columns:", paste(missing, collapse = ", ")))
  columns <- c(required, intersect("N_obs", names(header)))
  result <- tryCatch(fread(path, select = columns, colClasses = "character",
    na.strings = NULL, check.names = FALSE, showProgress = FALSE), error = identity)
  if (inherits(result, "error")) structural_failure(conditionMessage(result))
  if (!nrow(result)) structural_failure("Empty GWAMA result")
  add_issue <- function(invalid, reason) {
    rows <- which(invalid)
    if (length(rows)) issues <<- rbindlist(list(issues,
      data.table(Row = rows, SNPID = result$SNPID[rows], Reason = reason)))
  }
  add_issue(is.na(result$SNPID) | trimws(result$SNPID) == "" |
    duplicated(result$SNPID) | duplicated(result$SNPID, fromLast = TRUE),
    "missing_or_duplicate_SNPID")
  numeric_columns <- setdiff(columns, c("SNPID", "Direction"))
  for (column in numeric_columns) {
    values <- suppressWarnings(as.numeric(result[[column]]))
    valid <- is.finite(values)
    if (column %in% c("SE", "N_eff", "N_obs")) valid <- valid & values > 0
    if (column == "PVAL") valid <- valid & values >= 0 & values <= 1
    add_issue(!valid, if (column == "PVAL") "PVAL_not_finite_or_outside_0_1" else
      paste0(column, if (column %in% c("SE", "N_eff", "N_obs"))
        "_not_finite_positive" else "_not_finite"))
  }
  direction_valid <- !is.na(result$Direction) &
    grepl("^[+?-]+$", result$Direction) & nchar(result$Direction) == length(loadings)
  add_issue(!direction_valid, "invalid_Direction_or_trait_count")
  # Direction records actual post-alignment availability. Check each loading
  # separately: signed weights that cancel in sum may still define a valid Z.
  has_weight <- rep(FALSE, nrow(result))
  for (i in seq_along(loadings)) {
    if (loadings[i] == 0) next
    available <- substr(result$Direction, i, i) %in% c("+", "-")
    n <- suppressWarnings(as.numeric(dat[[i]]$N[match(result$SNPID, dat[[i]]$SNPID)]))
    has_weight <- has_weight | (available & is.finite(n) & n > 0)
  }
  add_issue(direction_valid & !has_weight, "zero_available_weight")
  publish(nrow(result))
}
