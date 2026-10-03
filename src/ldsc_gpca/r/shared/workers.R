# Shared bounded retries and defensive worker-result reporting.
normalize_worker_result <- function(result, chromosome) {
  valid <- is.list(result) && is.logical(result$success) &&
    length(result$success) == 1L && !is.na(result$success) &&
    is.character(result$output) && length(result$output) == 1L &&
    is.character(result$error) && length(result$error) == 1L
  if (!valid) {
    detail <- if (inherits(result, "try-error")) as.character(result)[1L] else
      "Parallel worker did not return a valid result (it may have exited or been killed)"
    result <- list(success = FALSE, output = NA_character_, error = detail)
  }
  if (isTRUE(result$success) && (is.na(result$output) || !nzchar(result$output))) {
    result$success <- FALSE
    result$error <- "Worker returned success without an output name"
  }
  if (!result$success && (is.na(result$error) || !nzchar(result$error)))
    result$error <- "Worker failed without an error message"
  if (!is.null(result$attempts) && (!is.numeric(result$attempts) || length(result$attempts) != 1L ||
      !is.finite(result$attempts) || result$attempts < 0)) result$attempts <- 1L
  result$chromosome <- chromosome
  result
}

run_parallel_workers <- function(chromosomes, worker, n_cores, outdir,
                                 stage = "GWAMA", label_column = "Chromosome") {
  if (n_cores == 0L) {
    available <- parallel::detectCores(logical = FALSE)
    n_cores <- if (is.na(available)) 1L else max(1L, available - 1L)
  }
  results <- lapply(chromosomes, function(chr) {
    result <- normalize_worker_result(NULL, chr)
    result$attempts <- 0L
    result$error <- "Job has not completed"
    result
  })
  status_file <- file.path(outdir, paste0(stage, "_Run_Status.csv"))
  attempt_file <- file.path(outdir, paste0(stage, "_Worker_Attempts.csv"))
  write.csv(make_run_status(results, label_column), status_file, row.names = FALSE)
  history <- list()
  for (attempt in 1:2) {
    pending <- which(!vapply(results, function(x) x$success, logical(1)))
    if (!length(pending)) break
    cores <- if (attempt == 1L) min(n_cores, length(pending)) else 1L
    message(stage, " attempt ", attempt, "/2: ", length(pending), " job(s), ", cores, " worker(s)")
    call_worker <- function(i) tryCatch(worker(chromosomes[[i]]),
      error = function(e) structure(conditionMessage(e), class = "try-error"))
    run <- function() {
      if (.Platform$OS.type == "windows") return(lapply(pending, call_worker))
      if (cores > 1L) return(parallel::mclapply(pending, call_worker,
        mc.cores = cores, mc.preschedule = FALSE))
      # mclapply(mc.cores=1) runs in the parent; keep even serial retries isolated.
      lapply(pending, function(i) {
        child <- parallel::mcparallel(call_worker(i))
        value <- parallel::mccollect(child)
        if (is.null(value)) NULL else value[[1L]]
      })
    }
    returned <- tryCatch(run(), error = function(e)
      rep(list(structure(conditionMessage(e), class = "try-error")), length(pending)))
    for (j in seq_along(pending)) {
      i <- pending[j]
      result <- normalize_worker_result(returned[[j]], chromosomes[[i]])
      result$attempts <- attempt
      results[[i]] <- result
      history[[length(history) + 1L]] <- make_run_status(list(result), label_column)
    }
    write.csv(do.call(rbind, history), attempt_file, row.names = FALSE)
    write.csv(make_run_status(results, label_column), status_file, row.names = FALSE)
  }
  results
}

make_run_status <- function(gwama_results, label_column = "Chromosome") {
  gwama_results <- lapply(seq_along(gwama_results), function(i) {
    result <- gwama_results[[i]]
    chromosome <- if (is.list(result) && length(result$chromosome) <= 1L) result$chromosome else i
    normalize_worker_result(result, chromosome)
  })
  status <- data.frame(
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
    Attempts = vapply(gwama_results, function(result)
      if (is.null(result$attempts)) 1L else as.integer(result$attempts), integer(1)),
    stringsAsFactors = FALSE
  )
  names(status)[1L] <- label_column
  status
}
