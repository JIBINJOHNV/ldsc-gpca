# GenomicSEM two-pass LDSC; injected functions permit deterministic orchestration tests.
run_genomicsem_ldsc <- function(manifest, outdir, ld, wld, hm3, mode, cores,
                               info_filter, maf_filter, chromosomes, n_blocks,
                               chisq_max, invalid_h2_action,
                               munge_fun = GenomicSEM::munge, ldsc_fun = GenomicSEM::ldsc) {
  dat <- read.csv(manifest, colClasses = "character", check.names = FALSE)
  for (column in c("sampleprevalence", "populationprevalence", "N"))
    dat[[column]] <- as.numeric(dat[[column]])
  old <- setwd(outdir)
  on.exit(setwd(old), add = TRUE)
  events <- data.frame(stage = character(), action = character(), reason = character())
  stage <- "initialization"
  record <- function(action, reason) {
    events[nrow(events) + 1L, ] <<- list(stage, action, reason)
    write.csv(events, file.path(outdir, "GenomicSEM_LDSC_Events.csv"), row.names = FALSE)
  }
  write.csv(events, "GenomicSEM_LDSC_Events.csv", row.names = FALSE)
  on.exit(writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt")), add = TRUE)
  validate_munged <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    lines <- tryCatch(readLines(con, n = 2L, warn = FALSE), finally = close(con))
    if (length(lines) < 2L || !all(c("SNP", "A1", "A2", "N", "Z") %in% strsplit(lines[1], "\t", fixed = TRUE)[[1]]))
      stop(paste("Munged file needs tab-separated SNP,A1,A2,N,Z and data:", path))
  }
  tryCatch(withCallingHandlers({
    if (mode == "munge") {
      stage <- "munge"
      dir.create("munge_output")
      worker <- function(name) {
        before <- setwd(file.path(outdir, "munge_output"))
        on.exit(setwd(before))
        row <- dat[match(name, dat$traitname), , drop = FALSE]
        paths <- file.path(outdir, "munge_output", paste0(name, c(".sumstats.gz", ".sumstats")))
        previous <- file.info(paths)[, c("size", "mtime", "ctime")]
        # The package owns scheduling/retries; the upstream single-trait algorithm is unchanged.
        munge_fun(files = row$source_file, hm3 = hm3, trait.names = name,
          N = row$N, info.filter = info_filter, maf.filter = maf_filter,
          parallel = FALSE, cores = 1L, overwrite = TRUE)
        candidates <- paths[file.exists(paths)]
        if (length(candidates) != 1L) stop(paste("Missing or ambiguous munged output:", name))
        if (identical(file.info(candidates)[, c("size", "mtime", "ctime")], previous[candidates, , drop = FALSE]))
          stop(paste("No fresh munged output:", name))
        validate_munged(candidates)
        list(success = TRUE, chromosome = name, output = candidates, error = NA_character_)
      }
      runs <- run_parallel_workers(as.list(dat$traitname), worker, cores, outdir,
                                   stage = "GenomicSEM_Munging", label_column = "Trait")
      status <- make_run_status(runs, "Trait")
      if (any(!status$Success)) stop("Analysis stopped: GenomicSEM munging failed after 2 attempts. ",
        "See GenomicSEM_Munging_Worker_Attempts.csv.\n",
        paste(paste0(status$Trait[!status$Success], ": ", status$Error[!status$Success]), collapse = "\n"))
      dat$traits <- status$Output
    } else dat$traits <- dat$source_file
    stage <- "munged_input_validation"
    for (path in dat$traits) validate_munged(path)
    estimate <- function(selected, stand, log_name) ldsc_fun(
      traits = selected$traits, sample.prev = selected$sampleprevalence,
      population.prev = selected$populationprevalence, trait.names = selected$traitname,
      ld = ld, wld = wld, sep_weights = !identical(ld, wld), chr = chromosomes,
      n.blocks = n_blocks, chisq.max = chisq_max, stand = stand, ldsc.log = log_name)
    stage <- "raw_ldsc"
    LDSCoutput_raw <- estimate(dat, FALSE, "GenomicSEM_raw")
    save(LDSCoutput_raw, file = "genomicsem_LDSC_raw.RData")
    if (!is.matrix(LDSCoutput_raw$S) || !identical(dim(LDSCoutput_raw$S), rep(nrow(dat), 2L)) ||
        !identical(colnames(LDSCoutput_raw$S), dat$traitname))
      stop("Raw LDSC S dimensions/trait order do not match the manifest.")
    h2 <- diag(LDSCoutput_raw$S)
    bad <- !is.finite(h2) | h2 <= 0
    qc <- data.frame(traitname = dat$traitname, h2_estimate = as.numeric(h2),
      scale = ifelse(is.na(dat$populationprevalence), "observed", "liability"),
      QC_action = ifelse(bad, if (invalid_h2_action == "drop") "removed" else "error", "retained"),
      QC_reason = ifelse(!is.finite(h2), "Non-finite h2", ifelse(h2 <= 0, "Non-positive h2", "")))
    write.csv(qc, "GenomicSEM_LDSC_Trait_QC.csv", row.names = FALSE)
    if (any(bad) && invalid_h2_action == "error") stop("Invalid heritability; see trait QC audit.")
    dat <- dat[!bad, , drop = FALSE]
    write.csv(dat, "Selected_Traits.csv", row.names = FALSE)
    if (nrow(dat) < 2L) stop("Fewer than two traits remain after heritability QC.")
    stage <- "standardized_ldsc"
    LDSCoutput <- estimate(dat, TRUE, "GenomicSEM_standardized")
    LDSCoutput <- validate_genomicsem_structure(LDSCoutput)
    if (!identical(colnames(LDSCoutput$S), dat$traitname) ||
        any(diag(LDSCoutput$S) <= 0) ||
        any(vapply(LDSCoutput[c("S", "I", "V", "S_Stand", "V_Stand")], function(x) any(!is.finite(x)), logical(1))))
      stop("Final LDSC has invalid estimates or trait ordering; final RData not published.")
    save(LDSCoutput, file = "genomicsem_LDSC.RData")
    record("completed", paste(nrow(dat), "traits; GenomicSEM standardized output saved."))
    message("Saved ", file.path(outdir, "genomicsem_LDSC.RData"))
    invisible(LDSCoutput)
  }, warning = function(w) record("warning", conditionMessage(w))), error = function(e) {
    record("error", conditionMessage(e))
    stop(e)
  })
}
