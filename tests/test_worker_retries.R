root <- commandArgs(TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
out <- tempfile(); dir.create(out)
controller_pid <- Sys.getpid()
ok <- function(chr) list(success=TRUE, chromosome=chr, output=paste0("Chr", chr), error=NA_character_)
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error=identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}

# Defensive reporting cannot itself crash on bad child results.
bad <- list(NULL, structure("upstream failure", class="try-error"), list(success=FALSE),
            list(success=TRUE, output=character(), error=NA_character_),
            list(success=c(TRUE, FALSE), output="x", error=NA_character_))
for(value in bad) stopifnot(!make_run_status(list(ok(1), value))$Success[2])

for(kind in c("NULL", "try-error", "malformed", "error", "exit")) {
  if(kind == "exit" && .Platform$OS.type == "windows") next
  folder <- file.path(out, kind); dir.create(folder)
  worker <- function(chr) {
    path <- file.path(folder, paste0("calls",chr))
    calls <- if(file.exists(path)) as.integer(readLines(path)) + 1L else 1L
    writeLines(as.character(calls), path)
    if(chr == 7L && calls == 1L) switch(kind,
      "NULL" = return(NULL), "try-error" = return(structure("temporary failure", class="try-error")),
      "malformed" = return(list(success=FALSE)), "error" = stop("temporary failure"),
      "exit" = {
        stopifnot(Sys.getpid() != controller_pid)
        tools::pskill(Sys.getpid(), signal=9L)
      })
    ok(chr)
  }
  status <- make_run_status(suppressWarnings(run_parallel_workers(as.list(c(7L, 2L, 9L)), worker, 2L, folder)))
  stopifnot(all(status$Success), identical(status$Chromosome, c("7", "2", "9")),
            identical(status$Attempts, c(2L, 1L, 1L)), readLines(file.path(folder,"calls2")) == "1")
  attempts <- read.csv(file.path(folder, "GWAMA_Worker_Attempts.csv"))
  stopifnot(nrow(attempts) == 4L, sum(!attempts$Success) == 1L)
}
folder <- file.path(out, "permanent"); dir.create(folder)
status <- make_run_status(run_parallel_workers(list(NULL), function(chr) stop("persistent failure"), 1L, folder))
stopifnot(identical(status$Chromosome, "whole_genome"), !status$Success, status$Attempts == 2L,
          grepl("persistent failure", status$Error))

# Both analysis backends must stop after retry exhaustion, retaining the status CSV.
traits <- c("A", "B", "C")
manifest <- file.path(out, "traits.csv"); write.csv(data.frame(traitname=traits), manifest, row.names=FALSE)
rows <- data.frame(p1=c("A","A","A","B","B","C"), p2=c("A","B","C","B","C","C"),
  rg=c(1,.2,.3,1,.4,1), se=.1, h2_obs=.2, h2_obs_se=.03, h2_int=1.05, h2_int_se=.01,
  gcov_int=c(1.05,.1,.2,1.05,.3,1.05), gcov_int_se=.01)
rows$z <- rows$rg/rows$se; rows$p <- 2*pnorm(-abs(rows$z))
ldsc <- file.path(out,"ldsc.csv"); write.csv(rows, ldsc, row.names=FALSE)
inputs <- file.path(out,"inputs"); dir.create(inputs)
for(chr in 1:22) for(trait in traits) write.table(data.frame(SNPID=paste0("rs",chr),CHR=chr,BP=100,
  EA="A",OA="G",EAF=.2,N=50000,Z=1,P=.3), file.path(inputs,paste0(trait,"_chr",chr,"_GenomicPCA_inputs.tsv")),
  sep="\t",row.names=FALSE,quote=FALSE)
bad_source <- file.path(out,"bad_gwama.R")
writeLines('my_GWAMA <- function(...) stop("simulated GWAMA execution failure")',bad_source)
for(backend in c("python_ldsc", "genomicsem")) {
  folder <- file.path(out, backend)
  if(backend == "python_ldsc") {
    args <- list(input=manifest, ldsc_results=ldsc, outdir=folder, gpca_input_folder=inputs,
      source_path=bad_source, splitby_chr="split", n_cores=2L, allow_missing_traits=FALSE,
      failed_ldsc_action="error", pca_matrix="correlation", heritability_scale="auto",
      rg_normalization="pair", pc1_orientation="tutorial", validate_only=FALSE, ldsc_chunk_size=250000L,
      duplicate_tolerance=1e-6, duplicate_z_tolerance=1e-4, self_rg_tolerance=1e-6,
      comparison_epsilon=1e-12, h2_z_warn_threshold=2, z_consistency_tolerance=.05,
      z_consistency_action="warn", rg_out_of_range_action="warn", negative_eigen_action="warn",
      matrix_eigen_tolerance=1e-8)
    dir.create(folder)
    fail(gpsca_analysis(args), "after 2 attempts")
  } else {
    source(file.path(root,"src/ldsc_gpca/r/gpsca_gwama_v2.r"))
    # Inject only the validated LDSC read/QC boundary; exercise the real GWAMA workflow.
    S <- matrix(.04,3,3,dimnames=list(traits,traits)); diag(S) <- .2
    I <- diag(1.05,3); dimnames(I) <- list(traits,traits)
    LDSCoutput <- list(S=S, S_Stand=cov2cor(S), I=I)
    rdata <- file.path(out,"ldsc.RData"); save(LDSCoutput,file=rdata)
    qc_genomicsem <- function(x, ...) x
    write_genomicsem_diagnostics <- function(...) NULL
    fail(genomicsem_main(c("--input",manifest,"--ldsc_results",rdata,"--outdir",folder,
      "--gpca_input_folder",inputs,"--source_path",bad_source,"--n_cores","2")), "after 2 attempts")
  }
  status <- read.csv(file.path(folder,"GWAMA_Run_Status.csv"))
  stopifnot(nrow(status)==22L, all(!status$Success), all(status$Attempts==2L),
            all(grepl("simulated GWAMA",status$Error)))
}

# Native GenomicSEM munging uses the same retry runner; exhausted failures never reach LDSC.
source(file.path(root,"src/ldsc_gpca/r/genomicsem/ldsc_pipeline.R"))
folder <- file.path(out,"native_munge"); dir.create(folder)
resolved <- file.path(out,"resolved.csv")
write.csv(data.frame(traitname=traits,source_file="fixture",N=50000,
  sampleprevalence=NA,populationprevalence=NA),resolved,row.names=FALSE)
fail(run_genomicsem_ldsc(resolved,folder,"ref","ref","hm3","munge",2L,.9,.01,22L,200L,NA,"error",
  munge_fun=function(...) stop("simulated native munge failure"),
  ldsc_fun=function(...) stop("LDSC MUST NOT START")), "munging failed after 2 attempts")
status <- read.csv(file.path(folder,"GenomicSEM_Munging_Run_Status.csv"))
stopifnot(identical(status$Trait,traits),all(!status$Success),all(status$Attempts==2L))

# A recovered trait reaches LDSC in manifest order; successful traits run only once.
folder <- file.path(out,"native_munge_retry"); dir.create(folder)
munge_once <- function(files, hm3, trait.names, N, info.filter, maf.filter, parallel, cores, overwrite) {
  stopifnot(!parallel, cores==1L, overwrite, N==50000, info.filter==.9, maf.filter==.01)
  marker <- paste0(trait.names,".calls")
  n <- if(file.exists(marker)) as.integer(readLines(marker))+1L else 1L
  writeLines(as.character(n),marker)
  if(trait.names=="B" && n==1L) stop("temporary native munge failure")
  con <- gzfile(paste0(trait.names,".sumstats.gz"),"wt")
  writeLines(c("SNP\tA1\tA2\tN\tZ","rs1\tA\tG\t50000\t1"),con); close(con)
}
fail(run_genomicsem_ldsc(resolved,folder,"ref","ref","hm3","munge",2L,.9,.01,22L,200L,NA,"error",
  munge_fun=munge_once, ldsc_fun=function(traits, trait.names, ...) {
    stopifnot(identical(trait.names,c("A","B","C")),
              identical(basename(traits),paste0(trait.names,".sumstats.gz")),all(file.exists(traits)))
    stop("REACHED LDSC AFTER SUCCESSFUL MUNGING")
  }), "REACHED LDSC AFTER SUCCESSFUL MUNGING")
status <- read.csv(file.path(folder,"GenomicSEM_Munging_Run_Status.csv"))
stopifnot(identical(status$Trait,traits),all(status$Success),identical(status$Attempts,c(1L,2L,1L)))
for(i in seq_along(traits)) stopifnot(as.integer(readLines(file.path(folder,"munge_output",
  paste0(traits[i],".calls"))))==status$Attempts[i])
unlink(out,recursive=TRUE)
cat("Worker retry, actual child exit, both GWAMA backends and native munging orchestration tests passed\n")
