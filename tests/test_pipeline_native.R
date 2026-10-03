# Native two-pass orchestration with deterministic injected LDSC estimates.
# This tests cutoff forwarding and retained-trait handoff, not real regressions.
root <- commandArgs(TRUE)[1]
source(file.path(root, "src/ldsc_gpca/r/genomicsem/reader.R"))
source(file.path(root, "src/ldsc_gpca/r/genomicsem/ldsc_pipeline.R"))
temporary <- tempfile("pipeline-native-")
dir.create(temporary)
names <- c("002", "001", "003")
paths <- file.path(temporary, paste0(names, ".sumstats"))
for (path in paths) writeLines(c("SNP\tA1\tA2\tN\tZ", "rs1\tG\tA\t1000\t2"), path)
manifest <- file.path(temporary, "manifest.csv")
write.csv(data.frame(traitname=names, source_file=paths, sampleprevalence="NA",
                    populationprevalence="NA", N="NA"), manifest, row.names=FALSE)
out <- file.path(temporary, "out")
dir.create(out)
calls <- list()
fake_ldsc <- function(traits, sample.prev, population.prev, trait.names, ld, wld,
                      sep_weights, chr, n.blocks, chisq.max, stand, ldsc.log) {
  calls[[length(calls)+1L]] <<- list(traits=trait.names, cutoff=chisq.max, stand=stand)
  n <- length(trait.names)
  S <- diag(if (stand) c(.2,.3) else c(.2,.3,-.1))
  dimnames(S) <- list(trait.names,trait.names)
  I <- diag(n);dimnames(I) <- dimnames(S)
  m <- n*(n+1L)/2L
  list(S=S, I=I, S_Stand=I, V=diag(.001,m), V_Stand=diag(.001,m))
}
run_genomicsem_ldsc(manifest, out, "ld", "weights", "", "existing", 2,
                   .9, .01, 22, 200, 80, "drop", ldsc_fun=fake_ldsc)
stopifnot(length(calls)==2L, identical(calls[[1]]$stand,FALSE),
          identical(calls[[2]]$stand,TRUE), calls[[1]]$cutoff==80,
          calls[[2]]$cutoff==80, identical(calls[[2]]$traits,c("002","001")))
selected <- read.csv(file.path(out,"Selected_Traits.csv"),colClasses="character")
stopifnot(identical(selected$traitname,c("002","001")))
env <- new.env()
load(file.path(out,"genomicPCA_LDSC.RData"),envir=env)
stopifnot(identical(colnames(env$LDSCoutput$S),c("002","001")),
          all(c("S","V","I","S_Stand","V_Stand") %in% names(env$LDSCoutput)))
unlink(temporary,recursive=TRUE)
cat("Native two-pass cutoff=80 and retained-trait order checks passed\n")
