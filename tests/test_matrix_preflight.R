root <- commandArgs(TRUE)[1L]
source(file.path(root, "src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r"))
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
local({
  folder <- tempfile(); dir.create(folder); on.exit(unlink(folder, recursive = TRUE))
  traits <- c("C", "A", "B")
  named <- function(m) {dimnames(m) <- list(traits, traits); m}
  G <- named(matrix(.2,3,3)); diag(G) <- 1
  I <- named(matrix(c(1,.9,.9,.9,1,-.9,.9,-.9,1),3,3))
  args <- list(outdir=folder, negative_eigen_action="warn", matrix_eigen_tolerance=1e-8,
    pc1_orientation="tutorial", pca_matrix="correlation", validate_only=FALSE,
    pc1_only=FALSE, cti_action="error", max_cti_drop_fraction=0)
  fail(validate_analysis_matrices(G,I,traits,args), "not positive definite")
  audit <- read.csv(file.path(folder,"GenomicPCA_Matrix_Validation.csv"))
  stopifnot(audit$PC1_Status=="PASS",audit$CTI_Status=="INDEFINITE",!audit$GWAMA_Matrix_Eligible)
  args$pc1_only <- TRUE
  result <- suppressWarnings(validate_analysis_matrices(G,I,traits,args))
  stopifnot(identical(result$traits,traits),nrow(result$excluded)==0,
    identical(result$pc1$loadings,compute_pc1(G,traits)$loadings))
  args$pc1_only <- FALSE; args$cti_action <- "explore_drop"; args$max_cti_drop_fraction <- .34
  result <- suppressWarnings(validate_analysis_matrices(G,I,traits,args))
  stopifnot(length(result$traits)==2,nrow(result$excluded)==1,result$cti$pass,
    identical(result$CTI,I[result$traits,result$traits]))
  args$max_cti_drop_fraction <- .1
  fail(validate_analysis_matrices(G,I,traits,args),"exclusion limit")
  args$cti_action <- "error"; args$max_cti_drop_fraction <- 0
  fail(validate_analysis_matrices(-diag(3),I,traits,args),"order")
  fail(validate_analysis_matrices(named(-diag(3)),I,traits,args),"positive first eigenvalue")
  args$negative_eigen_action <- "error"
  fail(suppressWarnings(validate_analysis_matrices(I,named(diag(3)),traits,args)),"negative eigenvalues")
  audit <- read.csv(file.path(folder,"GenomicPCA_Matrix_Validation.csv"))
  stopifnot(audit$PC1_Status=="PASS",!audit$PCA_Matrix_Policy_Pass)
  args$negative_eigen_action <- "warn"
  good <- validate_analysis_matrices(G,named(diag(3)),traits,args)
  stopifnot(identical(good$pc1$loadings,compute_pc1(G,traits)$loadings))
  for (f in c(-1,1,Inf,NaN)) {a<-args;a$max_cti_drop_fraction<-f;fail(validate_matrix_options(a),'fraction')}
  a<-args;a$cti_action<-"explore_drop";fail(validate_matrix_options(a),'positive')
  a<-args;a$pc1_only<-TRUE;a$validate_only<-TRUE;fail(validate_matrix_options(a),'Choose')
  for (sep in c(",", "\t", " ")) for (compressed in c(FALSE, TRUE)) {
    lines <- c(paste(c("p1","p2","extra"),collapse=sep),
      paste(c("C","A","0.01000"),collapse=sep), paste(c("B","A","1e-99"),collapse=sep))
    path <- file.path(folder, if (compressed) "source.gz" else "source.txt")
    con <- if (compressed) gzfile(path,"wt") else file(path,"wt")
    writeLines(lines,con);close(con)
    target <- export_retained_ldsc(path,c("C","A"),folder,chunk_size=1L)
    stopifnot(identical(readLines(target),lines[1:2]))
  }
  for (m in list(named(matrix(1,3,3)),named(diag(c(1,1,1e-12))))) stopifnot(!assess_cti(m)$pass)
  bad<-G;bad[1,2]<-NA;fail(validate_analysis_matrices(bad,I,traits,args),'non-finite')
  bad<-I;bad[1,2]<-0;fail(validate_analysis_matrices(G,bad,traits,args),'asymmetric')
})
cat('Separate PC1/CTI checks, selection, boundaries, invalid inputs, and policy reporting passed\n')
