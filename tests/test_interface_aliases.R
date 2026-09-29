root <- commandArgs(TRUE)[1]
library(argparse); library(data.table)
source(file.path(root,'src/ldsc_gpca/r/shared/input.R'))
bundled_gwama_path <- 'unchanged'
source(file.path(root,'src/ldsc_gpca/r/genomicsem/cli.R'))
source(file.path(root,'src/ldsc_gpca/r/python_ldsc/cli.R'))
fail <- function(expr,pattern) {
  msg <- tryCatch({force(expr);''},error=conditionMessage)
  stopifnot(grepl(pattern,msg))
}
base <- c('--input','traits.csv','--ldsc_results','results','--outdir','out','--n_cores','22')
p <- genomicsem_parser(); x <- p$parse_args(base)
stopifnot(x$n_cores==22L,x$ldsc_results=='results')
for (flag in c('--ldsc_path','--input_file','--output_folder','--cores','--ld'))
  fail(p$parse_args(c(base,flag,'x')),'unrecognized arguments')
arguments <- base
commandArgs <- function(trailingOnly=FALSE) arguments
x <- parse_command_line();stopifnot(x$n_cores==22L,x$ldsc_results=='results')
for (flag in c('--python_ldsc','--ldsc_path','--cores','--tolerance','--ld')) {
  arguments <- c(base,flag,'x');fail(parse_command_line(),'unrecognized arguments')
}
tmp <- tempfile(fileext='.csv')
writeLines(c('traitname,vcf_files','002,a','001,b'),tmp)
stopifnot(identical(read_trait_manifest(tmp),c('002','001')))
for (old in c('gwas_name','sampleprevalence','populationprevalence','pop_prevalence','munge_inputs','traits')) {
  writeLines(c(paste0('traitname,',old),'A,A','B,B'),tmp)
  fail(read_trait_manifest(tmp),'Unsupported manifest headers')
}
writeLines(c('gwas_name','A','B'),tmp);fail(read_trait_manifest(tmp),'Unsupported manifest headers')
writeLines(c('traitname,traitname','A,A','B,B'),tmp);fail(read_trait_manifest(tmp),'unique')
writeLines(c('traitname','A','A'),tmp);fail(read_trait_manifest(tmp),'unique')
unlink(tmp);cat('R single-name checks passed\n')
