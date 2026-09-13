nextflow.enable.dsl=2

process CHECK_IMAGE {
    publishDir params.outdir, mode: 'copy'
    output:
    path 'tools.txt'
    script:
    '''
    {
      ldsc-gpca --version
      ldsc-gpca gpca --help
      ldsc-gpca genomicsem gpca --help
      bcftools --version
      python -c 'import os; from ldsc_gpca.ldsc_runtime import check_runtime; check_runtime(conda=os.environ["CONDA_EXE"], prefix=os.environ["LDSC_GPCA_LDSC_PREFIX"]); print("LDSC runtime OK")'
      Rscript -e 'library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "\\n")'
    } > tools.txt
    '''
}

workflow {
    CHECK_IMAGE()
}
