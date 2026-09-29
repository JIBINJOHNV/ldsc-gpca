nextflow.enable.dsl=2

process GPCA_QC {
    publishDir params.outdir, mode: 'copy'
    input:
    path manifest, stageAs: 'manifest.csv'
    path ldsc, stageAs: 'ldsc/*'
    output:
    path 'qc'
    script:
    """
    ldsc-gpca gpca \\
      --input '${manifest}' \\
      --ldsc_results '${ldsc}' \\
      --outdir qc --validate_only
    """
}

workflow {
    if (!params.input || !params.ldsc_results)
        error 'Supply --input MANIFEST.csv and --ldsc_results PAIRWISE.csv[.gz] (local paths or gs:// URLs).'
    GPCA_QC(file(params.input, checkIfExists: true), file(params.ldsc_results, checkIfExists: true))
}
