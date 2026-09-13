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
      --python_ldsc '${ldsc}' \\
      --outdir qc --validate_only
    """
}

workflow {
    if (!params.manifest || !params.ldsc)
        error 'Supply --manifest MANIFEST.csv and --ldsc PAIRWISE.csv[.gz] (local paths or gs:// URLs).'
    GPCA_QC(file(params.manifest, checkIfExists: true), file(params.ldsc, checkIfExists: true))
}
