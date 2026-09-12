# python_ldsc/reporting.R: function bodies preserved from the original workflow.

write_analysis_outputs <- function(outdir, ldsc, pc1, correlation_matrix,
                                   pca_matrix, CTI,
                                   missing_traits, failed_traits,
                                   self_pair_qc,
                                   covariance_matrix = NULL,
                                   genetic_covariance_results = NULL) {
  write.csv(
    ldsc$validation_summary,
    file.path(outdir, "Python_LDSC_Input_Validation_Summary.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$heritability_results,
    file.path(outdir, "Global_Heritability_Estimates.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$genetic_correlation_results,
    file.path(outdir, "Global_Genetic_Correlations.csv"),
    row.names = FALSE
  )
  write.csv(
    ldsc$rg_se_matrix,
    file.path(outdir, "Python_LDSC_RG_SE_Matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    ldsc$intercept_se_matrix,
    file.path(outdir, "Python_LDSC_Intercept_SE_Matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    CTI,
    file.path(outdir, "GenomicPCA_CTI_Used.csv"),
    row.names = TRUE
  )
  write.csv(
    correlation_matrix,
    file.path(outdir, "GenomicPCA_Correlation_Matrix_Used.csv"),
    row.names = TRUE
  )
  write.csv(
    pca_matrix,
    file.path(outdir, "GenomicPCA_PCA_Matrix_Used.csv"),
    row.names = TRUE
  )
  if (!is.null(covariance_matrix)) {
    write.csv(
      covariance_matrix,
      file.path(outdir, "GenomicPCA_Covariance_Matrix_Used.csv"),
      row.names = TRUE
    )
  }
  if (!is.null(genetic_covariance_results)) {
    write.csv(
      genetic_covariance_results,
      file.path(outdir, "Global_Genetic_Covariances_Derived.csv"),
      row.names = FALSE
    )
  }
  write.csv(
    pc1$loading_results,
    file.path(outdir, "GenomicPCA_PC1_Weights_Used.csv"),
    row.names = FALSE
  )
  write.csv(
    pc1$eigenvalue_results,
    file.path(outdir, "GenomicPCA_Selected_Traits_Eigenvalues.csv"),
    row.names = FALSE
  )
  write.csv(
    missing_traits,
    file.path(outdir, "Python_LDSC_Dropped_Missing_Traits.csv"),
    row.names = FALSE
  )
  write.csv(
    failed_traits,
    file.path(outdir, "Python_LDSC_Dropped_Failed_Traits.csv"),
    row.names = FALSE
  )
  write.csv(
    self_pair_qc,
    file.path(outdir, "Python_LDSC_Self_Pair_QC.csv"),
    row.names = FALSE
  )
}
