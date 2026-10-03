# Publish the Python LDSC matrices and analysis audit tables.

write_analysis_outputs <- function(outdir, ldsc, pc1, correlation_matrix,
                                   pca_matrix, CTI,
                                   missing_traits, failed_traits,
                                   self_pair_qc,
                                   covariance_matrix = NULL,
                                   genetic_covariance_results = NULL) {
  tables <- list(Python_LDSC_Input_Validation_Summary = ldsc$validation_summary,
    Global_Heritability_Estimates = ldsc$heritability_results,
    Global_Genetic_Correlations = ldsc$genetic_correlation_results,
    GenomicPCA_PC1_Weights_Used = pc1$loading_results,
    GenomicPCA_Selected_Traits_Eigenvalues = pc1$eigenvalue_results,
    Python_LDSC_Dropped_Missing_Traits = missing_traits,
    Python_LDSC_Dropped_Failed_Traits = failed_traits,
    Python_LDSC_Self_Pair_QC = self_pair_qc)
  matrices <- list(Python_LDSC_RG_SE_Matrix = ldsc$rg_se_matrix,
    Python_LDSC_Intercept_SE_Matrix = ldsc$intercept_se_matrix,
    GenomicPCA_CTI_Used = CTI, GenomicPCA_Correlation_Matrix_Used = correlation_matrix,
    GenomicPCA_PCA_Matrix_Used = pca_matrix)
  if (!is.null(covariance_matrix)) matrices$GenomicPCA_Covariance_Matrix_Used <- covariance_matrix
  if (!is.null(genetic_covariance_results)) tables$Global_Genetic_Covariances_Derived <- genetic_covariance_results
  for (name in names(tables)) write.csv(tables[[name]],
    file.path(outdir, paste0(name, ".csv")), row.names = FALSE)
  for (name in names(matrices)) write.csv(matrices[[name]],
    file.path(outdir, paste0(name, ".csv")), row.names = TRUE)
}
