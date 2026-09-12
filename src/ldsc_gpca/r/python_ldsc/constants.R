required_python_ldsc_base_columns <- c(
  "p1", "p2", "rg", "se", "z", "p",
  "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"
)

heritability_column_sets <- list(
  liability = c("h2_liab", "h2_liab_se"),
  observed = c("h2_obs", "h2_obs_se")
)

internal_heritability_columns <- c("h2", "h2_se")

# All processing after input normalization uses scale-neutral h2/h2_se names.
required_python_ldsc_columns <- c(
  required_python_ldsc_base_columns,
  internal_heritability_columns
)

python_ldsc_numeric_columns <- setdiff(
  required_python_ldsc_columns,
  c("p1", "p2")
)

duplicate_comparison_columns <- c(
  "rg", "se", "z", "p", "gcov_int", "gcov_int_se"
)
