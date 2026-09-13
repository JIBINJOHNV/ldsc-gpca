# Native GenomicSEM argument definitions and file contracts.

genomicsem_parser <- function() {
  p <- argparse::ArgumentParser(prog = "ldsc-gpca genomicsem gpca", description = paste(
    "GPCA/GWAMA from GenomicSEM LDSCoutput RData. Strict QC is the default.",
    "Warnings, removals and errors are saved in GenomicSEM_QC_Events.csv."),
    formatter_class = "argparse.RawDescriptionHelpFormatter",
    epilog = paste0(
      "INPUT FILE CONTRACT\n",
      "  --input: comma-separated CSV with unique, non-empty traitname header.\n",
      "    Manifest row order defines matrix, loading and GWAMA input order.\n",
      "    Package auto-preparation also requires vcf_files if the input folder is omitted.\n",
      "  --ldsc_path: binary RData containing LDSCoutput, with S,V,I,S_Stand,V_Stand.\n",
      "    This is NOT a delimited table. Generate it upstream with stand=TRUE.\n",
      "  --gpca_input_folder: tab-separated TSV files; exact column order:\n",
      "    SNPID,CHR,BP,EA,OA,EAF,N,Z,P (A1/A2/p aliases accepted for EA/OA/P).\n",
      "    Split: {traitname}_chr{CHR}_GenomicPCA_inputs.tsv (chromosomes 1-22).\n",
      "    Whole genome: {traitname}_GenomicPCA_inputs.tsv.\n",
      "  --source_path: optional custom modified GWAMA R source; default is bundled v1.2.6.\n\n",
      "AVAILABILITY\n",
      "  genomicsem ldsc runs native munge/LDSC or accepts existing munged files.\n",
      "  Both PCA matrix choices write all-PC variance and PC1 protein contributions.\n",
      "  GWAMA uses PC1 only. Native covariance uses S on its supplied scales.\n"))
  p$add_argument("--input", required = TRUE, metavar = "MANIFEST.csv", help = "Comma-separated manifest; required header traitname. Required; no default.")
  p$add_argument("--ldsc_path", required = TRUE, metavar = "LDSC.RData", help = "Native LDSCoutput RData; required objects below. Required; no default.")
  p$add_argument("--outdir", required = TRUE, metavar = "DIRECTORY", help = "Output directory; use a fresh directory. Required; no default.")
  p$add_argument("--gpca_input_folder", metavar = "DIRECTORY", help = "Tab-separated GWAMA input files. Default: unset; package CLI prepares from VCF unless --validate_only.")
  p$add_argument("--source_path", default = bundled_gwama_path, metavar = "GWAMA.R", help = "Optional custom modified GWAMA R source. Default: bundled N_weighted_GWAMA.function.1_2_6.R.")
  p$add_argument("--splitby_chr", choices = c("split", "nosplit"), default = "split", help = "Input naming mode. Default: split (chromosomes 1-22).")
  p$add_argument("--cores", type = "integer", default = 0L, help = "Chromosome workers. Default: 0 = auto; Windows runs sequentially.")
  p$add_argument("--validate_only", action = "store_true", default = FALSE, help = "Write QC and PCA diagnostics without running GWAMA. Default: false.")
  p$add_argument("--allow_missing_traits", action = "store_true", default = FALSE, help = "Remove absent manifest traits with reasons. Default: stop if absent.")
  p$add_argument("--failed_ldsc_action", choices = c("error", "drop_traits"), default = "error", help = "Invalid trait/pair estimates: error (default), or audited deterministic trait removal; never impute.")
  p$add_argument("--h2_z_warn_threshold", type = "double", default = 2, help = "Warn for retained h2/SE below this value; never removes traits. Default: 2; 0 disables.")
  p$add_argument("--rg_out_of_range_action", choices = c("warn", "error"), default = "warn", help = "Finite off-diagonal abs(rg)>1: warn (default) or error; never clamp.")
  p$add_argument("--pca_matrix", choices = c("correlation", "covariance"), default = "correlation", help = "correlation uses S_Stand (default); covariance uses native S, on its original scales.")
  p$add_argument("--pc1_orientation", choices = c("tutorial", "as_computed"), default = "tutorial", help = "tutorial (default): flip ALL PC1 loadings if median<0; as_computed keeps eigenvector sign.")
  p$add_argument("--negative_eigen_action", choices = c("warn", "error"), default = "warn", help = "Substantive negative PCA eigenvalues: warn (default) or error; no matrix repair.")
  p$add_argument("--matrix_eigen_tolerance", type = "double", default = 1e-8, help = "Relative negative-eigenvalue threshold. Default: 1e-8.")
  p
}

# One event per reason; pair events name BOTH traits. Run/matrix events have
# no Trait: they must not imply a specific trait caused a matrix-wide problem.
