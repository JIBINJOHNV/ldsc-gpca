# GenomicSEM argument definitions and file contracts.

genomicsem_parser <- function() {
  p <- argparse::ArgumentParser(prog = "ldsc-gpca genomicsem gpca", description = paste(
    "GPCA/GWAMA from GenomicSEM LDSCoutput RData. Strict QC is the default.",
    "Warnings, removals and errors are saved in GenomicSEM_QC_Events.csv."),
    usage = "%(prog)s --input MANIFEST.csv --ldsc_results genomicsem_LDSC.RData --outdir DIRECTORY [options]",
    formatter_class = get0("gpca_help_formatter", ifnotfound = "argparse.RawDescriptionHelpFormatter"),
    allow_abbrev = FALSE,
    epilog = paste0(
      "INPUT FILE CONTRACT\n",
      "  --input: comma-separated CSV with unique, non-empty traitname header.\n",
      "    Manifest row order defines matrix, loading and GWAMA input order.\n",
      "    Package auto-preparation also requires vcf_files if the input folder is omitted.\n",
      "  --ldsc_results: binary RData containing LDSCoutput, with S,V,I,S_Stand,V_Stand.\n",
      "    Use genomicsem_LDSC.RData from genomicsem ldsc, generated with stand=TRUE.\n",
      "    This is NOT a delimited table; the first-pass _raw.RData is not GPCA input.\n",
      "  --gpca_input_folder: tab-separated TSV files; required columns (extra/reordered columns allowed):\n",
      "    SNPID,CHR,BP,EA,OA,EAF,N,Z,P (A1/A2/p aliases accepted for EA/OA/P).\n",
      "    Split: {traitname}_chr{CHR}_GenomicPCA_inputs.tsv (chromosomes 1-22).\n",
      "    Whole genome: {traitname}_GenomicPCA_inputs.tsv.\n",
      "  --source_path: optional custom modified GWAMA R source; default is bundled v1.2.6.\n\n",
      "GWAMA OUTPUT QC\n",
      "  Zero available weight or invalid BETA/Z/P/SE/N fails before success/export.\n",
      "  P=0 underflow is valid. Raw results and *.GWAMA_QC_* audits are retained.\n",
      "  No undefined row is silently excluded or assigned a fabricated Z.\n\n",
      "AVAILABILITY\n",
      "  genomicsem ldsc runs GenomicSEM munge/LDSC or accepts existing munged files.\n",
      "  Both PCA matrix choices write all-PC variance and PC1 protein contributions.\n",
      "  GWAMA uses PC1 only. GenomicSEM covariance uses S on its supplied scales.\n"))
  required <- p$add_argument_group("Required inputs")
  inputs <- p$add_argument_group("GWAMA inputs (not required with --validate_only)")
  traits <- p$add_argument_group("Trait and failed-result handling")
  pca <- p$add_argument_group("PCA settings")
  diagnostics <- p$add_argument_group("Diagnostic warning/error controls")
  execution <- p$add_argument_group("Execution and performance")
  required$add_argument("--input", required = TRUE, metavar = "MANIFEST.csv", help = "Comma-separated manifest; required header traitname. Required; no default.")
  required$add_argument("--ldsc_results", required = TRUE, metavar = "genomicsem_LDSC.RData", help = "GenomicSEM LDSCoutput RData; required objects below. Required; no default.")
  required$add_argument("--outdir", required = TRUE, metavar = "DIRECTORY", help = "Output directory; use a fresh directory. Required; no default.")
  inputs$add_argument("--gpca_input_folder", metavar = "DIRECTORY", help = "Tab-separated GWAMA input files. Default: unset; package CLI prepares from VCF unless --validate_only.")
  inputs$add_argument("--source_path", default = bundled_gwama_path, metavar = "GWAMA.R", help = "Optional custom modified GWAMA R source. Default: bundled N_weighted_GWAMA.function.1_2_6.R.")
  inputs$add_argument("--splitby_chr", choices = c("split", "nosplit"), default = "split", help = "Input naming mode. Default: split (chromosomes 1-22).")
  execution$add_argument("--n_cores", dest = "n_cores", type = "integer", default = 0L, help = "Chromosome workers. Default: 0 = auto; Windows runs sequentially. Failed jobs retry once, one at a time; failure after 2 attempts stops analysis.")
  execution$add_argument("--validate_only", action = "store_true", default = FALSE, help = "Write QC and PCA diagnostics without running GWAMA. Default: off.")
  traits$add_argument("--allow_missing_traits", action = "store_true", default = FALSE, help = "Remove absent manifest traits with reasons. Default: off; stop if absent.")
  traits$add_argument("--failed_ldsc_action", choices = c("error", "drop_traits"), default = "error", help = "Stop on invalid trait/pair estimates or remove traits with an audit; never impute.")
  diagnostics$add_argument("--h2_z_warn_threshold", type = "double", default = 2, help = "Warn for retained h2/SE below this value; never removes traits. Default: 2; 0 disables.")
  diagnostics$add_argument("--rg_out_of_range_action", choices = c("warn", "error"), default = "warn", help = "Warn or stop for finite off-diagonal abs(rg)>1; never clamp.")
  pca$add_argument("--pca_matrix", choices = c("correlation", "covariance"), default = "correlation", help = "correlation uses S_Stand; covariance uses GenomicSEM S, on its original scales.")
  pca$add_argument("--pc1_orientation", choices = c("tutorial", "as_computed"), default = "tutorial", help = "tutorial flips ALL PC1 loadings if median<0; as_computed keeps the eigenvector sign.")
  diagnostics$add_argument("--negative_eigen_action", choices = c("warn", "error"), default = "warn", help = "Warn or stop for substantive negative PCA eigenvalues; no matrix repair.")
  diagnostics$add_argument("--matrix_eigen_tolerance", type = "double", default = 1e-8, help = "Relative negative-eigenvalue threshold. Default: 1e-8.")
  p
}

# One event per reason; pair events name BOTH traits. Run/matrix events have
# no Trait: they must not imply a specific trait caused a matrix-wide problem.
