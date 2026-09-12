# shared/input.R: function bodies preserved from the original workflow.

read_trait_manifest <- function(path) {
  manifest <- fread(path, data.table = FALSE, check.names = FALSE)

  if (!"traitname" %in% names(manifest)) {
    stop("The input CSV must contain a column named 'traitname'.", call. = FALSE)
  }

  traits <- trimws(as.character(manifest$traitname))

  if (anyNA(traits) || any(traits == "")) {
    stop("The traitname column contains missing or empty names.", call. = FALSE)
  }
  if (anyDuplicated(traits)) {
    duplicates <- unique(traits[duplicated(traits)])
    stop(
      "The traitname column must be unique. Duplicates:\n",
      paste(duplicates, collapse = "\n"),
      call. = FALSE
    )
  }
  if (length(traits) < 2L) {
    stop("Genomic PCA/GWAMA requires at least two traits.", call. = FALSE)
  }

  traits
}
