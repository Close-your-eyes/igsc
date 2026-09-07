.igsc_bioconductor_packages <- c(
  "Biostrings",
  "DECIPHER",
  "GenomicRanges",
  "IRanges",
  "Rsamtools",
  "biomaRt",
  "pwalign"
)

.igsc_pak_packages <- c(
  brathering = "close-your-eyes/brathering",
  colrr = "close-your-eyes/colrr",
  fcexpr = "close-your-eyes/fcexpr"
)

.igsc_cran_packages <- c(
  "Matrix",
  "Peptides",
  "RColorBrewer",
  "TSP",
  "collapse",
  "crayon",
  "forcats",
  "fst",
  "ggbeeswarm",
  "ggplot2",
  "ggrepel",
  "igraph",
  "janitor",
  "knitr",
  "openxlsx",
  "patchwork",
  "pbapply",
  "randomNames",
  "rentrez",
  "reticulate",
  "rextendr",
  "rmarkdown",
  "scales",
  "seriation",
  "stringdist",
  "tidyjson",
  "xml2"
)

.igsc_optional_packages <- c(
  .igsc_bioconductor_packages,
  names(.igsc_pak_packages),
  .igsc_cran_packages
)

# Install an optional package when first needed, then verify that it loaded.
.ensure_package <- function(package) {
  if (!package %in% .igsc_optional_packages) {
    stop("Optional package '", package, "' is not registered by igsc.", call. = FALSE)
  }

  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Installing optional package '", package, "'.")

  tryCatch(
    {
      if (package %in% .igsc_bioconductor_packages) {
        BiocManager::install(package, ask = FALSE, update = FALSE)
      } else if (package %in% names(.igsc_pak_packages)) {
        pak::pak(unname(.igsc_pak_packages[[package]]))
      } else {
        utils::install.packages(package)
      }
    },
    error = function(error) {
      stop(
        "Installation of optional package '", package, "' failed: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Optional package '", package,
      "' is required for this operation but could not be installed.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.ensure_packages <- function(packages) {
  invisible(lapply(unique(packages), igsc:::.ensure_package))
}
