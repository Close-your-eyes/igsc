#' Load amino acid substitution matrices
#'
#' Loads BLOSUM and PAM matrices from the `pwalign` package into the global
#' environment, then reports which matrices are available.
#'
#' @details Existing objects in the global environment with the same names
#'   will be replaced.
#' @return Called for its side effect of loading matrices.
#' @export
#'
#' @examples
#' \dontrun{
#' get_substitution_matrices()
#' BLOSUM62["A", "G"]
#' }
get_substitution_matrices <- function() {

  igsc:::.ensure_package("pwalign")

  matrices <- c(
    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
    "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"
  )

  data(list = matrices, package = "pwalign", envir = .GlobalEnv)
  message("Matrices now available in the global environment: ",
          paste(matrices, collapse = ", "))
}
