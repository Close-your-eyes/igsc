#' Generate random DNA sequences
#'
#' Generate one or more DNA strings by independently sampling the canonical
#' bases `A`, `T`, `C`, and `G` with equal probability and replacement.
#'
#' This function uses R's global random-number generator. Call [base::set.seed()]
#' before generation when reproducible sequences are required.
#'
#' @param n Positive integer giving the number of sequences to generate.
#' @param length Non-negative integer giving the number of bases in each
#'   sequence.
#'
#' @return A character vector of `n` random DNA sequences, each containing
#'   `length` bases.
#' @export
#'
#' @examples
#' set.seed(42)
#' generate_random_dna(n = 3, length = 12)
generate_random_dna <- function(n, length) {
  replicate(n, paste0(sample(c("A", "T", "C", "G"), length, replace = TRUE), collapse = ""))
}
