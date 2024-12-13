#' Title
#'
#' @param n
#' @param length
#'
#' @return
#' @export
#'
#' @examples
generate_random_dna <- function(n, length) {
  replicate(n, paste0(sample(c("A", "T", "C", "G"), length, replace = TRUE), collapse = ""))
}
