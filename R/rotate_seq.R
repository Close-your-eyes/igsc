#' Rotate character strings by splitting at a given position
#'
#' This function rotates each string in a character vector by splitting at position \code{cut} (1-indexed) and swapping the two parts.
#' For example, "abcdef" with \code{cut = 3} becomes "cdefab".
#'
#' The function will stop if \code{cut} is less than 1 or greater than the length of any string in \code{seq}.
#'
#' @param seq A character vector where each element is a string to be rotated.
#' @param cut A single integer between 1 and the length of each string in \code{seq} (inclusive). Defaults to 1.
#'
#' @return A character vector of the same length and names as \code{seq}, with each string rotated.
#'
#' @examples
#' # Basic rotation
#' rotate_seq("abcdef", cut = 3)
#' # [1] "cdefab"
#'
#' # Vector of strings
#' rotate_seq(c("abc", "defg", "hijkl"), cut = 2)
#' # [1] "cabc" "fgde" "ijkhl"
#'
#' # Named input preserves names
#' rotate_seq(c(a = "abcdef", b = "ghijkl"), cut = 3)
#' #        a        b
#' # "cdefab" "ijklgh"
#'
#' @export
rotate_seq <- function(seq, cut = 1) {

  lens <- nchar(seq)

  if (cut == 1) {
    return(seq)  # no rotation needed
  }

  if (cut < 1 || any(cut > lens)) {
    stop("cut must be between 1 and string length for all elements in 'seq'")
  }

  part1 <- substr(seq, cut, lens)
  part2 <- substr(seq, 1, cut - 1)

  return(stats::setNames(paste0(part1, part2), names(seq)))
}
