#' Pick n random names from the randomNames package
#'
#' The randomNames function from randomNames does not necessarily provide unique names (https://github.com/CenterForAssessment/randomNames/issues/55).
#' Hence this function has been written to fix that issue as elegant as possible.
#'
#' @param n numeric (integer) of how many unique names to return
#' @param names_to_avoid a vector of names to avoid
#' @param max_iter maximum number of while-loop iterations to avoid an infinite loop
#' @param ... arguments passed to randomNames::randomNames
#'
#' @return a vector of n unique random names
#' @export
#'
#' @examples
#' pick_randomNames(n = 1000, names_to_avoid = c("Chris", "Diana", "Leonie), which.names = "first")
pick_randomNames <- function(n, names_to_avoid = NULL, max_iter = 500, ...) {
  names_to_avoid <- names_to_avoid[which(!is.na(names_to_avoid))]
  names_to_avoid <- trimws(names_to_avoid)
  names <- trimws(unique(randomNames::randomNames(n = n, ...)))
  names <- names[which(!names %in% names_to_avoid)]
  iters <- 0
  while (length(names) < n && iters <= max_iter) {
    temp <- trimws(unique(randomNames::randomNames(n = n - length(names), ...)))
    temp <- temp[which(!temp %in% c(names_to_avoid, names))]
    names <- c(names, temp)
    iters <- iters + 1
  }
  if (length(names) < n) {
    stop(paste0("max_iter reached before n unique random names could be picked. ", length(names), " unique names were reached. Either increase max_iter or reduce stringency of name selection (more genders, more ethnicities etc.)."))
  }
  return(names)
}
