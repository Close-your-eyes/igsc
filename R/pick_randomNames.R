#' Pick n random names from the randomNames package
#'
#' The randomNames function from randomNames does not necessarily provide unique names (https://github.com/CenterForAssessment/randomNames/issues/55).
#' Hence this function has been written to fix that issue as elegant as possible.
#'
#' @param n numeric (integer) of how many unique names to return
#' @param names_to_avoid a vector of names to avoid
#' @param max_iter maximum number of while-loop iterations to avoid an infinite loop
#' @param min_name_nchar min length of name returned by randomNames
#' @param avoid_chars character vector of forbidden characters in names or NULL;
#' reject names with one or more of these characters
#' @param ... arguments passed to randomNames::randomNames
#'
#' @return a vector of n unique random names
#' @export
#'
#' @examples
#' pick_randomNames(n = 1000, names_to_avoid = c("Chris", "Diana", "Leonie), which.names = "first")
pick_randomNames <- function(n,
                             names_to_avoid = NULL,
                             max_iter = 500,
                             min_name_nchar = 3,
                             avoid_chars = c(" ", "-", "'", ",", ";"),
                             randomNames_args = list(which.names = "first")) {

  if ("which.names" %in% names(randomNames_args)) {
    if (randomNames_args$which.names == "both") {
      if ("names.sep" %in% names(randomNames_args)) {
        sep <- randomNames_args$name.sep
      } else {
        sep <- formals(randomNames::randomNames)$name.sep
      }
      if (length(intersect(trimws(sep), avoid_chars)) > 0) {
        collision <- intersect(trimws(sep), avoid_chars)
        stop("name.sep and avoid_chars are colliding: ", collision, ". Please change name.sep or remove this character from avoid_chars.")
      }
    }
  }

  names_to_avoid <- names_to_avoid[which(!is.na(names_to_avoid))]
  names_to_avoid <- trimws(names_to_avoid)
  names <- trimws(unique(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args, list(n = n)))))
  #names <- trimws(unique(randomNames::randomNames(n = n, ...)))
  names <- names[which(!names %in% names_to_avoid)]
  names <- names[sapply(names, nchar, simplify = T) >= min_name_nchar]
  if (!is.null(avoid_chars)) {
    names <- names[which(!grepl(paste(avoid_chars, collapse = "|"), names))]
  }

  iters <- 0
  while (length(names) < n && iters <= max_iter) {
    temp <- trimws(unique(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args, list(n = n - length(names))))))
    #temp <- trimws(unique(randomNames::randomNames(n = n - length(names), ...)))
    temp <- temp[which(!temp %in% c(names_to_avoid, names))]
    temp <- temp[sapply(temp, nchar, simplify = T) >= min_name_nchar]
    if (!is.null(avoid_chars)) {
      temp <- temp[which(!grepl(paste(avoid_chars, collapse = "|"), temp))]
    }
    names <- c(names, temp)
    iters <- iters + 1
  }
  if (length(names) < n) {
    stop(paste0("max_iter reached before n unique random names could be picked. ", length(names), " unique names were reached. Either increase max_iter or reduce stringency of name selection (more genders, more ethnicities etc.)."))
  }
  return(names)
}
