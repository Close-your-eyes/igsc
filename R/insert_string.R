#' Insert a string at a character position
#'
#' Inserts a string before the specified character. Positions start at 1.
#' Use position 1 to prepend, or one position past the final character
#' to append. For an empty string, use position 1.
#'
#' @param text A single, non-missing character string to modify.
#' @param insert A single, non-missing character string to insert.
#' @param x A single, non-missing, whole-number position within the string
#'   or one position past its end.
#'
#' @return A single character string containing the inserted text.
#' @export
#'
#' @examples
#' insert_string("Hello world", "beautiful ", 7)
#' insert_string("Hello", "Say: ", 1)
#' insert_string("Hello", "!", 6)
#' insert_string("", "Hello", 1)
insert_string <- function(text, insert, x) {
  stopifnot(
    is.character(text), length(text) == 1L, !is.na(text),
    is.character(insert), length(insert) == 1L, !is.na(insert),
    is.numeric(x), length(x) == 1L, !is.na(x),
    x >= 1, x <= nchar(text) + 1, x == floor(x)
  )

  paste0(
    substr(text, 1, x - 1),
    insert,
    substring(text, x)
  )
}
