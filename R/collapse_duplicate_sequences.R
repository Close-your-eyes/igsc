#' Have a set of sequences checked for duplicates
#'
#' Duplicate sequences in the set are collapse to one sequences and named
#' as 'unique seq i'.
#'
#' @param seq_set XStringSet or character vector of sequences
#'
#' @return seq_set with unique sequences only as character vector
#' @export
#'
#' @examples
collapse_duplicate_sequences <- function(seq_set) {

 # seq_set <- as.character(seq_set)
  if (length(unique((as.character(seq_set)))) != length(as.character(seq_set)))  {
    unique_seq_set <- unique(seq_set)
    n.dup <- sapply(unique_seq_set, function(x) {length(which(x == seq_set))})
    names(unique_seq_set) <- paste("unique seq ", seq(1,length(unique_seq_set)), " (", n.dup, ")", sep = "")
    seq_set <- unique_seq_set
  } else {
    message("All sequences found unique. Will return seq_set as it is.")
  }

  return(seq_set)
}
