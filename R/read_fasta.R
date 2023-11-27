#' Read sequences from a fasta-formatted file
#'
#' vroom::vroom_lines is used to quickly read lines of a text file (as a fasta file is one).
#'
#'
#' @param file path to a file with fasta-formatted sequences
#' @param rm_comments logical whether to remove comment lines
#' @param comment_indicator first characters of a line which indicate a comment
#'
#' @return a named character vector of sequences in file
#' @export
#'
#' @examples
#' \dontrun{
#' seqs <- igsc::read_fasta(file = "my/fasta/file.fasta")
#' }
read_fasta <- function(file,
                       rm_comments = T,
                       comment_indicator = c(";", "#")) {

  if (missing(file)) {
    stop("Please provide a path to a file in 'file'.")
  }

  if (!is.character(file)) {
    stop("file has to be the path to a fasta-formatted file.")
  }

  if (!file.exists(file)) {
    stop("file not found.")
  }

  # how to handle very large fasta files?
  # offer to avoid concatenation ?

  lines <- vroom::vroom_lines(file)

  if (rm_comments) {
    rm_regex <- paste(paste0("^", comment_indicator), collapse = "|")
    rm <- !grepl(rm_regex, trimws(lines))
    if (any(!rm)) {
      message(sum(which(!rm)), " comment lines removed")
      lines <- lines[rm]
    }
  }

  ind <- grep("^>", trimws(lines))
  if (length(ind) == 0) {
    stop("fasta-formated sequences (names starting with '>') not found.")
  }

  # this is done to compensate linebreaks in sequences
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  seqs <- mapply(x = start, y = end, function(x,y) paste(lines[x:y], collapse = ""), SIMPLIFY = T)
  names(seqs) <- gsub("^>", "", trimws(lines[ind]))

  return(seqs)
}
