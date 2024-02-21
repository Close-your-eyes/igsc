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
                       trimws = F,
                       rm_comments = F,
                       comment_indicator = c(";", "#"),
                       rm_leading_arrow = T,
                       make.names = F) {

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

  if (grepl("\\.fastq\\.gz$", path)) {
    lines <- vroom::vroom_lines(gzfilefile)
  } else {
    lines <- vroom::vroom_lines(file)
  }


  if (trimws) {
    lines <- trimws(lines)
  }

  # this is tooo slow !!
  if (rm_comments) {
    message("looking for comment lines.")
    rm_inds <- unique(unlist(purrr::map(comment_indicator, function(x) {
      which(stringi::stri_startswith_fixed(pattern = x, from = 1, str = lines))
    })))
    if (length(rm_inds) > 0) {
      message(length(rm_inds), " comment lines removed")
      lines <- lines[!which(seq_along(lines) %in% rm_inds)]
    }
  }

  ind <- which(stringi::stri_startswith_fixed(pattern = ">", from = 1, str = lines))
  if (length(ind) == 0) {
    stop("fasta-formated sequences (names starting with '>') not found.")
  }
  # this is done to compensate linebreaks in sequences
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  seqs <- purrr::map_chr(igsc:::seq2(start, end), function(x) paste(lines[x], collapse = ""), .progress = T)

  names(seqs) <- gsub("^>", "", lines[ind])
  if (make.names) {
    names(seqs) <- make.names(names(seqs))
  }
  return(seqs)
}
