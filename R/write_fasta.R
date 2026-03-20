#' Write sequences to a file in fasta format
#'
#'
#'
#' @param seqs named character vector or list of sequences; if a list only provide one sequence per index
#' @param file full path to the output file to be written; recommended file extension: .fa or .fasta
#' @param append append? if F will not overwrite.
#' @param linewidth number of character per line (only valid for sequence, not names)
#'
#' @return no return; file written to disk
#' @export
#'
#' @examples
write_fasta <- function(seqs,
                        file,
                        linewidth = 60,
                        mc.cores = floor(parallel::detectCores()/4)) {

  if (missing(file)) {
    stop("Output file path has to be provided in file.")
  }

  if (is.list(seqs) && any(lengths(seqs) > 1)) {
    stop("If seqs is a list, every list entry should only contain one sequence. lengths(seqs) should be 1 for every entry.")
  }

  if (is.null(names(seqs))) {
    message("seqs have no names. Will name as 'seq_i' in order provided.")
    names(seqs) <- paste0("seq_", seq_along(seqs))
  }

  # lines <- purrr::map2(
  #   names(seqs),
  #   seqs,
  #   ~ c(paste0(">", .x), split_chunks(.y, linewidth))
  # ) |>
  #   unlist(use.names = FALSE)

  lines <- parallel::mcmapply(FUN = function(name, seq) c(paste0(">", name), split_chunks(seq, linewidth)),
    name = names(seqs),
    seq  = seqs,
    SIMPLIFY = FALSE,
    mc.cores = mc.cores)
  lines <- unlist(lines, use.names = FALSE)

  vroom::vroom_write_lines(lines, file = file)
  message(file)
}

split_chunks <- function(x, n) {
  starts <- seq(1L, stringi::stri_length(x), by = n)
  ends   <- pmin(starts + n - 1L, stringi::stri_length(x))
  return(stringi::stri_sub(x, starts, ends))
}
