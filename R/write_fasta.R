#' Write sequences to a file in fasta format
#'
#'
#'
#' @param seqs named character vector or list of sequences; if a list only provide one sequence per index
#' @param file full path to the output file to be written; recommended file extension: .fa or .fasta
#' @param mode how to write file: 'w' to initiate or overwrite; 'a' to append an existing file
#' @param linewidth number of character per line (only valid for sequence, not names)
#'
#' @return no return; file written to disk
#' @export
#'
#' @examples
write_fasta <- function(seqs,
                        file,
                        mode = c("w", "a"),
                        linewidth = 60) {

  if (missing(file)) {
    stop("Output file path has to be provided in file.")
  }

  mode <- match.arg(mode, c("w", "a"))
  if (mode == "a" && !file.exists(file)) {
    stop("file not found. set mode = 'w' or check the file path in file.")
  }

  if (is.list(seqs) && any(lengths(seqs) > 1)) {
    stop("If seqs is a list, every list entry should only contain one sequence. lengths(seqs) should be 1 for every entry.")
  }

  if (is.null(names(seqs))) {
    message("seqs have no names. Will name as 'seq_i' in order provided.")
    names(seqs) <- paste0("seq_", seq_along(seqs))
  }

  out <- file(description = file, open = mode)
  for (i in seq_along(seqs)) {
    writeLines(paste0(">", names(seqs[i])), out)

    seq <- strsplit(seqs[[i]], "")[[1]]
    seq <- split(seq, ceiling(seq_along(seq)/linewidth))
    for (j in seq_along(seq)) {
      writeLines(paste(seq[[j]], collapse = ""), out)
    }
  }
  close(out)
}
