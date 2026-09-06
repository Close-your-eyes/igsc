#' Write sequences to a FASTA file
#'
#' Write character sequences as FASTA records, with one header per sequence and
#' sequence text wrapped to a configurable line width. Record formatting is
#' parallelized across sequences with [parallel::mcmapply()].
#'
#' Names of `seqs` become FASTA headers. If `seqs` is unnamed, headers are
#' generated as `seq_1`, `seq_2`, and so on, and a message is emitted. Existing
#' files at `file` may be replaced by [vroom::vroom_write_lines()]. The parent
#' directory is not created automatically.
#'
#' When `gzip = TRUE`, the external `gzip` command compresses the completed
#' file. On success, the original uncompressed file is removed and the output is
#' available at `paste0(file, ".gz")`.
#'
#' @param seqs Character vector or list containing the sequences to write. If a
#'   list is supplied, every element must contain no more than one sequence.
#'   Names, when present, are used as FASTA headers.
#' @param file Output-file path. Extensions such as `.fa` or `.fasta` are
#'   recommended but not required.
#' @param linewidth Positive integer giving the maximum number of sequence
#'   characters per output line. Header lines are not wrapped.
#' @param mc.cores Positive integer giving the number of worker processes used
#'   by [parallel::mcmapply()] while formatting records.
#' @param gzip Logical; compress the written file with the external `gzip`
#'   command.
#' @param verbose Logical; report the requested output path after writing.
#'
#' @return `NULL`, invisibly. This function is called for its file-writing side
#'   effect.
#' @export
#'
#' @examples
#' sequences <- c(
#'   alpha = "AACCGGTTAACCGGTT",
#'   beta = "TTTTCCCCAAAAGGGG"
#' )
#' fasta_file <- tempfile(fileext = ".fa")
#'
#' write_fasta(
#'   seqs = sequences,
#'   file = fasta_file,
#'   linewidth = 8,
#'   mc.cores = 1,
#'   verbose = FALSE
#' )
#' readLines(fasta_file)
#'
#' unlink(fasta_file)
write_fasta <- function(seqs,
                        file,
                        linewidth = 60,
                        mc.cores = floor(parallel::detectCores()/4),
                        gzip = F,
                        verbose = T) {

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

  lines <- parallel::mcmapply(FUN = function(name, seq) c(paste0(">", name),
                                                          split_chunks(seq, linewidth)),
                              name = names(seqs),
                              seq  = seqs,
                              SIMPLIFY = FALSE,
                              mc.cores = mc.cores)
  lines <- unlist(lines, use.names = FALSE)

  vroom::vroom_write_lines(lines, file = file)
  if (gzip) {
    system(paste0("gzip ", file))
  }
  if (verbose) {
    message(file)
  }
}

split_chunks <- function(x, n) {
  starts <- seq(1L, stringi::stri_length(x), by = n)
  ends   <- pmin(starts + n - 1L, stringi::stri_length(x))
  return(stringi::stri_sub(x, starts, ends))
}
