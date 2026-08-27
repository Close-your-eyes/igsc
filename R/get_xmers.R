#' Generate overlapping subsequences of fixed length
#'
#' Splits an amino acid sequence into fixed-length subsequences (x-mers)
#' with a specified overlap between consecutive subsequences.
#'
#' @param aa_sequence A single character string containing an amino acid
#'   sequence. Whitespace is removed before processing.
#' @param x A positive integer specifying the maximum length of each
#'   subsequence. Defaults to `15`.
#' @param overlap A non-negative integer specifying the number of residues
#'   shared by consecutive complete subsequences. It must be smaller than
#'   `x`. Defaults to `11`.
#' @param include_incomplete Logical. If `TRUE`, include trailing subsequences
#'   shorter than `x`. If `FALSE`, only complete x-mers are returned.
#'   Defaults to `FALSE`.
#'
#' @return A character vector containing the overlapping subsequences.
#'   If `include_incomplete = FALSE` and the sequence is shorter than `x`,
#'   `character(0)` is returned.
#'
#' @details
#' Consecutive starting positions are separated by `x - overlap` residues.
#' When `include_incomplete = TRUE`, subsequences are generated from every
#' regular starting position through the end of the sequence. Consequently,
#' multiple trailing subsequences shorter than `x` may be returned.
#'
#' @examples
#' aa_seq <- "MKTIIALSYIFCLVFADYKDDDDKGGGGS"
#'
#' # Complete 15-mers with an 11-residue overlap
#' get_xmers(aa_seq)
#'
#' # Also return trailing subsequences shorter than 15 residues
#' get_xmers(aa_seq, include_incomplete = TRUE)
#'
#' # A short sequence is returned when incomplete x-mers are enabled
#' get_xmers("MKTII", x = 10, include_incomplete = TRUE)
#'
#' @export
get_xmers <- function(
    aa_sequence,
    x = 15L,
    overlap = 11L,
    include_incomplete = FALSE
) {
  if (!is.character(aa_sequence) ||
      length(aa_sequence) != 1L ||
      is.na(aa_sequence)) {
    stop("`aa_sequence` must be a single non-missing character string.")
  }

  if (length(x) != 1L || is.na(x) || x <= 0L || x != as.integer(x)) {
    stop("`x` must be a single positive integer.")
  }

  if (length(overlap) != 1L ||
      is.na(overlap) ||
      overlap < 0L ||
      overlap != as.integer(overlap)) {
    stop("`overlap` must be a single non-negative integer.")
  }

  if (!is.logical(include_incomplete) ||
      length(include_incomplete) != 1L ||
      is.na(include_incomplete)) {
    stop("`include_incomplete` must be `TRUE` or `FALSE`.")
  }

  x <- as.integer(x)
  overlap <- as.integer(overlap)

  if (overlap >= x) {
    stop("`overlap` must be smaller than `x`.")
  }

  aa_sequence <- gsub("\\s+", "", aa_sequence)
  sequence_length <- nchar(aa_sequence)

  if (sequence_length == 0L) {
    return(character(0))
  }

  if (!include_incomplete && sequence_length < x) {
    return(character(0))
  }

  step <- x - overlap
  final_start <- if (include_incomplete) {
    sequence_length
  } else {
    sequence_length - x + 1L
  }

  starts <- seq.int(1L, final_start, by = step)

  substring(
    aa_sequence,
    first = starts,
    last = pmin(starts + x - 1L, sequence_length)
  )
}
