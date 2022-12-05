#' Get the longest consecutive stretch of distinct sequence
#'
#' Use case: Have a consensus sequence calculated from a set of DNA or AA sequences. Then
#' get the longest consecutive stretch of distinct sequence from the consensus.
#'
#' @param seq Sequence as character, DNAStringSet, DNAString, RNAStringSet, RNAString,
#' AAStringSet or AAString
#' @param seq_type set sequence type to AA or NT if necessary; if NULL
#' it is attempted to guess the type
#'
#' @return a character sequence
#' @export
#'
#' @examples
consecutive_distinct_seq <- function (seq,
                                      seq_type = NULL) {
  if (length(seq) > 1) {
    stop("Please provide one sequence, XString, XStringSet or character.")
  }
  if (is.null(seq_type)) {
    if (methods::is(seq, "DNAStringSet") || methods::is(seq, "RNAStringSet") || methods::is(seq, "AAStringSet") || methods::is(seq, "DNAString") || methods::is(seq, "RNAString") || methods::is(seq, "AAString")) {
      if (methods::is(seq, "DNAStringSet") || methods::is(seq, "RNAStringSet") || methods::is(seq, "DNAString") || methods::is(seq, "RNAString")) {
        seq_type <- "NT"
      }
      if (methods::is(seq, "AAStringSet") || methods::is(seq, "AAString")) {
        seq_type <- "AA"
      }
    } else {
      seq_type <- guess_type(seq)
    }
  }

  if (seq_type == "NT") {
    ref_symbols <- Biostrings::IUPAC_CODE_MAP[-c(1:4)]
  } else if (seq_type == "AA") {
    ref_symbols <- setdiff(strsplit(seq, "")[[1]], Biostrings::AA_PROTEINOGENIC)
  }

  out_seq <- names(which.max(sapply(strsplit(as.character(seq), paste(c(names(ref_symbols)), collapse = "|"))[[1]], nchar)))
  pa <- Biostrings::pairwiseAlignment(subject = seq, pattern = out_seq, type = "local")
  pa_limits <- c(pa@subject@range@start, pa@subject@range@start + pa@subject@range@width - 1)

  return(list(seq = out_seq, limits = pa_limits))
}
