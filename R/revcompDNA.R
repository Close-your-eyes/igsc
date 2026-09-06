#' Reverse and/or Complement DNA Sequences
#'
#' Reverse, complement, or reverse-complement a character vector of DNA
#' sequences using an Rcpp, Biostrings, or base R implementation.
#'
#' @param x A character vector containing DNA sequences. Sequences are
#'   converted to uppercase before processing.
#' @param fun Character string selecting the implementation. One of
#'   \code{"rcpp"} (the default), \code{"Biostrings"}, or \code{"r"}.
#' @param rev A single logical value indicating whether each sequence should
#'   be reversed.
#' @param comp A single logical value indicating whether each sequence should
#'   be complemented. The canonical substitutions are \code{A <-> T} and
#'   \code{C <-> G}.
#'
#' @return A character vector containing the transformed DNA sequences.
#'
#' @details
#' When both \code{rev} and \code{comp} are \code{TRUE}, the reverse
#' complement is returned. When only \code{rev} is \code{TRUE}, the sequences
#' are reversed. When only \code{comp} is \code{TRUE}, the sequences are
#' complemented.
#'
#' The \code{"Biostrings"} implementation requires the Biostrings package.
#' Handling of non-canonical nucleotide symbols may differ between
#' implementations.
#'
#' @export
#'
#' @examples
#' revcompDNA("ATGC")
#'
#' revcompDNA(c("ATGC", "AATT"))
#'
#' # Reverse without complementing
#' revcompDNA("ATGC", rev = TRUE, comp = FALSE)
#'
#' # Complement without reversing
#' revcompDNA("ATGC", rev = FALSE, comp = TRUE)
revcompDNA <- function(x,
                       fun = c("rcpp", "Biostrings", "r"),
                       rev = T,
                       comp = T) {
  if (!is.character(x)) {
    stop("x has to be a character vector of DNA sequences.")
  }

  fun <- rlang::arg_match(fun)

  x <- toupper(x)

  if (fun == "Biostrings") {
    if (rev && comp)  {
      x <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
    } else if (rev) {
      x <- as.character(Biostrings::reverse(Biostrings::DNAStringSet(x)))
    } else if (comp) {
     x <- as.character(Biostrings::complement(Biostrings::DNAStringSet(x)))
    }
    return(x)
  }

  if (fun == "r") {
    unlist(lapply(x, function(y) {
      # N remains N, as any other character
      if (rev) {
        y <- rev(strsplit(y, "")[[1]])
      } else {
        y <- strsplit(y, "")[[1]]
      }
      if (comp) {
        AtoT <- which(y == "A")
        TtoA <- which(y == "T")
        CtoG <- which(y == "C")
        GtoC <- which(y == "G")
        y[AtoT] <- "T"
        y[TtoA] <- "A"
        y[CtoG] <- "G"
        y[GtoC] <- "C"
      }
    }))
    return(paste(y, collapse = ""))
  }

  if (fun == "rcpp") {
    if (rev && comp)  {
      mode <- "both"
    } else if (rev) {
      mode <- "reverse"
    } else if (comp) {
      mode <- "complement"
    }
    xnames <- names(x)
    return(stats::setNames(igsc:::revcomp_rcpp2(dna_strings = x, mode = mode), xnames))
  }

  # examples:
  # https://github.com/r-lib/bench/issues/59
  # random_dna <- generate_random_dna(n = 5, length = 10000)
  # bench::mark(bio = revcompDNA(random_dna, fun = "Biostrings"),
  #             R = revcompDNA(random_dna, fun = "r"),
  #             Rcpp = revcompDNA(random_dna, fun = "rcpp"), iterations = 200, check = F)

}
