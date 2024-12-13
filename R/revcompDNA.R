#' Title
#'
#' @param x
#' @param fun
#' @param rev
#' @param comp
#'
#' @return
#' @export
#'
#' @examples
revcompDNA <- function(x,
                       fun = c("rcpp", "Biostrings", "r"),
                       rev = T,
                       comp = T) {
  if (!is.character(x)) {
    stop("x has to be a character vector of DNA sequences.")
  }

  fun <- match.arg(fun, c("rcpp", "Biostrings", "r"))

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
