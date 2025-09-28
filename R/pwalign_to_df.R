#' Convert pairwiseAlignment object to XStringSet
#'
#' E.g. for printing. See DECIPHER::BrowseSeqs.
#'
#' @param pa pairwiseAlignment from pwalign::pairwiseAlignment
#' @param verbose print messages?
#' @param subject_width include whole subject width or pattern-aligned range only
#'
#' @returns XStringSet
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
pwalign_to_xstringset <- function(pa,
                                  verbose = T,
                                  subject_width = c("whole", "aligned")) {


  if (length(pa) > 1) {
    message("pwalign_to_xstringset: Using first pairwiseAlignment only.")
    pa <- pa[[1]]
  }

  subject_width <- rlang::arg_match(subject_width)

  type <- pwalign_guess_type(pa)
  xstringfun <- switch(type,
                       DNA = Biostrings::DNAStringSet,
                       RNA = Biostrings::RNAStringSet,
                       AA = Biostrings::AAStringSet)


  if (any(c(is.null(pa@pattern@unaligned@ranges@NAMES), is.null(pa@subject@unaligned@ranges@NAMES))) && verbose) {
    message("No names provided for pattern and/or subject. If desired, provide named XStringSets, each of length 1, as subject/pattern to pwalign::pairwiseAlignment.")
  }

  pattern_name <- pa@pattern@unaligned@ranges@NAMES %||% "pattern"
  subject_name <- pa@subject@unaligned@ranges@NAMES %||% "subject"

  if (subject_width == "whole") {
    subject <- pwalign::alignedSubject(pa)
    pattern <- pwalign::alignedPattern(pa)
  } else if (subject_width == "aligned") {
    pattern <- pa@pattern
    subject <- pa@subject
  }

  vec <- stats::setNames(c(as.character(pattern), as.character(subject)),
                         nm = c(pattern_name, subject_name))

  xstr <- xstringfun(vec)
  attr(xstr, "subject_name") <- subject_name

  return(xstr)
}



#' Convert pairwiseAlignment object to long data frame
#'
#' Wrapper around pwalign_to_xstringset and xstringset_to_df.
#' df may be used for plotting. Or passed to compare_seq_df_long.
#'
#' @param pa pairwiseAlignment from pwalign::pairwiseAlignment
#' @param verbose print messages?
#' @param subject_width include whole subject width or pattern-aligned range only
#'
#' @returns long data frame
#' @export
#'
#' @examples
pwalign_to_df <- function(pa,
                          verbose = T,
                          subject_width = c("whole", "aligned")) {

  xstr <- pwalign_to_xstringset(pa = pa,
                                verbose = verbose,
                                subject_width = subject_width)

  return(xstringset_to_df(
    xstringset = xstr,
    subject_name = attr(xstr, "subject_name")
  ))
}



#' Title
#'
#' @param df
#' @param type
#' @param name_col
#' @param seq_col
#' @param sym_internal_gap
#' @param sym_terminal_gap
#'
#' @returns
#' @export
#'
#' @examples
df_to_xstringset <- function(df,
                             type = c("DNA", "RNA", "AA"),
                             name_col = "seq.name",
                             seq_col = "seq",
                             sym_internal_gap = "-",
                             sym_terminal_gap = "-") {

  type <- rlang::arg_match(type)
  if (!name_col %in% names(df)) {
    stop("name_col not found in df.")
  }
  if (!seq_col %in% names(df)) {
    stop("seq_col not found in df.")
  }
  df <-
    df %>%
    dplyr::mutate({{seq_col}} := ifelse(is.na(!!rlang::sym(seq_col)), sym_internal_gap, !!rlang::sym(seq_col))) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col), values_fill = sym_terminal_gap)
  df <- df[,-1]
  df <- unlist(lapply(df, paste, collapse = ""))

  if (type == "DNA") {
    return(Biostrings::DNAStringSet(x = df))
  } else if (type == "RNA") {
    return(Biostrings::RNAStringSet(x = df))
  } else if (type == "AA") {
    return(Biostrings::AAStringSet(x = df))
  }
}


guess_type2 <- function(seq_vector) {
  seq_vector <- unique(toupper(as.character(seq_vector)))
  seq_vector <- gsub("[^A-Z]", "", seq_vector)

  if (all(gsub("[ATCG]", "", seq_vector) == "")) {
    return("DNA")
  }
  if (all(gsub("[AUCG]", "", seq_vector) == "")) {
    return("RNA")
  }
  if (all(gsub(paste0("[",gsub("[^A-Z]", "", paste(Biostrings::DNA_ALPHABET, collapse = "")), "]"), "", seq_vector) == "")) {
    return("DNA")
  }
  if (all(gsub(paste0("[",gsub("[^A-Z]", "", paste(Biostrings::RNA_ALPHABET, collapse = "")), "]"), "", seq_vector) == "")) {
    return("RNA")
  }
  if (all(gsub(paste0("[",gsub("[^A-Z]", "", paste(Biostrings::AA_ALPHABET, collapse = "")), "]"), "", seq_vector) == "")) {
    return("AA")
  }

  stop("Alignment type could not be determined unambigously. Please define it: 'NT' or 'AA'.")
}

pwalign_guess_type <- function(pa) {
  if (methods::is(pa@pattern@unaligned, "QualityScaledDNAStringSet")) {
    type <- "DNA"
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledRNAStringSet")) {
    type <- "RNA"
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledAAStringSet")) {
    type <- "AA"
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledBStringSet")) {
    type <- guess_type2(c(pa@subject, pa@pattern))
  }
  return(type)
}
