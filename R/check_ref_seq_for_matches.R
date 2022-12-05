#' Check a vector of sequences for equality on each position
#'
#' A vector of sequences is compared to a reference sequence. In case
#' of equality all sequences except for ref_seq are changed to 'match' or in case
#' of inequality to 'mismatch'.
#'
#' @param seq_set XStringSet or character vector of sequences
#' @param ref_seq_name name of reference sequence which will not be altered
#' @param pos_col name of position column
#' @param seq_col name of sequence column
#' @param name_col name of column which holds sequence names
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
check_ref_seq_for_matches <- function(seq_set,
                                      ref_seq_name,
                                      pos_col = "position",
                                      seq_col = "seq",
                                      name_col = "seq.name") {

  if (missing(ref_seq_name)) {
    stop("Please provide the name of the ref_seq in seq_set.")
  }

  if (!is.character(seq_set)) {
    seq_set <- as.character(seq_set)
  }
  if (length(unique(sapply(seq_set, nchar))) != 1) {
    stop("All sequences in seq_set have to have equal lengths.")
  }
  if (is.null(names(seq_set))) {
    stop("seq_set has to have names.")
  }
  if (!ref_seq_name %in% names(seq_set)) {
    stop("ref_seq_name not found in seq_set.")
  }

  out <- utils::stack(strsplit(seq_set, ""))
  names(out) <- c(seq_col, name_col)
  out[,pos_col] <- rep(1:nchar(seq_set[1]), length(seq_set))

  ## best way here is to use group_by from dplyr, even though a loop would also work
  out <-
    out %>%
    dplyr::group_by(!!rlang::sym(pos_col)) %>%
    dplyr::mutate(case = dplyr::case_when(nlevels(as.factor(!!rlang::sym(seq_col))) == 1 ~ "match",
                                          (nlevels(as.factor(!!rlang::sym(seq_col))) > 1 & !"-" %in% !!rlang::sym(seq_col)) ~ "mismatch",
                                          (nlevels(as.factor(!!rlang::sym(seq_col))) > 1 & "-" %in% !!rlang::sym(seq_col)) ~ "gap")) %>%
    dplyr::ungroup()
  out <- as.data.frame(out)

  out$case <- ifelse(out[,name_col] == ref_seq_name, out[,seq_col], out$case)
  out$case <- ifelse(out$case == "gap" & out[,name_col] != ref_seq_name, out[,seq_col], out$case)
  #out$case <- ifelse(out$case %in% c(names(Biostrings::IUPAC_CODE_MAP[-c(1:4)]), "+"), "ambiguous", out$case)
  out[,seq_col] <- out$case
  out[,name_col] <- factor(out[,name_col], levels = c(ref_seq_name, names(seq_set)[-length(seq_set)]))
  out <- out[,which(names(out) != "case")]
## unstack out?
  return(out)
}
