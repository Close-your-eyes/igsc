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
                                      ref_seq_name = NULL,
                                      pos_col = "position",
                                      seq_col = "seq",
                                      name_col = "seq.name") {

  if (is.null(ref_seq_name)) {
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

  seq_df <- utils::stack(strsplit(seq_set, ""))
  names(seq_df) <- c(seq_col, name_col)
  seq_df[,pos_col] <- rep(1:nchar(seq_set[1]), length(seq_set))
  seq_df <- as.data.frame(tidyr::pivot_wider(seq_df, names_from = seq.name, values_from = seq))

  ## best way here is to use group_by from dplyr, even though a loop would also work
  for (x in names(seq_df)[which(!names(seq_df) %in% c("position", ref_seq_name))]) {
    seq_df[,x] <- ifelse(seq_df[,x] == seq_df[,ref_seq_name], "match", ifelse(seq_df[,x] == "-", "-", "mismatch"))
    seq_df[,x] <- ifelse(seq_df[,x] == "mismatch" & seq_df[,ref_seq_name] == "-", "insertion", seq_df[,x])
    seq_df[,x] <- ifelse(seq_df[,x] == "-" & seq_df[,ref_seq_name] != "-", "gap", seq_df[,x])
  }
  seq_df[,ref_seq_name] <- ifelse(seq_df[,ref_seq_name] == "-", "gap", seq_df[,ref_seq_name])
  seq_df <- as.data.frame(tidyr::pivot_longer(seq_df, cols = dplyr::all_of(c(names(seq_df)[which(!names(seq_df) %in% c("position", ref_seq_name))], ref_seq_name)), names_to = name_col, values_to = seq_col))


  return(seq_df)
}
