#' Title
#'
#' @param df
#' @param ref
#' @param pos_col
#' @param seq_col
#' @param name_col
#' @param match_symbol
#' @param mismatch_symbol
#' @param keep_gaps
#' @param nonref_mismatch_as
#'
#' @returns
#' @export
#'
#' @examples
compare_nonref_r <- function(df,
                             ref,
                             pos_col = "position",
                             seq_col = "seq",
                             name_col = "seq.name",
                             match_symbol = ".",
                             mismatch_symbol = "x",
                             keep_gaps = T,
                             nonref_mismatch_as = "base") {

  for (i in setdiff(unique(df[[name_col]]), ref)) {
    idx <- which(df[[name_col]] == i)
    for (j in df[idx,pos_col]) {
      v1 <- df[intersect(which(df[[pos_col]] == j), which(df[[name_col]] == i)),seq_col]
      v2 <- df[intersect(which(df[[pos_col]] == j), which(df[[name_col]] == ref)),seq_col]

      if (v1 == v2) {
        ## match branch
        if (v1 != "-" || !keep_gaps) {
          df[intersect(which(df[[pos_col]] == j), which(df[[name_col]] == i)),seq_col] <- match_symbol
        }
      } else {
        ## mismatch branch
        # v1==- and v2!=-
        # v1=!- and v2==- or whatever
        if (nonref_mismatch_as == "mismatch_symbol") {
          if (v1 != "-" || !keep_gaps) { # does make sense I think
            df[intersect(which(df[[pos_col]] == j), which(df[[name_col]] == i)),seq_col] <- mismatch_symbol
          }
        }
      }
    }
  }
  return(df)
}
