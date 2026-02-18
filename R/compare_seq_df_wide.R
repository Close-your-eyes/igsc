#' Compare sequences in columns of a data frame
#'
#' Prepare for plotting with algmnt_plot.
#'
#' @param df data frame
#' @param ref reference or subject to compare non-ref or patterns to
#' (col name); can be guessed if NULL
#' @param non_ref non-ref or pattern col names; can be guessed if NULL
#' @param pos_col column name of position column
#' @param seq_col sequence column name when return_as_long
#' @param name_col sequence name column name when return_as_long
#' @param match_symbol symbol for matching positions
#' @param mismatch_symbol symbol for mismatching positions
#' @param change_nonref alter non-refs to match, mismatch or insertion
#' @param nonref_mismatch_as what should mismatches become in non-ref
#' @param change_ref alter ref to match, mismatch or insertion
#' @param ref_mismatch_as what should mismatches become in ref
#' @param insertion_as what should insertion become
#' @param keep_match_gaps keep matching gaps as gaps (T)? of match (F)?
#' @param return_as_long return as long data frame?
#' @param rm_pure_NA_non_ref if all non-refs are NA at one position, remove this
#' position? intended to make subsequent plot tighter.
#'
#' @returns data frame
#' @export
#'
#' @examples
compare_seq_df_wide <- function(df,
                                ref = NULL,
                                non_ref = NULL,
                                pos_col = "position",
                                seq_col = "seq",
                                name_col = "seq.name",
                                #seq_original = "seq_original",
                                match_symbol = ".",
                                mismatch_symbol = "x",
                                change_nonref = F,
                                nonref_mismatch_as = c("base", "mismatch_symbol"),
                                change_ref = F,
                                ref_mismatch_as = c("base", "mismatch_symbol"),
                                insertion_as = "base",
                                keep_match_gaps = T,
                                return_as_long = F,
                                rm_pure_NA_non_ref = F) {


  nonref_mismatch_as <- rlang::arg_match(nonref_mismatch_as)
  ref_mismatch_as <- rlang::arg_match(ref_mismatch_as)

  if (!pos_col %in% names(df)) {
    stop(pos_col, " not found in df.")
  }

  if (is.null(ref)) {
    ref <- get_ref_as_longest_seq(
      df = df,
      shape = "wide"
    )
  }

  name_order <- names(which(purrr::map_lgl(df, is.character)))
  if (is.null(non_ref)) {
    non_ref <- name_order[which(name_order != ref)]
  } else {
    non_ref <- intersect(name_order, non_ref)
  }
  name_order <- name_order[which(name_order %in% c(ref, non_ref))]


  if (anyDuplicated(df[[pos_col]])) {
    message("Found duplicate rows in df. Will use dplyr::distict to fix this but you should check your data frame.")
    df <- dplyr::distinct(df, !!rlang::sym(pos_col), .keep_all = T)
  }
  #
  # df <- df %>%
  #   dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
  #   tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col))

  # matches to pattern
  if (change_nonref) {
    # super quick version of rowwise pipeline below
    for (i in non_ref) {
      df[[i]] <- igsc:::mutate_value_cpp(
        value               = df[[i]],
        ref_col             = df[[ref]],
        match_symbol        = match_symbol,
        mismatch_symbol     = mismatch_symbol,
        keep_match_gaps     = keep_match_gaps,
        nonref_mismatch_as  = nonref_mismatch_as
      )
      df[[i]][which(df[[i]] == "NA")] <- NA
      if (insertion_as != "base") {
        df[[i]] <- ifelse(df[[ref]] == "-" & df[[i]] != "-", insertion_as, df[[i]])
      }
    }
  }

  # matches to ref
  if (change_ref) {
    if (ref_mismatch_as == "base") {
      mismatch_replace <- df[[ref]]
    } else if (ref_mismatch_as == "mismatch_symbol") {
      mismatch_replace <- rep(mismatch_symbol, nrow(df))
    }

    if (change_nonref) {
      fun1 <- function(x) length(unique(x[intersect(which(!is.na(x)), which(x != match_symbol))]))
    } else {
      fun1 <- function(x) length(unique(x[which(!is.na(x))]))
    }
    fun2 <- function(x) all(x == "-")
    df[[ref]] <- ifelse(
      apply(
        df[,c(non_ref)],
        MARGIN = 1,
        FUN = fun1
      ) == 0,
      # keep "-" when gaps in all seq - keep_match_gaps decides, tested a bit
      ifelse(
        apply(
          df[,c(ref, non_ref)],
          MARGIN = 1,
          FUN = fun2
        ) & keep_match_gaps, "-", match_symbol),
      mismatch_replace
    )
    if (insertion_as != "base") {
      # any or all? any(x[-1] == "-")
      df[[ref]] <- ifelse(apply(df[,c(ref, non_ref)], 1,
                                function(x) (!x[1] %in% c("-", match_symbol, mismatch_symbol) && any(x[-1][which(!is.na(x[-1]))] == "-"))), insertion_as, df[[ref]])
    }
  }

  if (rm_pure_NA_non_ref) {
    na_sum <- apply(df[,non_ref], 1, function(x) sum(is.na(x)))
    df <- df[which(!dplyr::near(na_sum/length(non_ref), 1)),]
    df <- dplyr::arrange(df, !!rlang::sym(pos_col))
    df[[pos_col]] <- seq(1, nrow(df), 1)
  }

  if (return_as_long) {
    df <- tidyr::pivot_longer(df,cols = dplyr::all_of(c(ref, non_ref)),
                              names_to = name_col,
                              values_to = seq_col)
  }

  # if (!is.null(seq_original)) {
  #   df2 <- dplyr::left_join(df2,
  #                           dplyr::rename(df, {{seq_original}} := !!rlang::sym(seq_col)),
  #                           by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  # }
  # df2 <- dplyr::mutate(df2, {{name_col}} := factor(!!rlang::sym(name_col), levels = name_order))
  #
  # join other columns from initial data frame
  # if ("pattern" %in% names(df2)) {
  #   ## try brathering::get_unique_level_columns instead of coalesce_other_cols_with_unique_vals
  #   ## when is pattern there?
  #   df2 <- df2 %>%
  #     brathering::coalesce_join(coalesce_other_cols_with_unique_vals(df = df,
  #                                                                    name_col = name_col,
  #                                                                    pos_col = pos_col,
  #                                                                    seq_col = seq_col), by = "pattern")
  # }
  #
  # if (any(!names(df)[-which(names(df) == seq_col)] %in% names(df2))) {
  #   names(df)[!names(df) %in% names(df2)]
  #   df2 <- dplyr::left_join(df2, dplyr::select(df, -!!rlang::sym(seq_col)), by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  # }

  attr(df, "subject_name") <- ref
  attr(df, "ref") <- ref

  return(df)

}
