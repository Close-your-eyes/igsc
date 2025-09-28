#' Compare sequences to reference in a long data frame format
#'
#' @param df data frame in long format
#' @param ref reference name (in name_col), e.g. the alignment subject
#' @param pos_col name of position column
#' @param seq_col name of sequence column
#' @param name_col name of sequence name column
#' @param seq_original keep original seq_col in this column to be appended
#' @param match_symbol symbol for matching positions
#' @param mismatch_symbol symbol for mismatching positions
#' @param change_nonref change non reference sequences
#' @param nonref_mismatch_as replacement for nonref mismatches
#' @param change_ref change reference sequences
#' @param ref_mismatch_as replacement for ref mismatches
#' @param insertion_as replacement for insertions: base or any symbol
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' df <- data.frame(
#'   stringsAsFactors = FALSE,
#'   sub = rep(c("A","C","T","G","A","C","T","C"),
#'             times = c(4,4,4,4,2,2,2,2)),
#'   position = 1:24,
#'   subject.position = 1:24,
#'   pat5 = c(rep("A",4), rep(NA,20)),
#'   pat2 = c(rep(NA,4), rep("C",4), rep(NA,16)),
#'   pat3 = c(rep(NA,8), rep("T",4), rep(NA,12)),
#'   pat4 = c(rep(NA,12), rep("G",4), rep(NA,8)),
#'   pat1 = c(rep(NA,20), "T","T","C","C")
#' )
#' df = tidyr::pivot_longer(df,
#'                          cols = c(sub, dplyr::starts_with("pat")),
#'                          names_to = "seq.name",
#'                          values_to = "seq")
#' df2 <- compare_seq_df_long(df,
#'                            change_ref = T)
#' algnmt_plot(df2)
#' df3 <- compare_seq_df_long(df,
#'                            change_ref = T,
#'                            change_nonref = T)
#' algnmt_plot(df3)
#' algnmt_plot(df)
#' df$seq.name <- factor(df$seq.name, levels = c("sub", paste0("pat", 1:5)))
#' algnmt_plot(df)
compare_seq_df_long <- function(df,
                                ref = NULL,
                                pos_col = "position",
                                seq_col = "seq",
                                name_col = "seq.name",
                                seq_original = "seq_original",
                                match_symbol = ".",
                                mismatch_symbol = "x",
                                change_nonref = F,
                                nonref_mismatch_as = c("base", "mismatch_symbol"),
                                change_ref = F,
                                ref_mismatch_as = c("base", "mismatch_symbol"),
                                insertion_as = "base") {

  # insertion_as can be base or any other character, like "x"

  if (!change_nonref && !change_ref) {
    if (!is.null(seq_original)) {
      if (!seq_col %in% names(df)) {
        stop(seq_col, " not found in df.")
      }
      df[[seq_original]] <- df[[seq_col]]
    }
    return(df)
  }

  if (!seq_col %in% names(df)) {
    stop(seq_col, " not found in df.")
  }
  if (!pos_col %in% names(df)) {
    stop(pos_col, " not found in df.")
  }
  if (!name_col %in% names(df)) {
    stop(name_col, " not found in df.")
  }
  if (is.null(ref)) {
    temp <- df |>
      dplyr::group_by(!!rlang::sym(name_col)) |>
      dplyr::summarise(n_non_na = sum(!is.na(!!rlang::sym(seq_col)))) |>
      dplyr::slice_max(n_non_na, n = 1)
    if (nrow(temp) == 1) {
      ref <- as.character(temp[[name_col]])
      message("ref: ", ref)
    } else {
      stop("ref could not be guessed.")
    }
  }
  if (!ref %in% df[[name_col]]) {
    stop(ref, " not in ",  name_col, " column of df.")
  }



  nonref_mismatch_as <- rlang::arg_match(nonref_mismatch_as)
  ref_mismatch_as <- rlang::arg_match(ref_mismatch_as)
  #insertion_as <- rlang::arg_match(insertion_as)

  name_order <- levels(df[,name_col,drop=T])
  if (is.null(name_order)) {
    name_order <- unique(df[,name_col,drop=T])
  }
  non_ref <- name_order[which(name_order != ref)]

  check_dups <-
    df %>%
    dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
    dplyr::summarise(n = dplyr::n(), .by = c(!!rlang::sym(pos_col), !!rlang::sym(name_col))) %>%
    dplyr::filter(n > 1L)
  if (nrow(check_dups) > 1) {
    message("Found duplicate rows in df. Will use dplyr::distict to fix this but you should check your data frame.")
    df <- dplyr::distinct(df)
  }

  df2 <- df %>%
    dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col))

  # matches to pattern
  if (change_nonref) {
    df2 <- df2 %>% tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(!names(.) %in% c(ref, pos_col))]))
    if (nonref_mismatch_as == "base") {
      df2 <- dplyr::mutate(df2, value = ifelse(value == !!rlang::sym(ref), match_symbol, value))
    } else if (nonref_mismatch_as == "mismatch_symbol") {
      df2 <- dplyr::mutate(df2, value = ifelse(value == !!rlang::sym(ref), match_symbol, mismatch_symbol))
    }
    if (insertion_as != "base") {
      df2 <- dplyr::mutate(df2, value = ifelse(!!rlang::sym(ref) == "-" & value != "-", insertion_as, value))
    }
    df2 <- tidyr::pivot_wider(df2, names_from = name, values_from = value)
  }

  # matches to ref
  if (change_ref) {
    if (ref_mismatch_as == "base") {
      mismatch_replace <- df2[[ref]]
    } else if (ref_mismatch_as == "mismatch_symbol") {
      mismatch_replace <- rep(mismatch_symbol, nrow(df2))
    }

    if (change_nonref) {
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) length(unique(x[intersect(which(!is.na(x)), which(x != match_symbol))]))) == 1, match_symbol, mismatch_replace)
    } else {
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) length(unique(x[which(!is.na(x))]))) == 1, match_symbol, mismatch_replace)
    }
    if (insertion_as != "base") {
      # any or all? any(x[-1] == "-")
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) (!x[1] %in% c("-", match_symbol, mismatch_symbol) && any(x[-1][which(!is.na(x[-1]))] == "-"))), insertion_as, df2[[ref]])
    }
  }

  df2 <- df2 %>% tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(names(.) != pos_col)]),
                                     names_to = name_col,
                                     values_to = seq_col)

  if (!is.null(seq_original)) {
    df2 <- dplyr::left_join(df2,
                            dplyr::rename(df, {{seq_original}} := !!rlang::sym(seq_col)),
                            by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  }
  df2 <- dplyr::mutate(df2, {{name_col}} := factor(!!rlang::sym(name_col), levels = name_order))

  # join other columns from initial data frame
  if ("pattern" %in% names(df2)) {
    ## try brathering::get_unique_level_columns instead of coalesce_other_cols_with_unique_vals
    ## when is pattern there?
    df2 <- df2 %>%
      brathering::coalesce_join(coalesce_other_cols_with_unique_vals(df = df,
                                                                     name_col = name_col,
                                                                     pos_col = pos_col,
                                                                     seq_col = seq_col), by = "pattern")
  }


  if (any(!names(df)[-which(names(df) == seq_col)] %in% names(df2))) {
    names(df)[!names(df) %in% names(df2)]
    df2 <- dplyr::left_join(df2, dplyr::select(df, -!!rlang::sym(seq_col)), by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  }

  attr(df2, "subject_name") <- ref

  return(df2)
}

coalesce_other_cols_with_unique_vals <- function(df,
                                                 name_col,
                                                 pos_col,
                                                 seq_col) {
  unique_add_cols <-
    df %>%
    dplyr::group_by(!!rlang::sym(name_col)) %>%
    dplyr::summarise(dplyr::across(-c(!!rlang::sym(pos_col), !!rlang::sym(seq_col)), dplyr::n_distinct)) %>%
    tidyr::pivot_longer(cols = -pattern) %>%
    dplyr::filter(value == 1)
  add_cols_vals <-
    df %>%
    dplyr::select(pattern, unique(unique_add_cols$name)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(cols = - pattern)
  add_cols_vals_unique_wide <-
    unique_add_cols %>%
    dplyr::select(-value) %>%
    dplyr::left_join(add_cols_vals, by = dplyr::join_by(pattern, name)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = name, values_from = value)
  return(add_cols_vals_unique_wide)
}
