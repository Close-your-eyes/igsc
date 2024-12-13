#' Title
#'
#' @param df_long
#' @param ref
#' @param pos_col
#' @param seq_col
#' @param name_col
#' @param seq_original
#' @param match_symbol
#' @param mismatch_symbol
#' @param change_pattern
#' @param pattern_mismatch_as
#' @param change_ref
#' @param ref_mismatch_as
#' @param insertion_as
#'
#' @return
#' @export
#'
#' @examples
compare_seq_df_long <- function(df_long,
                                ref,
                                pos_col = "position",
                                seq_col = "seq",
                                name_col = "seq.name",
                                seq_original = "seq_original",
                                match_symbol = ".",
                                mismatch_symbol = "x",
                                change_pattern = F,
                                pattern_mismatch_as = c("base", "mismatch_symbol"),
                                change_ref = F,
                                ref_mismatch_as = c("base", "mismatch_symbol"),
                                insertion_as = c("base", "insertion")) {

  # insertion_as can be base or any other character, like "x"

  if (!change_pattern && !change_ref) {
    if (!is.null(seq_original)) {
      if (!seq_col %in% names(df_long)) {
        stop(seq_col, " not found in df_long.")
      }
      df_long[[seq_original]] <- df_long[[seq_col]]
    }
    return(df_long)
  }

  if (!seq_col %in% names(df_long)) {
    stop(seq_col, " not found in df_long.")
  }
  if (!pos_col %in% names(df_long)) {
    stop(pos_col, " not found in df_long.")
  }
  if (!name_col %in% names(df_long)) {
    stop(name_col, " not found in df_long.")
  }
  if (missing(ref)) {
    stop("ref is missing.")
  }
  if (!ref %in% df_long[,name_col,drop=T]) {
    stop(ref, " not in ",  name_col, " column of df_long.")
  }

  pattern_mismatch_as <- match.arg(pattern_mismatch_as, c("base", "mismatch_symbol"))
  ref_mismatch_as <- match.arg(ref_mismatch_as, c("base", "mismatch_symbol"))
  insertion_as <- match.arg(insertion_as, c("base", "insertion"))

  name_order <- levels(df_long[,name_col,drop=T])
  if (is.null(name_order)) {
    name_order <- unique(df_long[,name_col,drop=T])
  }
  non_ref <- name_order[which(name_order != ref)]

  check_dups <-
    df_long %>%
    dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
    dplyr::summarise(n = dplyr::n(), .by = c(!!rlang::sym(pos_col), !!rlang::sym(name_col))) %>%
    dplyr::filter(n > 1L)
  if (nrow(check_dups) > 1) {
    message("Found duplicate rows in df_long. Will use dplyr::distict to fix this but you should check your data frame.")
    df_long <- dplyr::distinct(df_long)
  }

  df2 <-
    df_long %>%
    dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col))

  # matches to pattern
  if (change_pattern) {
    df2 <- df2 %>% tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(!names(.) %in% c(ref, pos_col))]))
    if (pattern_mismatch_as == "base") {
      df2 <- dplyr::mutate(df2, value = ifelse(value == !!rlang::sym(ref), match_symbol, value))
    } else if (pattern_mismatch_as == "mismatch_symbol") {
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
    if (change_pattern) {
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) length(unique(x[intersect(which(!is.na(x)), which(x != match_symbol))]))) == 1, match_symbol, mismatch_replace)
    } else {
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) length(unique(x[which(!is.na(x))]))) == 1, match_symbol, mismatch_replace)
    }
    if (insertion_as != "base") {
      # any or all? any(x[-1] == "-")
      df2[[ref]] <- ifelse(apply(df2[,c(ref, non_ref)], 1, function(x) (!x[1] %in% c("-", match_symbol, mismatch_symbol) && any(x[-1][which(!is.na(x[-1]))] == "-"))), insertion_as, df2[[ref]])
    }
  }

  df2 <- df2 %>% tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(names(.) != pos_col)]), names_to = name_col, values_to = seq_col)
  if (!is.null(seq_original)) {
    df2 <- dplyr::left_join(df2, df_long %>% dplyr::rename({{seq_original}} := !!rlang::sym(seq_col)), by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  }
  df2 <- dplyr::mutate(df2, {{name_col}} := factor(!!rlang::sym(name_col), levels = name_order))
  # join other columns from initial data frame
  df2 <-
    df2 %>%
    scexpr:::coalesce_join(coalesce_other_cols_with_unique_vals(df_long = df_long,
                                                                name_col = name_col,
                                                                pos_col = pos_col,
                                                                seq_col = seq_col), by = "pattern")

  if (any(!names(df_long)[-which(names(df_long) == seq_col)] %in% names(df2))) {
    names(df_long)[!names(df_long) %in% names(df2)]
    df2 <- dplyr::left_join(df2, dplyr::select(df_long, -!!rlang::sym(seq_col)), by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  }

  return(df2)
}

coalesce_other_cols_with_unique_vals <- function(df_long,
                                                 name_col,
                                                 pos_col,
                                                 seq_col) {
  unique_add_cols <-
    df_long %>%
    dplyr::group_by(!!rlang::sym(name_col)) %>%
    dplyr::summarise(dplyr::across(-c(!!rlang::sym(pos_col), !!rlang::sym(seq_col)), dplyr::n_distinct)) %>%
    tidyr::pivot_longer(cols = -pattern) %>%
    dplyr::filter(value == 1)
  add_cols_vals <-
    df_long %>%
    dplyr::select(pattern, unique(unique_add_cols$name)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(cols = - pattern)
  add_cols_vals_unique_wide <-
    unique_add_cols %>%
    dplyr::select(-value) %>%
    left_join(add_cols_vals, by = dplyr::join_by(pattern, name)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = name, values_from = value)
  return(add_cols_vals_unique_wide)
}

## other procedure:
## assign matches and mismatches
## similar procedure as in igsc::MultiplePairwiseAlignmentsToOneSubject
'if (any(unlist(matches_to_origin_and_feature))) {
  match_mismatch_list <- lapply(stats::setNames(names(df0)[-c(1,2)], names(df0)[-c(1,2)]), function(x) df0[,x] == df0[,"origin"]) # subject seq in df (not df.original) may already contain "insertion"; for overlapping pattern where one receives an insertion and the other not, this is relevant
  if (any(sapply(matches_to_origin_and_feature, "[", 1))) {
    # do it outside of loop below to avoid double calculation
    matches <- purrr::pmap_lgl(match_mismatch_list, function(...) any_false(unlist(list(...))))
  }
  df0 <- purrr::map(setNames(matches_to_origin_and_feature, sapply(matches_to_origin_and_feature, function(x) paste(ifelse(x, "match", "base"), collapse = "_"))), function(y) {
    if (any(y)) {
      df0_match <- df0
      if (y[2]) {
        # matches/mismatches to patterns
        for (i in names(match_mismatch_list)) {
          df0_match[,i] <- ifelse(match_mismatch_list[[i]], "match", "mismatch")
        }
      }
      if (y[1]) {
        # matches/mismatches to subject
        df0_match[which(matches),"origin"] <- "match"
        df0_match[which(!matches),"origin"] <- "mismatch"
      }
      return(df0_match)
    } else {
      return(df0)
    }
  })
} else {
  df0 <- list("base_base" = df0)
}'
