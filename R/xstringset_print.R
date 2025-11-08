#' Print XStringSet to console
#'
#' XStringSet may contain multiple sequence alignment from DECIPHER::AlignSeqs
#'
#' @param set XStringSet
#' @param ref name of reference sequence
#' @param linewidth print width
#' @param print_pos print position above
#' @param print_pos_end print position of each sequence at lineend (gaps excluded);
#' requires keep_match_gaps incompare_seq_df_long
#' @param pos_by_ref above positions in reference to ref (considering its gaps),
#' or just total position
#' @param col_out colored output?
#' @param call_compare_seq_df_long call igsc::compare_seq_df_long to optionally
#' replace matches and/or mismatches
#' @param compare_seq_df_long_args arguments to compare_seq_df_long
#'
#' @returns xstringset
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
#' granzymes <- c("GZMA","GZMB","GZMH","GZMK","GZMM")
#' out <- get_sequences_from_biomart(granzymes)
#' msaln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(setNames(out$coding,
#'                                                                out$hgnc_symbol)))
#' # print as is
#' xstringset_print(msaln)
#' # print with call_compare_seq_df_long anf GZMB as reference
#' xstringset_print(msaln, call_compare_seq_df_long = T,
#'                  ref = "GZMB", print_pos_end = T)
#' # run compare_seq_df_long  yourself
#' df <- igsc:::xstringset_to_df(msaln)
#' msalndf <- compare_seq_df_long(df, ref = "GZMB",
#'                                seq_original = NULL,
#'                                change_nonref = T)
#' msalnx <- df_to_xstringset(msalndf)
#' # print
#' xstringset_print(msalnx)
xstringset_print <- function(set,
                             ref = NULL,
                             linewidth = 100,
                             print_pos = T,
                             print_pos_end = F,
                             pos_by_ref = T,
                             col_out = T,
                             call_compare_seq_df_long = F,
                             compare_seq_df_long_args = list(seq_original = NULL,
                                                             change_nonref = T,
                                                             change_ref = T,
                                                             keep_match_gaps = T)) {

  if (call_compare_seq_df_long) {
    df <- igsc:::xstringset_to_df(set)
    if ("ref" %in% names(compare_seq_df_long_args)) {
      ref <- compare_seq_df_long_args[["ref"]]
      compare_seq_df_long_args[["ref"]] <- NULL
    }
    type <- guess_type2(df[["seq"]]) # known from igsc:::xstringset_to_df
    df <- Gmisc::fastDoCall(compare_seq_df_long,
                            args = c(list(df = df, ref = ref), compare_seq_df_long_args))

    set <- df_to_xstringset(df, type = type)
    ref <- attr(df, "ref")
  }

  set <- as.character(set)
  if (is.null(names(set))) {
    names(set) <- paste0("seq", stringr::str_pad(seq_along(set), width = 2, pad = "0"))
    ref <- NULL
  }
  if (!is.null(ref) && !ref %in% names(set)) {
    message("ref mot found in names of set.")
    ref <- NULL
  }

  ## order set: ref first
  if (!is.null(ref)) {
    set <- set[c(ref, setdiff(names(set), ref))]
  }

  seq_compared <- ifelse(any(grepl("\\.", set)), T, F)

  # in compare_seq_df_long gaps may be converted to "." if they match the ref
  # in this case they cannot be counted by get_symbols and the sym_line is wrong
  # if such sequence becomes the ref here.
  # that actually requires to call compare_seq_df_long appropriately
  message_printed <- F
  if (!is.null(ref) && pos_by_ref) {
    # when there is ref && pos_by_ref: make sym_line according to it
    if (seq_compared && !call_compare_seq_df_long) {
      message_printed <- T
      message("'.' found in reference seq. when compare_seq_df_long was called, keep_match_gaps = T may be required to allow counting of gaps in reference herein.")
    }
    sym_line <- paste(get_symbols(s = strsplit(set[[ref]], "")[[1]], start.pos = 1), collapse = "")
  } else {
    # just a dummy seq w/o gaps
    sym_line <- paste(get_symbols(s = rep("A", nchar(set[1])), start.pos = 1), collapse = "")
  }

  # make seq names equal in length
  names(set) <- pad_strings(names(set))
  linewidth <- min(c(linewidth, nchar(set[1])))

  # save names which are lost below
  pa_names <- names(set)

  # vectorized string split for every row
  starts <- seq(1, nchar(set[1]), by = linewidth)
  sym_line <- stringr::str_sub_all(sym_line, starts, starts+linewidth-1)[[1]]
  set <- stringr::str_sub_all(set, starts, starts+linewidth-1)

  # end pos after every line, only needed when print_pos_end
  res <- 0
  res <- purrr::map(set, ~res + cumsum(linewidth - stringr::str_count(.x, "-")))

  # add name spacing to symline
  sym_line2 <- purrr::map_chr(sym_line, ~paste0(paste(rep(" ", nchar(pa_names[1])), collapse = ""), .x))

  if (col_out) {
    type <- igsc:::guess_type(unname(unlist(purrr::map(set, strsplit, split = ""))))
    set <- purrr::map(set, ~purrr::map_chr(.x, col_letters, type = type, mc.cores = 1))
  }
  # add name before seq lines
  set <- purrr::map2(set, pa_names, ~paste0(.y, .x, "  "))

  # add position at end
  if (print_pos_end) {
    if (seq_compared && !call_compare_seq_df_long && !message_printed) {
      message("'.' found in reference seq. when compare_seq_df_long was called, keep_match_gaps = T may be required to allow counting of gaps in reference herein.")
    }
    set <- purrr::map2(set, res, ~paste0(.x, .y))
  }

  pos <- 0
  for (i in seq_along(sym_line)) {
    # did not know how to vectorize, so: loop
    c(sym_line[i], pos) %<-% .get_pos_line(name = pa_names[1], chunk = sym_line[i], pos = pos)
  }

  for (i in seq_along(sym_line)) {
    if (print_pos) {
      cat(sym_line[i], "\n")
      cat(sym_line2[i], "\n")
    }
    for (j in seq_along(set)) {
      cat(set[[j]][i], "\n")
    }
    cat(paste(" "), "\n")
  }
  cat(paste(" "), "\n")
  cat(paste(" "), "\n")
  cat(paste(" "), "\n")

  invisible(set)

}
