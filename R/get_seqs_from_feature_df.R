#' Get sequences based on data generated with webscrape_ncbi
#'
#' @param feature_df data frame of features as returned by webscrape_ncbi
#' @param origin reference sequence for features
#' @param concat concatenate sequences which have gaps in origin
#' @param make_revcomp return the reverse complement of sequences when complement column is TRUE,
#' for alignment to origin do not have revcomp return,
#' for writing features into a ref genome do have the it return
#' @param min_to_max return whole sequences from start to end including gaps (TRUE), or remove gaps (FALSE)
#' @param return which formats to return; restriction to what is needed saves memory and increases speed
#' @param matches_to_origin_and_feature mark matches and mismatches in origin and/or feature sequences
#' @param order_features order features by min position
#' @param compare_seq_df_long_args
#'
#' @return
#' @export
#'
#' @examples
get_seqs_from_feature_df <- function(feature_df,
                                     origin,
                                     concat = T,
                                     make_revcomp = F,
                                     min_to_max = F,
                                     return = c("sequences", "df_wide", "df_long"),
                                     matches_to_origin_and_feature = list(c(F,T),c(F,F)),
                                     order_features = F,
                                     compare_seq_df_long_args = list()) {

  return <- match.arg(return, c("sequences", "df_wide", "df_long"), several.ok = T)

  if (missing(origin)) {
    stop("origin sequence has to be provided.")
  }

  # how to handle those symbols from NCBI actually?
  if (any(grepl("<|>", feature_df$range))) {
    message("'<' or '>' found in range column. Those will be removed.")
  }

  # prepare boundaries
  boundaries <- lapply(stats::setNames(feature_df$range, make.unique(feature_df$value)), function(x) sapply(strsplit(strsplit(gsub("<|>", "", x), ",")[[1]], "\\.\\."), as.numeric, simplify = F))
  if (order_features) {
    temp_order <- order(sapply(boundaries, function(x) x[[1]][1])) # min(unlist(x)) # min position vs. first position (first exon in circular genome may be at a high position)
    boundaries <- boundaries[temp_order]
    feature_df <- feature_df[temp_order,]
  }

  if (min_to_max) {
    boundaries <- lapply(boundaries, function(x) list(c(x[[1]][1], unlist(x)[length(unlist(x))]))) # list(c(min(unlist(x)), max(unlist(x)))))
    for (i in seq_along(boundaries)) {
      feature_df$range[i] <- paste0(boundaries[[i]][[1]][1], "..", boundaries[[i]][[1]][2])
    }
  }

  sequences <- NULL
  if ("sequences" %in% return) {
    # concat the sequence from segments
    # browser() with "NC_006273.2", check index 878. length of x is 1 - how to handle?
    sequences <- purrr::pmap(list(x = boundaries[878], revcomp = feature_df$complement[878],
                                  value = feature_df$value[878], range = feature_df$range[878]), function(x,revcomp,value,range) {
                                    seq <- unlist(lapply(x, function(y) {
                                      if (y[1] > y[2]) {
                                        # exon 1 is at later position as first one
                                        return(paste0(substr(origin, y[1], nchar(origin)), substr(origin, 1, y[2])))
                                      } else {
                                        # exon 1 at smallest position
                                        return(substr(origin, y[1], y[2]))
                                      }
                                    }))
                                    if (revcomp && make_revcomp) {
                                      seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))
                                      seq <- rev(seq) # if no pasting above, separate seq have to be reversed here
                                    }
                                    if (concat) {
                                      seq <- paste(seq, collapse = "")
                                    } else {
                                      # w/o concat, name each subseq individually
                                      if (revcomp && make_revcomp) {
                                        # check ?!
                                        names(seq) <- paste0(value, "___", rev(strsplit(range, ",")[[1]]))
                                      } else {
                                        names(seq) <- paste0(value, "___", strsplit(range, ",")[[1]])
                                      }
                                    }
                                    return(seq)
                                  })

    if (anyDuplicated(feature_df$value) != 0 && "Feature" %in% names(feature_df)) {
      names(sequences) <- make.unique(paste0(feature_df$value, "___", feature_df$Feature))
      # range could be added in addition or subfeature
    } else {
      names(sequences) <- make.unique(feature_df$value)
    }
  }

  dfs <- purrr::pmap(list(x = boundaries, revcomp = feature_df$complement,
                          value = feature_df$value, range = feature_df$range), function(x,revcomp,value,range) {
                            seq <- unlist(lapply(x, function(y) substr(origin, y[1], y[2])))
                            if (revcomp && make_revcomp) {
                              seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))
                              seq <- rev(seq) # if no pasting above, separate seq have to be reversed here
                            }
                            seq <- paste(seq, collapse = "")
                            df <- data.frame(position = unlist(seq2(sapply(x, "[", 1), sapply(x, "[", 2)), use.names = F), seq = strsplit(seq, "")[[1]])
                            names(df)[2] <- value
                            return(df)
                          })
  df0 <- data.frame(position = seq(1, nchar(origin)), origin = strsplit(origin, "")[[1]]) # origin column
  dfs <- join_chunkwise(data_frame_list = dfs,
                        join_fun = dplyr::full_join,
                        join_by = "position",
                        max_final_size = 5)
  df0 <- purrr::reduce(c(list(df0), dfs), dplyr::left_join, by = "position")

  if ("change_pattern" %in% names(compare_seq_df_long_args) || "change_ref" %in% names(compare_seq_df_long_args)) {
    df0_long <- tidyr::pivot_longer(df0, -position, names_to = "seq.name", values_to = "seq")
    df0_long <- do.call(compare_seq_df_long, args = c(compare_seq_df_long_args, list(df_long = df0_long, ref = "origin", seq_original = NULL)))
    df0 <- tidyr::pivot_wider(df0_long, values_from = seq, names_from = seq.name)
  }
  df0 <- list(df0)

  df0_long <- NULL
  if ("df_long" %in% return) {
    df0_long <- purrr::map(df0, function(x) tidyr::pivot_longer(x, cols = -position, names_to = "seq.name", values_to = "seq"))
    # attach a column defining start and end position
    start_end_boundaries_df <- stack(lapply(boundaries, function(x) c(x[[1]][1], unlist(x)[length(unlist(x))])))
    names(start_end_boundaries_df) <- c("position", "seq.name")
    start_end_boundaries_df$seq.name <- as.character(start_end_boundaries_df$seq.name)
    start_end_boundaries_df <- rbind(start_end_boundaries_df, data.frame(seq.name = c("origin", "origin"), position = c(1,nchar(unname(origin)))))
    start_end_boundaries_df$start_end <- c("start","end") # recycling, every second position is start or end, respectively
    df0_long <- purrr::map(df0_long, function(x) dplyr::left_join(x, start_end_boundaries_df, by = c("position" = "position", "seq.name" = "seq.name")))
    if (order_features) {
      df0_long <- purrr::map(df0_long, function(x) {
        x$seq.name <- factor(x$seq.name, levels = names(df0[[1]])[-which(names(df0[[1]]) == "position")]) # do it like this, with names of df0[[1]], in case duplicate column names have been altered; order should have been maintained as defined at the beginning
        return(x)
      })
    } else {
      df0_long <- purrr::map(df0_long, function(x) {
        x$seq.name <- as.factor(x$seq.name)
        return(x)
      })
    }
  }

  if (!"df_wide" %in% return) {
    df0 <- NULL
  }

  # df0 and df0_long being list is remaining from previous procedure
  return(list(sequences = sequences, df_wide = df0[[1]], df_long = df0_long[[1]]))
}

any_false <- function(x) {
  if (all(is.na(x))) {
    return(T)
  } else if (any(!x[which(!is.na(x))])) {
    return(F)
  } else if (all(x[which(!is.na(x))])) {
    return(T)
  } else {
    stop("Logical error.")
  }
}

