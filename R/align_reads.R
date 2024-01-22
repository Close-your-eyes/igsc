#' Align reads r1 and r2 to a set of reference sequences
#'
#' @param r1 vector of r1 reads from paired sequencing, must be same length as r2
#' @param r2 see r1
#' @param r1_table optional, table of r1 generated with table(r1, useNA = "no")
#' @param r2_table see r1_table
#' @param ref_seq_list names list of vectors of reference seqences to align r1 and r2 to
#' @param mapply_fun
#' @param min_reads_to_plot
#' @param max_reads_to_plot
#' @param maxmis
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
align_reads <- function(r1,
                        r2,
                        r1_table = NULL,
                        r2_table = NULL,
                        ref_seq_list,
                        mapply_fun = parallel::mcmapply,
                        min_reads_to_plot = 3,
                        max_reads_to_plot = 200,
                        maxmis = 0,
                        ...) {

  if (missing(r1) || missing(r2)) {
    stop("r1 or r2 missing.")
  }

  if (length(r1) != length(r2)) {
    stop("r1 and r2 are of different length.")
  }

  # lapply_fun parallel::mclapply, lapply, purrr::map
  # ... args to lapply_fun like mc.cores
  lapply_fun <- match.fun(lapply_fun)


  if (is.null(r1_table)) {
    r1_table <- table(r1, useNA = "no")
  } else {
    if (!all(names(r1_table) %in% r1)) {
      message("some names of r1_table are not in r1.")
    }
    if (!all(r1[which(!is.na(r1))] %in% names(r1_table))) {
      message("some r1 are not in names of r1_table.")
    }
  }
  if (is.null(r2_table)) {
    r2_table <- table(r2, useNA = "no")
  } else {
    if (!all(names(r2_table) %in% r2)) {
      message("some names of r2_table are not in r2.")
    }
    if (!all(r1[which(!is.na(r1))] %in% names(r1_table))) {
      message("some r1 are not in names of r1_table.")
    }
  }

  if (is.null(names(ref_seq_list))) {
    stop("ref_seq_list does not have names.")
  }

  if (!is.numeric(min_reads_to_plot) || !is.numeric(max_reads_to_plot)) {
    stop("min_reads_to_plot and max_reads_to_plot have to be numeric.")
  }
  if (min_reads_to_plot < 0 || max_reads_to_plot < 0) {
    stop("min_reads_to_plot and max_reads_to_plot have be positive integers.")
  }
  if (max_reads_to_plot <= min_reads_to_plot) {
    stop("max_reads_to_plot should be greate than min_reads_to_plot")
  }

  if (!is.numeric(maxmis)) {
    stop("maxmis has to be numeric.")
  }
  if (maxmis < 0) {
    stop("maxmis has to be zero or a positive integer.")
  }

  # rrr <- r1
  # revcomp <- F
  # names of first argument become names of returned list
  match_df_list <- purrr::map2(list(r1 = r1, r2 = r2), c(F,T), function(rrr, revcomp, ...) {
    # here, do not make r2 revcomp but the subjects
    message("aligning reads.")
    match_df <- match_read_unpaired3(reads = rrr,
                                     ref_seq_list = ref_seq_list,
                                     revcomp_subject = revcomp,
                                     maxmis = maxmis,
                                     mapply_fun = mapply_fun, # purrr::map2 will not work due to argument passing, .x and .y, can this be fixed?
                                     ...) # ...
    message("  done.")
    if (is.null(match_df)) {
      return(NULL)
    }
    ## nesting needed?? - maybe not, if so, reduce complexity
    #match_df_nest <- purrr::map(match_df, tidyr::nest, .by = ref_seq_ind, .key = "read_inds")
    # return freq of unique reads?
    return(match_df)
  })

  # remove NULL from lists
  # TODO ??

  # add read sequences to data frame
  for (i in c("r1", "r2")) {
    # treat string as variable name: https://stackoverflow.com/questions/9057006/getting-strings-recognized-as-variable-names-in-r
    # eval(as.name(paste()))
    match_df_list[[i]] <- purrr::map(match_df_list[[i]], function(x) {
      x[,paste0("read_seq_", i)] <- eval(as.name(paste(i)))[x$read_ind]
      return(x)
    })
  }

  # get paired reads which aligned
  ref_names <- unique(c(names(match_df_list[["r1"]]), names(match_df_list[["r2"]])))
  ref_names <- setNames(ref_names, ref_names)
  r1_r2_intersect <- purrr::map(ref_names, function(ref_name) {
    if (!is.null(match_df_list[["r1"]][[ref_name]]) && !is.null(match_df_list[["r2"]][[ref_name]])) {
      paired <- intersect(match_df_list[["r1"]][[ref_name]]$read_ind, match_df_list[["r2"]][[ref_name]]$read_ind)
    } else {
      paired <- NULL
    }
    return(list(paired = paired,
                r1_unpaired = setdiff(match_df_list[["r1"]][[ref_name]]$read_ind, paired),
                r2_unpaired = setdiff(match_df_list[["r2"]][[ref_name]]$read_ind, paired)))
  })

  match_df_list[["r1_r2"]] <- purrr::map(ref_names, function(ref_name) {
    if (!is.null(match_df_list[["r1"]][[ref_name]]) && !is.null(match_df_list[["r2"]][[ref_name]])) {
      dplyr::full_join(match_df_list[["r1"]][[ref_name]],
                       match_df_list[["r2"]][[ref_name]], by = dplyr::join_by(ref_seq_ind, read_ind))
      # what if r1 and r2 aligned to different ref_seq_inds? - check with a join df of paired reads only; na NA should be there
    } else {
      return(NULL)
    }
  })
  ## now use r1_r2_match_df for igsc::MultiplePairwiseAlignmentsToOneSubject
  ## loop over all ref_seq which got reads aligned
  #ref_seq_ind <- 1485

  plot_list <- purrr::map(names(match_df_list[["r1_r2"]]), function(ref_name) {
    ref_seq_ind_table <- table(match_df_list[["r1_r2"]][[ref_name]][["ref_seq_ind"]])
    # min and max number of paired reads to plot per ref_seq_ind
    inds_to_plot <- intersect(names(which(ref_seq_ind_table >= min_reads_to_plot)), names(which(ref_seq_ind_table <= max_reads_to_plot)))
    plots <- purrr::map(setNames(inds_to_plot, inds_to_plot), function(ref_seq_ind) {
      #print(ref_seq_ind)
      temp <-
        match_df_list[["r1_r2"]][[ref_name]] %>%
        dplyr::filter(ref_seq_ind == !!ref_seq_ind) %>%
        dplyr::group_by(read_seq_r1, read_seq_r2) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop")

      # make r2 rev comp
      temp$read_seq_r2[which(!is.na(temp$read_seq_r2))] <- revcompDNA(temp$read_seq_r2[which(!is.na(temp$read_seq_r2))])
      reads_groups <- apply(temp, 1, function(z) c(z["read_seq_r1"], z["read_seq_r2"]), simplify = F)
      names(reads_groups) <- as.character(seq(1, length(reads_groups), 1))
      ## assign names to reads groups based on if any read is NA
      r1_na <- sapply(reads_groups, function(z) is.na(z[["read_seq_r1"]]))
      r2_na <- sapply(reads_groups, function(z) is.na(z[["read_seq_r2"]]))
      if (length(temp_ind <- which(!r1_na & !r2_na)) > 0) {
        names(reads_groups)[temp_ind] <- paste0("read_pair_", seq(1,length(temp_ind),1), "\n(n=", temp$n[temp_ind], ")")
      }
      if (length(temp_ind <- which(!r1_na & r2_na)) > 0) {
        names(reads_groups)[temp_ind] <- paste0("r1_", seq(1,length(temp_ind),1), "\n(n=", temp$n[temp_ind], ")")
      }
      if (length(temp_ind <- which(r1_na & !r2_na)) > 0) {
        names(reads_groups)[temp_ind] <- paste0("r2_", seq(1,length(temp_ind),1), "\n(n=", temp$n[temp_ind], ")")
      }

      # reorder
      reads_groups <- c(reads_groups[which(grepl("pair", names(reads_groups)))],
                        reads_groups[which(!grepl("pair", names(reads_groups)))])

      plot_data <- igsc::MultiplePairwiseAlignmentsToOneSubject(subject = Biostrings::DNAStringSet(setNames(ref_seq_list[[ref_name]][as.numeric(ref_seq_ind)], "ref_seq")),
                                                                patterns = reads_groups,
                                                                type = "local",
                                                                seq_type = "NT",
                                                                verbose = F)
      plot_data[["plot"]] <-
        plot_data[["plot"]] +
        labs(title = stringr::str_wrap(names(ref_seq_list[[ref_name]][as.numeric(ref_seq_ind)]), width = 40), y = NULL, x = NULL) +
        theme(title = element_text(size = 6))
      return(plot_data)
    })
  })

  # stats
  # attach summary statistic:
  # total reads
  # # if differentiate between pre- and post-QC, then read_paired_reads needs to be within the align function
  # freq of unique reads
  # number of total matching reads
  # number of unique matching
  # # the same for paired r1,r2 and unpaired

  r1_r2_unqiue_paired_n <- nrow(tidyr::drop_na(data.frame(r1 = r1, r2 = r2))) # length(which(!is.na(r1) & !is.na(r2)))
  r1_r2_unqiue_paired_distinct_n <- nrow(dplyr::distinct(tidyr::drop_na(data.frame(r1 = r1, r2 = r2))))
  ## NULL?
  r1_r2_match_df_paired <- purrr::map(match_df_list[["r1_r2"]], tidyr::drop_na)

  stat_df_wide <- purrr::map(ref_names, function(ref_name) {
    data.frame(r1__matches_total_n = ifelse(!is.null(match_df_list[["r1"]][[ref_name]]), nrow(match_df_list[["r1"]][[ref_name]]), NA),
               r1__matches_total_freq = ifelse(!is.null(match_df_list[["r1"]][[ref_name]]), nrow(match_df_list[["r1"]][[ref_name]])/length(r1), NA),
               r1__matches_unique_n = ifelse(!is.null(match_df_list[["r1"]][[ref_name]]), nrow(match_df_list[["r1"]][[ref_name]]), NA),
               r1__matches_unique_freq = ifelse(!is.null(match_df_list[["r1"]][[ref_name]]), nrow(match_df_list[["r1"]][[ref_name]])/length(r1_table), NA),
               r1__ref_targets_n = length(unique(match_df_list[["r1"]][[ref_name]][["ref_seq_ind"]])),

               r2__matches_total_n = ifelse(!is.null(match_df_list[["r2"]][[ref_name]]), nrow(match_df_list[["r2"]][[ref_name]]), NA),
               r2__matches_total_freq =ifelse(!is.null(match_df_list[["r2"]][[ref_name]]),  nrow(match_df_list[["r2"]][[ref_name]])/length(r2), NA),
               r2__matches_unique_n = ifelse(!is.null(match_df_list[["r2"]][[ref_name]]), nrow(match_df_list[["r2"]][[ref_name]]), NA),
               r2__matches_unique_freq = ifelse(!is.null(match_df_list[["r2"]][[ref_name]]), nrow(match_df_list[["r2"]][[ref_name]])/length(r2_table), NA),
               r2__ref_targets_n = length(unique(match_df_list[["r2"]][[ref_name]][["ref_seq_ind"]])),

               paired__matches_total_n = ifelse(!is.null(r1_r2_match_df_paired[[ref_name]]), nrow(r1_r2_match_df_paired[[ref_name]]), NA),
               paired__matches_total_freq = ifelse(!is.null(r1_r2_match_df_paired[[ref_name]]), nrow(r1_r2_match_df_paired[[ref_name]])/length(which(!is.na(r1) & !is.na(r2))), NA),
               paired__matches_unique_n = ifelse(!is.null(r1_r2_match_df_paired[[ref_name]]), nrow(dplyr::distinct(r1_r2_match_df_paired[[ref_name]])), NA),
               paired__matches_unique_freq = ifelse(!is.null(r1_r2_match_df_paired[[ref_name]]), nrow(dplyr::distinct(r1_r2_match_df_paired[[ref_name]]))/r1_r2_unqiue_paired_distinct_n, NA),
               paired__ref_targets_n = length(unique(r1_r2_match_df_paired[[ref_name]][["ref_seq_ind"]])))
  })
  r1_ref_seq_inds_matches <- purrr::map(ref_names, function(ref_name) table(match_df_list[["r1"]][[ref_name]][["ref_seq_ind"]]))
  r2_ref_seq_inds_matches <- purrr::map(ref_names, function(ref_name) table(match_df_list[["r2"]][[ref_name]][["ref_seq_ind"]]))
  paired_ref_seq_inds_matches <- purrr::map(ref_names, function(ref_name) table(r1_r2_match_df_paired[[ref_name]][["ref_seq_ind"]]))

  stat_df_long <- purrr::map(stat_df_wide, function(x) x %>% tidyr::pivot_longer(cols = names(.), names_to = "stat"))
  stat_df <- purrr::map(stat_df_long, function(x) x %>% tidyr::separate(col = stat, into = c("read", "stat"), sep = "__") %>% tidyr::pivot_wider(names_from = read, values_from = value))

  # what to return now?
  read_matches <- list(plot_data = plot_list,
                       match_df_nest = match_df_list,
                       read_indices = r1_r2_intersect,
                       match_stats = list(wide = stat_df_wide,
                                          long = stat_df_long,
                                          norm = stat_df),
                       ref_seq_ind_matches = list(r1 = r1_ref_seq_inds_matches,
                                                  r2 = r2_ref_seq_inds_matches,
                                                  paired = paired_ref_seq_inds_matches))

  read_stats_wide <- data.frame(r1__total_n = length(r1),
                                r1__unique_n = length(r1_table),
                                r1__unique_freq = length(r1_table)/length(r1),
                                r2__total_n = length(r2),
                                r2__unique_n = length(r2_table),
                                r2__unique_freq = length(r2_table)/length(r2),
                                paired__total_n = r1_r2_unqiue_paired_n,
                                paired__unique_n = r1_r2_unqiue_paired_distinct_n,
                                paired__unique_freq = r1_r2_unqiue_paired_distinct_n/r1_r2_unqiue_paired_n)

  read_stats_long <-
    read_stats_wide %>%
    tidyr::pivot_longer(cols = names(.), names_to = "stat")
  read_stats <-
    read_stats_long %>%
    tidyr::separate(col = stat, into = c("read", "stat"), sep = "__") %>%
    tidyr::pivot_wider(names_from = read, values_from = value)

  return(list(read_stats = list(wide = read_stats_wide,
                                long = read_stats_long,
                                norm = read_stats),
              read_matches = read_matches))

}


revcompDNA <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}
