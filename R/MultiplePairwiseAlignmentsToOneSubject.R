#' Separately align multiple pattern sequences to one subject
#'
#' This function is useful to visualise the alignment position of multiple patterns on one subject.
#' It uses Biostrings::pairwiseAlignment(), obtains the individual alignment boundaries and converts
#' the results to a ggplot object.
#' The function will fail if a gap is induced in the subject and at least two pattern alignments overlap at this gap.
#' Method = local-global avoids indels in subject. Local cuts subject and pattern to the best matching sequence of both.
#'
#' @param subject a named character or named DNAStringSet of one subject (only the DNAStringSet but not DNAString can hold a name)
#' @param patterns a named character vector or named DNAStringSet of patterns to align to the subject sequence
#' @param type the type of alignment passed to Biostrings::pairwiseAlignment; not every type may work well with this function (if there are overlapping ranges of the alignments to the subject for example)
#' @param order_patterns order pattern increasingly by alignment position (start)
#' @param max_mismatch only use patterns that have a maximum number of mismatches
#' with the subject
#' @param fix_subject_indels in case of overlapping indels and shared subject ranges, cut respective patterns to avoid indels
#' @param seq_type set sequence type to AA or NT if necessary; if NULL
#' it is attempted to guess the type
#' @param return_max_mismatch_info_only only return information on mismatches of patterns with the subject;
#' in this case no alignment is calculated
#'
#' @return a list:
#' base.plot ggplot object of alignment shows patterns colored by nt,
#' match.plot ggplot object of alignment shows patterns colored match, mismatch, etc,
#' base.df = df and match.df are the respective data.frames used for plotting,
#' seq min.max.subject.position indicates the outer limits of all aligned patterns (min = start position of first aligned pattern, max = end position of the last aligned pattern)
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' s <- stats::setNames("AAAACCCCTTTTGGGGAACCTTCC", "sub")
#' s <- Biostrings::DNAStringSet(s)
#' p <- stats::setNames(c("TTCC", "CCCC", "TTTT", "GGGG", "AAAA"), c("pat1", "pat2", "pat3", "pat4", "pat5"))
#' p <- Biostrings::DNAStringSet(p)
#' als <- igsc::MultiplePairwiseAlignmentsToOneSubject(subject = s, patterns = p)
#' als_ordered <- igsc::MultiplePairwiseAlignmentsToOneSubject(subject = s, patterns = p, order_patterns = T)
#' }
MultiplePairwiseAlignmentsToOneSubject <- function(subject,
                                                   patterns,
                                                   type = c("global-local", "global", "local", "overlap", "local-global"),
                                                   max_mismatch = NA,
                                                   order_patterns = F,
                                                   fix_subject_indels = F,
                                                   rm_indel_inducing_pattern = F,
                                                   seq_type = NULL,
                                                   return_max_mismatch_info_only = F,
                                                   pairwiseAlignment_args = list(),
                                                   algnmt_plot_args = list(add_length_suffix = T)) {

  # take algorithm from timeline to decided automatically what can be in one row

  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }

  '  if (fix_indels) {
    message("fix_indels does not work yet. Set to F.")
    fix_indels <- F
  }'

  type <- match.arg(type, choices = c("global-local", "global", "local", "overlap", "local-global"))

  # this function assigns new values via assign
  #report identical pattern?
  prep_subject_and_patterns(subject = subject,
                            patterns = patterns,
                            seq_type = seq_type)

  # check for non-DNA characters first
  check_for_invalid_chars(subject = subject,
                          patterns = patterns,
                          max_mismatch = max_mismatch)
  if (return_max_mismatch_info_only) {
    return(list(patterns_invalid = patterns_invalid,
                pattern_mismatching = pattern_mismatching_return))
  }


  # calculate all alignments
  # fastDoCall does not work here, maybe due to method written in C
  pa <- do.call(Biostrings::pairwiseAlignment, args = c(list(subject = subject, pattern = patterns, type = type),
                                                        pairwiseAlignment_args))

  # make pa a list once and then iterate over list entries with purrr/furrr which is quicker!
  # use multiple threads to speed up?!
  # pal <- stats::setNames(as.list(pa), patnames(patterns)terns.names)
  # pal <- stats::setNames(purrr::flatten(parallel::mclapply(split(c(1:length(pa)), ceiling(seq_along(c(1:length(pa)))/10)), function(x) as.list(pa[x]), mc.cores = parallel::detectCores()-1)), names(patterns))

  check_for_indel_induction(pa = pa,
                            patterns = patterns,
                            max_mismatch = max_mismatch,
                            rm_indel_inducing_pattern = rm_indel_inducing_pattern)

  check_for_overlapping_indels(pa = pa,
                               patterns = patterns,
                               subject_inds_indel = subject_inds_indel,
                               fix_subject_indels = fix_subject_indels)

  # subject.ranges are defined herein
  # here also overlapping subject.ranges within groups are checked for
  make_pa_unique_and_order_and_rm_subset_alignments(pa = pa,
                                                    pattern_groups = pattern_groups)

  df <- paste_subject_seq(subject = subject,
                          subject.ranges.unique = subject.ranges.unique,
                          pa.unique = pa.unique)


  # pattern_order is defined here
  paste_patterns_to_subject(subject_indels = subject_indels,
                            patterns = patterns,
                            pa = pa,
                            df = df,
                            pattern_original_order = pattern_original_order)


  df.match <- prep_df_for_algnmt_plot(df = df,
                                      subject_name = names(subject),
                                      pattern_names = names(patterns),
                                      original_names = original_names,
                                      order_patterns = order_patterns,
                                      subject.ranges = subject.ranges,
                                      pattern_order = pattern_order,
                                      pattern_groups = pattern_groups,
                                      matches_to_pattern = T,
                                      matches_to_subject = T)

  df <- prep_df_for_algnmt_plot(df = df,
                                subject_name = names(subject),
                                pattern_names = names(patterns),
                                original_names = original_names,
                                order_patterns = order_patterns,
                                subject.ranges = subject.ranges,
                                pattern_order = pattern_order,
                                pattern_groups = pattern_groups,
                                matches_to_pattern = F,
                                matches_to_subject = F)


  if (!is.null(pattern_groups)) {
    algnmt_plot_args <- c(algnmt_plot_args, list(y_group_col = "pattern.group"))
  }

  # write original names into alignments; when the object cycles through C-code (with altered names) certain symbols (maybe like asterisk (*)) may cause problems.
  pa@pattern@unaligned@ranges@NAMES <- original_names[pa@pattern@unaligned@ranges@NAMES]

  g1 <- Gmisc::fastDoCall(algnmt_plot, args = c(list(algnmt = df,
                                                     algnmt_type = seq_type,
                                                     pairwiseAlignment = pa),
                                                algnmt_plot_args))

  # many things to adjust when y_group_col is provided
  g2 <- Gmisc::fastDoCall(algnmt_plot, args = c(list(algnmt = df.match,
                                                     algnmt_type = seq_type,
                                                     pairwiseAlignment = pa),
                                                algnmt_plot_args))

  return(list(base.plot = g1,
              match.plot = g2,
              algnmt.df.base = df,
              algnmt.df.match = df.match,
              min.max.subject.position = c(df %>% dplyr::filter(seq.name != original_names[1]) %>% dplyr::filter(seq != "-") %>% dplyr::slice_min(order_by = position, n = 1, with_ties = F) %>% dplyr::pull(subject.position),
                                           df %>% dplyr::filter(seq.name != original_names[1]) %>% dplyr::filter(seq != "-") %>% dplyr::slice_min(order_by = -position, n = 1, with_ties = F) %>% dplyr::pull(subject.position)),
              pairwise_alignments = pa,
              #pairwise_alignment_list = pal,
              pattern = if(order_patterns) {patterns[order(purrr::map_int(subject.ranges, min))]} else {patterns},
              pattern_invalid = patterns_invalid,
              pattern_indel_inducing = pattern_indel_inducing,
              pattern_mismatching = pattern_mismatching_return))
}

prep_df_for_algnmt_plot <- function(df,
                                    subject_name,
                                    pattern_names,
                                    original_names,
                                    order_patterns,
                                    subject.ranges,
                                    pattern_order,
                                    pattern_groups,
                                    matches_to_pattern = F,
                                    matches_to_subject = F) {

  # gap
  # match
  # mismatch
  # insertion

  if (matches_to_subject && !matches_to_pattern) {
    stop("For matches_to_subject, matches_to_pattern has to be TRUE.")
  }

  # "-" in subject is a gap
  ## if is.na(subject.position) --> always gap in subject and always insertion in pattern

  if (matches_to_pattern) {
    # is.na(subject.position) --> always gap in subject and always insertion in pattern
    subject_gap <- which(is.na(df[,"subject.position"]))
    for (i in pattern_names) {
      df[subject_gap, i][which(!is.na(df[subject_gap, i]))] <- "insertion"
    }
    df[subject_gap, subject_name] <- "gap"
    # "-" in any pattern --> gap in pattern, insertion in subject
    pattern_gap <- apply(df[,pattern_names], 1, function(x) which(x == "-"))
    pattern_gap_rows <- which(lengths(pattern_gap) > 0)
    df[pattern_gap_rows, subject_name] <- "insertion"
    pattern_gap_cols <- pattern_gap[pattern_gap_rows]
    for (i in seq_along(pattern_gap_rows)) {
      df[pattern_gap_rows[i], names(pattern_gap_cols[[i]])] <- "gap"
    }

    match_mismatch_list <- lapply(stats::setNames(pattern_names, pattern_names), function(x) df[,x] == df[,subject_name])
    for (i in names(match_mismatch_list)) {
      df[which(!df[,i] %in% c("gap", "insertion")),i] <- ifelse(match_mismatch_list[[i]][which(!df[,i] %in% c("gap", "insertion"))], "match", "mismatch")
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

    test3 <- purrr::pmap_lgl(match_mismatch_list, function(...) {
      any_false(unlist(list(...)))
    })
    df[intersect(which(!df[,subject_name] %in% c("gap", "insertion")), which(test3)),subject_name] <- "match"
    df[intersect(which(!df[,subject_name] %in% c("gap", "insertion")), which(!test3)),subject_name] <- "mismatch"

    '    for (x in pattern_names) {
      df[,x] <- ifelse(df[,x] == df[,subject_name], "match", "mismatch")
      df[,x] <- ifelse(df[,x] == "mismatch" & df[,subject_name] == "-", "insertion", df[,x])
      df[,x] <- ifelse(df[,x] == "-" & df[,subject_name] != "-", "gap", df[,x])
    }
    df[,subject_name] <- ifelse(df[,subject_name] == "-", "gap", df[,subject_name])
'
    'if (matches_to_subject) {
      all_match_or_NA <- apply(df[,pattern_names,drop=F], 1, function(x) is.na(x) | all(x[which(!is.na(x))] == "match"), simplify = F)
      all_match_or_NA <- sapply(all_match_or_NA, all)
      df[,subject_name] <- ifelse(all_match_or_NA, "match", df[,subject_name])
      any_mismatch <- apply(df[,pattern_names,drop=F], 1, function(x) any(x[which(!is.na(x))] == "mismatch"), simplify = F)
      any_mismatch <- sapply(any_mismatch, any)
      df[,subject_name] <- ifelse(any_mismatch, "mismatch", df[,subject_name])
      any_insertion <- apply(df[,pattern_names,drop=F], 1, function(x) any(x[which(!is.na(x))] == "gap"), simplify = F)
      any_insertion <- sapply(any_insertion, any)
      df[,subject_name] <- ifelse(any_insertion, "insertion", df[,subject_name])
      df[,subject_name][which(df[,subject_name] == "-")] <- "gap"
    }'
  }

  df <-
    df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(unname(original_names)), names_to = "seq.name", values_to = "seq") %>%
    ## here, original names are restored
    dplyr::mutate(seq.name = original_names[seq.name]) %>%
    ## factor order with original names
    dplyr::mutate(seq.name = factor(seq.name, levels = unname(c(original_names[1], original_names[-1][ifelse(rep(order_patterns, length(subject.ranges)), order(purrr::map_int(subject.ranges, min)), pattern_order)]))))

  if (!is.null(pattern_groups)) {
    pattern_groups <- c(stats::setNames("subject", subject_name), pattern_groups)
    df$pattern.group <- pattern_groups[df$seq.name]
    df$pattern.group <- factor(df$pattern.group, levels = unique(pattern_groups))
  }
  return(df)
}



prep_subject_and_patterns <- function(subject,
                                      patterns,
                                      seq_type) {

  if (length(subject) > 1) {
    stop("Please provide only one subject as DNAString, DNAStringSet or character.")
  }
  if (is.list(subject)) {
    stop("Please provide one subject only.")
  }

  if (is.list(patterns)) {
    if (all(!sapply(patterns, methods::is, class2 = "RNAStringSet")) && all(!sapply(patterns, methods::is, class2 = "DNAStringSet")) && all(!sapply(patterns, methods::is, class2 = "AAStringSet")) && all(!sapply(patterns, methods::is, class2 = "character"))) {
      stop("patterns has to be a XStringSet or character vector.")
    }
  } else {
    if (!methods::is(patterns, "RNAStringSet") && !methods::is(patterns, "DNAStringSet") && !methods::is(patterns, "AAStringSet") && !methods::is(patterns, "character")) {
      stop("patterns has to be a XStringSet or character vector.")
    }
  }
  if (!methods::is(subject, "RNAStringSet") && !methods::is(subject, "DNAStringSet") && !methods::is(subject, "AAStringSet") && !methods::is(subject, "character")) {
    stop("subject has to be a XStringSet or character vector.")
  }

  # handel list of patterns
  patterns_list <- NULL
  pattern_groups <- NULL
  if (is.list(patterns) && length(patterns) > 1 && !all(lengths(patterns) == 1)) {
    patterns <- patterns[which(!sapply(patterns, is.null))]
    patterns <- patterns[which(!sapply(sapply(patterns, is.na, simplify = F), all))]

    patterns <- lapply(patterns, function(x) x[which(!is.na(x))])
    patterns_list <- patterns
    # will this be slow for many patterns? - is needed if list of DNAStringsSets - very hypothetical; maybe restrict what data types can be passed as patterns
    patterns <- unlist(lapply(patterns_list, as.character))

    if (is.null(names(patterns_list))) {
      names(patterns_list) <- paste0("Group_", seq_along(patterns_list))
    }
    pattern_groups <- lapply(patterns_list, names)
    pattern_groups <- stats::setNames(as.character(stack(pattern_groups)$ind), stack(pattern_groups)$values)
  } else if (is.list(patterns) && length(patterns) > 1 && all(lengths(patterns) == 1)) {
    # each list entry if length = 1
    patterns <- unlist(patterns)
  } else {
    patterns <- patterns[[1]]
    if (is.null(patterns) || all(is.na(patterns))) {
      stop("patterns is Null or NA.")
    }
  }

  # make this a separate fun somewhen?
  if (anyDuplicated(patterns)) {
    if (is.null(patterns_list)) {
      message("These pattern are duplicates: ", paste((unique(sapply(which(duplicated(patterns)), function(x) which(patterns[x] == patterns)))), collapse = ", "))
    } else {
      message("pattern at these indices are duplicates: ", paste(which(duplicated(patterns)), collapse = ", "))
      # similar fun as above but cylce through list and return list indices on top
      '      tt<-sapply(which(duplicated(patterns)), function(x) {
        unlist(purrr::discard(sapply(stats::setNames(seq_along(patterns_list), seq_along(patterns_list)), function(y) {
          which(patterns[x] == patterns_list[[y]])
        }), function(z) length(z) == 0))
      }, simplify = F)
      paste(tt)'

    }
    message("pattern at these indices are duplicates: ", paste(which(duplicated(patterns)), collapse = ", "))
  }

  ## pull seqs from subject and patterns, then run guess_type
  unique_letters <- unique(c(unlist(strsplit(as.character(subject), "")), unlist(strsplit(as.character(patterns), ""))))
  if (is.null(seq_type)) {
    seq_type <- guess_type(unique_letters)
  } else {
    seq_type <- match.arg(seq_type, c("NT", "AA"))
  }

  if (seq_type == "NT") {
    if ("U" %in% unique_letters) {
      if (!methods::is(patterns, "RNAStringSet")) {
        patterns <- Biostrings::RNAStringSet(patterns)
      }
      if (!methods::is(subject, "RNAStringSet")) {
        subject <- Biostrings::RNAStringSet(subject)
      }
    } else {
      if (!methods::is(patterns, "DNAStringSet")) {
        patterns <- Biostrings::DNAStringSet(patterns)
      }
      if (!methods::is(subject, "DNAStringSet")) {
        subject <- Biostrings::DNAStringSet(subject)
      }
    }
  } else if (seq_type == "AA") {
    if (!methods::is(patterns, "AAStringSet")) {
      patterns <- Biostrings::AAStringSet(patterns)
    }
    if (!methods::is(subject, "AAStringSet")) {
      subject <- Biostrings::AAStringSet(subject)
    }
  }

  if (is.null(names(subject))) {
    names(subject) <- "subject"
  }

  if (is.null(names(patterns))) {
    names(patterns) <- paste0("pattern_", seq(1,length(patterns)))
  }
  names(patterns) <- make.unique(names(patterns))

  ## save original names for replacement later
  ### avoid using names from pa then
  original_names <- c(names(subject), names(patterns))
  names(subject) <- make.names(names(subject))
  names(patterns) <- make.names(names(patterns))
  names(original_names) <- c(names(subject), names(patterns))

  # save original order in case filterings below shuffles it
  pattern_original_order <- names(patterns)


  # assigns in parent environment (https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r?rq=1)
  assign("subject", subject, envir = parent.frame())
  assign("patterns", patterns, envir = parent.frame())
  assign("seq_type", seq_type, envir = parent.frame())

  assign("original_names", original_names, envir = parent.frame())
  assign("pattern_original_order", pattern_original_order, envir = parent.frame())
  assign("pattern_groups", pattern_groups, envir = parent.frame())
}


check_for_invalid_chars <- function(subject,
                                    patterns,
                                    max_mismatch) {

  if (!is.na(max_mismatch) && max_mismatch < 0) {
    message("max_mismatch has to be NA or >= 0. Is set to NA now.")
    max_mismatch <- NA
  }

  patterns_invalid <- NULL
  pattern_mismatching_return <- NULL
  if (!is.na(max_mismatch) && methods::is(subject, "DNAStringSet") && methods::is(patterns, "DNAStringSet")) {
    if (grepl("[^ACTGU]", subject)) {
      message("subject contains non-DNA or non-RNA characters. To check for max_mismatch currenly only ACTGU are allowed. max_mismatch is now set to NA.")
      max_mismatch <- NA
    }
    inds <- grepl("[^ACTGU]", patterns)
    if (any(inds)) {
      message(length(which(inds)), " patterns with non-DNA or non-RNA characters detected. Those patterns are removed in order to allow checking for max_mismatch. They are returned as patterns_invalid.")
      patterns_invalid <- patterns[which(inds)]
      patterns <- patterns[which(!inds)]

      if (length(patterns) == 0) {
        stop("No patterns left after filtering for ones with valid DNA/RNA characters only. Please fix the sequences or set max_mismatch = NA.")
      }
    }

    ####
    ####
    pattern_mismatching <- purrr::map(stats::setNames(0:max_mismatch, 0:max_mismatch), function(x) {
      ## different lengths of patterns not allowed - split patterns by length; then have the names returned
      ## if too slow, think of other procedure
      pattern_names <- purrr::map(unique(nchar(as.character(patterns))), function(y) {
        patterns_temp <- patterns[which(nchar(as.character(patterns)) == y)]
        inds <- Biostrings::vwhichPDict(subject = subject,
                                        pdict = Biostrings::PDict(x = patterns_temp, max.mismatch = x),
                                        max.mismatch = x)[[1]]
        names(patterns_temp[inds])
      })
      return(unlist(pattern_names))
    })
    # last index of pattern_mismatching contains all patterns with the amount of mismatches at this index or less
    message(length(pattern_mismatching[[length(pattern_mismatching)]]), " of ", length(patterns), " patterns found to have less or equal to ",  max_mismatch, " mismatches with the subject.")

    print_pattern_mismatching <- utils::stack(lengths(pattern_mismatching))
    names(print_pattern_mismatching) <- c("n patterns", "mismatches")
    print(print_pattern_mismatching)

    pattern_mismatching_return <- purrr::map(stats::setNames(pattern_mismatching, paste0("max_mis_", names(pattern_mismatching))), function(x) patterns[x])
    patterns <- patterns[pattern_mismatching[[length(pattern_mismatching)]]]
    if (length(patterns) == 0) {
      stop("No pattern left after filtering for max_mismatch.")
    }
  } else if (!is.na(max_mismatch)) {
    message("!is.na(max_mismatch) can only be applied for DNA.")
    max_mismatch <- NA
  }

  assign("patterns", patterns, envir = parent.frame())
  assign("pattern_mismatching_return", pattern_mismatching_return, envir = parent.frame())
  assign("patterns_invalid", patterns_invalid, envir = parent.frame())
  assign("max_mismatch", max_mismatch, envir = parent.frame())
}

check_for_indel_induction <- function(pa,
                                      patterns,
                                      max_mismatch,
                                      rm_indel_inducing_pattern) {

  # caution: leading gaps are not considered as indels!
  #indel_list <- 0
  pattern_indel_inducing <- NULL
  if (is.na(max_mismatch) || max_mismatch > 0) {
    # based on pa, not pal!!
    #indel_lists <- Biostrings::nindel(pa) # this does not discriminate pattern and subject
    #indel_list <- stats::setNames(apply(cbind(indel_lists@insertion, indel_lists@deletion), 1, sum), names(patterns))
    # loop over pa elements to create list of matrices, then tell which pattern had gap in subject or pattern
    indel_mat_pattern <- lapply(seq_along(pa), function(z) as.matrix(pa[z]@pattern@indel@unlistData))
    indel_mat_subject <- lapply(seq_along(pa), function(z) as.matrix(pa[z]@subject@indel@unlistData))

    # this if FYI
    pattern_inds_indel <- which(unlist(lapply(indel_mat_pattern, nrow)) > 0)
    if (length(pattern_inds_indel) > 0) {
      message(sum(pattern_inds_indel > 0), " patterns got indels. Indices: ", paste(pattern_inds_indel, collapse = ", "))
    }

    subject_inds_indel <- which(unlist(lapply(indel_mat_subject, nrow)) > 0)
    subject_inds_indel_pass <- which(unlist(lapply(indel_mat_subject, nrow)) == 0)
    if (any(subject_inds_indel > 0)) {
      message(sum(subject_inds_indel > 0), " patterns caused indels in subject. Indices: ", paste(subject_inds_indel, collapse = ", "))
      if (rm_indel_inducing_pattern) {
        message("Those are removed as rm_indel_inducing_pattern = T. They are returned as pattern_indel_inducing.")
        pattern_indel_inducing <- patterns[subject_inds_indel]
        patterns <- patterns[subject_inds_indel_pass]
        pa <- pa[subject_inds_indel_pass]
        #subject_inds_indel <- subject_inds_indel[which(subject_inds_indel == 0)] # needed?
      }
    }
  }

  assign("pa", pa, envir = parent.frame())
  assign("patterns", patterns, envir = parent.frame())
  assign("pattern_indel_inducing", pattern_indel_inducing, envir = parent.frame())
  assign("subject_inds_indel", subject_inds_indel, envir = parent.frame())
}

check_for_overlapping_indels <- function(pa,
                                         patterns,
                                         subject_inds_indel,
                                         fix_subject_indels) {


  if (length(patterns) > 1 && any(subject_inds_indel > 0)) { # min 2 pattern and min 1 indel in subject
    # find out if any pattern alignment overlap with gaps from another pattern alignment. this would cause problem in the alignment.
    subject_indels <- as.data.frame(pa@subject@range)
    names(subject_indels) <- c("al_start", "al_end", "al_width")
    subject_indels$group <- 1:nrow(subject_indels)
    subject_indels <- dplyr::left_join(subject_indels, as.data.frame(pa@subject@indel)[,-2], by = "group")
    subject_indels$start <- subject_indels$start + subject_indels$al_start - 1
    subject_indels$end <- subject_indels$start + subject_indels$width - 1
    subject_indels$corr_end <- NA

    subject.ranges <- seq2(subject_indels$al_start, subject_indels$al_end) # this is the same as subject.ranges below; this is like seq2
    indel_ranges <- seq2(subject_indels$start[!is.na(subject_indels$start)], subject_indels$end[!is.na(subject_indels$end)])

    do_fix <- F
    if (!fix_subject_indels) {
      for (i in seq_along(subject.ranges)) {
        if (any(subject.ranges[[i]] %in% unlist(indel_ranges[-i]))) {
          ## allow for indel at same position
          # change order around %in% ?
          if (!all(subject.ranges[[i]][which(subject.ranges[[i]] %in% unlist(indel_ranges[-i]))] %in% unlist(indel_ranges[-i]))) {
            stop("Overlapping indel and subject alignment range found at index ", i, ". This cannot be handled yet, except for shortening respective sequences to just before the indel insertion.
                 To do so, set fix_subject_indels = T.")
          }
        }
      }
    }
    if (fix_subject_indels) {
      for (i in seq_along(subject.ranges)) {
        for (j in seq_along(indel_ranges)) {
          if (i != j) {
            if (length(intersect(subject.ranges[[i]],indel_ranges[[j]])) > 0) {
              do_fix <- T
              subject_indels[j,"corr_end"] <- subject_indels[j,"start"] - 1
            }
          }
        }
      }
    }


    if (do_fix && fix_subject_indels) {
      for (k in seq_along(patterns)) {
        if (any(!is.na(subject_indels[which(ind$group == k), "corr_end"]))) {
          message(names(patterns)[k], " is cut at position ", min(subject_indels[which(subject_indels$group == k), "corr_end"], na.rm = T), " to avoid indel overlap with another's pattern range on the subject. Experimental, yet.")
          patterns[k] <- Biostrings::subseq(patterns[k], start = 1, end = min(subject_indels[which(subject_indels$group == k), "corr_end"], na.rm = T))
        }
      }
      pa <- do.call(Biostrings::pairwiseAlignment, args = c(list(subject = subject, pattern = patterns, type = type),
                                                            pairwiseAlignment_args))
    }

    assign("pa", pa, envir = parent.frame())
    assign("patterns", patterns, envir = parent.frame())
    assign("subject_indels", subject_indels, envir = parent.frame())
    assign("indel_ranges", indel_ranges, envir = parent.frame())
  } else {
    assign("subject_indels", NULL, envir = parent.frame())
  }
}


make_pa_unique_and_order_and_rm_subset_alignments <- function(pa,
                                                              pattern_groups) {
  #subject.ranges <- brathering::seq2(pa@subject@range@start, pa@subject@range@start+pa@subject@range@width-1)
  #subject.ranges <- mapply("seq", pa@subject@range@start, pa@subject@range@start+pa@subject@range@width-1)
  subject.ranges <- seq2(pa@subject@range@start, pa@subject@range@start+pa@subject@range@width-1)
  names(subject.ranges) <- pa@pattern@unaligned@ranges@NAMES

  ## check if subject ranges within groups are overlapping

  if (!is.null(pattern_groups)) {
    subject.ranges.split <- split(names(subject.ranges), pattern_groups[names(subject.ranges)])
    overlap_subject_ranges <- unlist(lapply(subject.ranges.split, function(y) {
      ranges_comb <- combn(y, 2, simplify = F)
      if (any(unlist(lapply(ranges_comb, function(x) length(intersect(subject.ranges[[x[1]]], subject.ranges[[x[2]]])) > 1)))) {
        return(T)
      } else {
        return(F)
      }
    }))
    if (any(overlap_subject_ranges)) {
      stop("These groups have patterns with overlapping alignment ranges in subject: ", paste(names(which(overlap_subject_ranges)), collapse = ", "))
    }
  }

  non_dups <- which(!duplicated(subject.ranges))
  subject.ranges.unique <- subject.ranges[non_dups]
  pa.unique <- pa[non_dups]

  if (length(subject.ranges.unique) > 1) {
    is.subset <- purrr::map_lgl(seq_along(subject.ranges.unique), function(i) {
      any(purrr::map_lgl(seq_along(subject.ranges.unique), function (j) {
        if (i == j) {
          return(F)
        } else {
          #all(subject.ranges.unique[[i]] %in% subject.ranges.unique[[j]])
          identical(subject.ranges.unique[[i]], subject.ranges.unique[[j]])
        }
      }))
    })
  } else {
    is.subset <- rep(F, length(subject.ranges.unique))
  }
  not.subset <- which(!is.subset)

  subject.ranges.unique <- subject.ranges.unique[not.subset]
  pa.unique <- pa.unique[not.subset]

  # order alignment and subject ranges increasingly
  al_order <- order(purrr::map_int(subject.ranges.unique, min))
  pa.unique <- pa.unique[al_order]
  subject.ranges.unique <- subject.ranges.unique[al_order]

  pa <- pa[order(purrr::map_int(subject.ranges, min))]

  assign("subject.ranges", subject.ranges, envir = parent.frame()) # needed? - yes
  assign("subject.ranges.unique", subject.ranges.unique, envir = parent.frame())
  assign("pa.unique", pa.unique, envir = parent.frame())
  assign("pa", pa, envir = parent.frame())

}

paste_subject_seq <- function(subject,
                              subject.ranges.unique,
                              pa.unique) {

  # paste together the complete subject
  total.subject.seq <- stringr::str_sub(as.character(subject), 1, (min(subject.ranges.unique[[1]]) - 1))
  for (i in seq_along(subject.ranges.unique)) {
    # test if there is a gap between the i-1th and the ith alignment; if so, fill with original sequence
    if (i != 1 && max(subject.ranges.unique[[i-1]])+1 < min(subject.ranges.unique[[i]])) {
      # if yes use original subject
      total.subject.seq <- paste0(total.subject.seq, substr(as.character(subject), max(subject.ranges.unique[[i-1]])+1, min(subject.ranges.unique[[i]])-1))
    }
    if (i < max(seq_along(subject.ranges.unique))) {
      r <- min(subject.ranges.unique[[i]])
      total.subject.seq <- paste0(total.subject.seq, substr(pa.unique@subject[i], min(subject.ranges.unique[[i]])-r+1, min(subject.ranges.unique[[i+1]])-r))
    }
    if (i == max(seq_along(subject.ranges.unique))) {
      total.subject.seq <- paste0(total.subject.seq, pa.unique@subject[i])
    }
  }

  ## attach the remaining sequence from subject
  total.subject.seq <- paste0(total.subject.seq, substr(as.character(subject), max(subject.ranges.unique[[length(subject.ranges.unique)]])+1, nchar(as.character(subject))))
  ## create data frame for plotting
  df <- dplyr::mutate(data.frame(seq = stats::setNames(strsplit(total.subject.seq, ""), names(subject))), position = dplyr::row_number())
  df[df[,names(subject)] != "-", "subject.position"] <- seq(1:nrow(df[df[,names(subject)] != "-", ]))

  return(df)
}


paste_patterns_to_subject <- function(subject_indels,
                                      patterns,
                                      pa,
                                      df,
                                      pattern_original_order) {
  # gaps only account for the next sequence, respectively, hence add 0 at beginning, and delete last index
  # these lines assumed that indels are not overlapping with the subsequent pattern
  #gaps <- c(0, Biostrings::nindel(pa)@insertion[,"WidthSum"])
  #gap_corr <- purrr::accumulate(gaps[-length(gaps)], `+`)

  # check if the subject_indel is at the same position as previous, then do not increment gap_corr at next index
  if (!is.null(subject_indels)) {
    subject_indels$gap_insert <- subject_indels$width
    for (i in 1:(nrow(subject_indels)-1)) {
      subject_indels$gap_insert[i] <- ifelse((subject_indels$start[i+1] == subject_indels$start[i] && subject_indels$end[i+1] == subject_indels$end[i]), 0, subject_indels$gap_insert[i])
    }
    gaps <- c(0, subject_indels$gap_insert)
    gap_corr <- purrr::accumulate(gaps[-length(gaps)], `+`)
  } else {
    gap_corr <- 0
  }

  start <- pa@subject@range@start
  alPa <- stats::setNames(as.character(pa@pattern), pa@pattern@unaligned@ranges@NAMES)
  seq <- stack(strsplit(alPa, ""))
  names(seq) <- c("seq", "pattern")
  position_var <- seq2((start + gap_corr), (start+nchar(alPa) - 1 + gap_corr))
  position = stack(stats::setNames(position_var, names(alPa)))
  names(position) <- c("position", "pattern")

  dfs <- cbind(seq[,"seq",drop=F], position)
  dfs <- split(dfs, dfs$pattern)
  dfs <- purrr::map(dfs, function(x) {
    names(x)[1] <- unique(as.character(x[,"pattern",drop=T]))
    x <- x[-which(names(x) == "pattern")]
    return(x)
  })

  dfs <- join_chunkwise(data_frame_list = dfs, join_by = "position")
  # then left join with complete seq in index 1
  df <- purrr::reduce(c(list(df), dfs), dplyr::left_join, by = "position")
  ## order patterns in data.frames below by factor order
  pattern_order <- match(names(patterns), pattern_original_order)

  assign("df", df, envir = parent.frame())
  assign("subject_indels", subject_indels, envir = parent.frame())
  assign("pattern_order", pattern_order, envir = parent.frame())

}

seq2_default <- Vectorize(seq.default, vectorize.args = c("from", "to"))
seq2 <- function(from = 1, to = 1) {
  x <- seq2_default(from = from, to = to)
  # make sure always a list is returned
  if (is.matrix(x)) {
    x <- unname(as.list(as.data.frame(x)))
  }
  return(x)
}


join_chunkwise <- function(data_frame_list,
                           join_fun = dplyr::full_join,
                           join_by,
                           chunks = 10,
                           max_final_size = 20) {
  # chunk wise joining can be much faster
  join_fun <- match.fun(join_fun)
  while (length(data_frame_list) > max_final_size) {
    data_frame_list <- purrr::map(split(c(1:length(data_frame_list)), ceiling(seq_along(c(1:length(data_frame_list)))/chunks)),
                                  function(x) purrr::reduce(data_frame_list[x], join_fun, by = join_by))
  }
  return(data_frame_list)
}
