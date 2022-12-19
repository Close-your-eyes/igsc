#' Separately align multiple pattern sequences to one subject
#'
#' This function is useful to visualise the alignment position of multiple patterns on one subject.
#' It uses Biostrings::pairwiseAlignment(), obtains the individual alignment boundaries and converts
#' the results to a ggplot object.
#' The function will fail if a gap is induced in the subject and at least two pattern alignments overlap at this gap.
#'
#' @param subject a named character or named DNAStringSet of one subject (only the DNAStringSet but not DNAString can hold a name)
#' @param patterns a named character vector or named DNAStringSet of patterns to align to the subject sequence
#' @param type the type of alignment passed to Biostrings::pairwiseAlignment; not every type may work well with this function (if there are overlapping ranges of the alignments to the subject for example)
#' @param nt_suffix add the length of the string to the name on the axis
#' @param order_patterns order pattern increasingly by alignment position (start)
#' @param max_mismatch only use patterns that have a maximum number of mismatches
#' with the subject
#' @param fix_indels in case of overlapping indels and shared subject ranges, cut respective patterns to avoid indels
#' @param ... additional arguments to Biostrings::pairwiseAlignment apart from subject, pattern and type
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
                                                   nt_suffix = T,
                                                   fix_indels = F,
                                                   rm_indel_inducing_pattern = F,
                                                   seq_type = NULL,
                                                   return_max_mismatch_info_only = F,
                                                   ...) {

  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }

'  if (fix_indels) {
    message("fix_indels does not work yet. Set to F.")
    fix_indels <- F
  }'

  type <- match.arg(type, choices = c("global-local", "global", "local", "overlap", "local-global"))

  if (!is.na(max_mismatch) && max_mismatch < 0) {
    message("max_mismatch has to be NA or >= 0. Is set to NA now.")
    max_mismatch <- NA
  }

  ## pull seqs from subject and patterns, then run guess_type

  if (is.null(seq_type)) {
    unique_seq_el <- unique(c(unlist(strsplit(as.character(subject), "")), unlist(strsplit(as.character(patterns), ""))))
    seq_type <- guess_type(unique_seq_el)
    if (seq_type == "NT") {
      if ("U" %in% unique_seq_el) {
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
    } else {
      stop("AA or NT not determined.")
    }
  }

  if (is.null(names(subject))) {
    names(subject) <- "subject"
  }
  names(subject) <- make.names(names(subject))
  if (is.null(names(patterns))) {
    names(patterns) <- paste0("pattern_", seq(1,length(patterns)))
  }
  names(patterns)[which(is.na(names(patterns)))] <- paste0("pattern_", which(is.na(names(patterns))))
  names(patterns) <- make.unique(names(patterns))
  names(patterns) <- make.names(names(patterns))

  if (nt_suffix) {
    names(subject) <- paste0(names(subject), "_", nchar(as.character(subject)), "nt")
    names(patterns) <- paste0(names(patterns), "_", sapply(as.character(patterns), function(x) nchar(x)), "nt")
  }


  # check for non-DNA characters first
  patterns_invalid <- NULL
  pattern_mismatching_return <- NULL
  if (!is.na(max_mismatch) && methods::is(subject, "DNAStringSet") && methods::is(patterns, "DNAStringSet")) {
    if (grepl("[^ACTGU]", subject)) {
      message("subject contains non-DNA or non-RNA characters. To check for max_mismatch currenly only ACTGU are allowed. max_mismatch is now set to NA.")
      max_mismatch <- NA
    }
    if (any(grepl("[^ACTGU]", patterns))) {
      message(length(which(grepl("[^ACTGU]", patterns))), " patterns with non-DNA or non-RNA characters detected. Those patterns are removed in order to allow checking for max_mismatch. They are returned as patterns_invalid.")
      #message(paste(names(patterns[which(grepl("[^ACTGU]", patterns))]), collapse = "\n"))
      patterns_invalid <- patterns[which(grepl("[^ACTGU]", patterns))]
      patterns <- patterns[which(!grepl("[^ACTGU]", patterns))]
      if (return_max_mismatch_info_only) {
        return(list(patterns_invalid = patterns_invalid,
                    pattern_mismatching = pattern_mismatching_return))
      }
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

  # calculate all alignments
  pa <- Biostrings::pairwiseAlignment(subject = subject, pattern = patterns, type = type, ...)

  # make pa a list once and then iterate over list entries with purrr/furrr which is quicker!
  # use multiple threads to speed up?!
  # pal <- stats::setNames(as.list(pa), patnames(patterns)terns.names)
  # pal <- stats::setNames(purrr::flatten(parallel::mclapply(split(c(1:length(pa)), ceiling(seq_along(c(1:length(pa)))/10)), function(x) as.list(pa[x]), mc.cores = parallel::detectCores()-1)), names(patterns))

  # check for indels induced in the subject; this is slow
  # caution: leading gaps are not considered as indels!
  pattern_indel_inducing <- NULL
  indel_list <- 0
  if (is.na(max_mismatch) || max_mismatch > 0) {
    # based on pa, not pal!!
    indel_lists <- Biostrings::nindel(pa)
    indel_list <- stats::setNames(apply(cbind(indel_lists@insertion, indel_lists@deletion), 1, sum), names(patterns))
    # based on list
    #indel_list <- purrr::map_int(stats::setNames(pal, names(patterns)), function(x) length(x@subject@indel@unlistData@start))

    if (any(indel_list > 0)) {
      message(sum(indel_list > 0), " patterns caused indels in the subject.")
      indel_df <- utils::stack(indel_list)
      names(indel_df) <- c("n indel", "pattern name")
      #print(indel_df[which(indel_df[,"n indel"] > 0), ])
      if (rm_indel_inducing_pattern) {
        message("Those are removed as rm_indel_inducing_pattern = T. They are returned as pattern_indel_inducing.")
        pattern_indel_inducing <- patterns[which(names(patterns) %in% names(which(indel_list > 0)))]
        patterns <- patterns[which(!names(patterns) %in% names(which(indel_list > 0)))]
        ## filter by name
        #pal <- pal[which(names(pal) %in% names(which(indel_list == 0)))]
        ## filter by index
        #pal <- pal[which(indel_list == 0)]
        pa <- pa[which(indel_list == 0)]
        indel_list <- indel_list[which(indel_list == 0)]
      }
    }
  }

  if (any(indel_list > 0)) {
    # find out if any pattern alignment overlap with gaps from another pattern alignment. this would cause problem in the alignment.
    ind <- as.data.frame(pa@subject@range)
    names(ind) <- c("al_start", "al_end", "al_width")
    ind$group <- 1:nrow(ind)
    ind <- dplyr::left_join(ind, as.data.frame(pa@subject@indel)[,-2], by = "group")
    ind$indel_start <- ind$start + ind$al_start - 1
    ind$indel_end <- ind$indel_start + ind$width - 1
    ind$corr_end <- NA

    als <- mapply("seq", ind$al_start, ind$al_end, SIMPLIFY = F)
    inds <- mapply("seq", ind$indel_start[!is.na(ind$indel_start)], ind$indel_end[!is.na(ind$indel_end)], SIMPLIFY = F)

    do_fix <- F
    for (i in seq_along(als)) {
      for (j in seq_along(inds)) {
        if (i != j) {
          if (length(intersect(als[[i]],inds[[j]])) > 0) {
            if (!fix_indels) {
              warning("Overlapping indel and subject alignment range found. This cannot be handled yet, except for shortening respective sequences to just before the indel insertion.
                 To do so, set fix_indels = T.")
              return(NULL)
            }
            do_fix <- T
            ind[j,"corr_end"] <- ind[j,"start"] - 1
          }
        }
      }
    }

    if (do_fix && fix_indels) {
      for (k in seq_along(patterns)) {
        if (any(!is.na(ind[which(ind$group == k), "corr_end"]))) {
          message(names(patterns)[k], " is cut at position ", min(ind[which(ind$group == k), "corr_end"], na.rm = T), " to avoid indel overlap with another's pattern range on the subject. Experimental, yet.")
          patterns[k] <- Biostrings::subseq(patterns[k], start = 1, end = min(ind[which(ind$group == k), "corr_end"], na.rm = T))
        }
      }
      pa <- Biostrings::pairwiseAlignment(subject = subject, pattern = patterns, type = type, ...)
      #pal <- stats::setNames(purrr::flatten(parallel::mclapply(split(c(1:length(pa)), ceiling(seq_along(c(1:length(pa)))/10)), function(x) as.list(pa[x]), mc.cores = parallel::detectCores()-1)), names(patterns))
    }
  }

  # get ranges
  subject.ranges <- purrr::map(split(data.frame(pa@subject@range), seq(nrow(data.frame(pa@subject@range)))), function(x) x$start:x$end)
  subject.ranges.unique <- subject.ranges[which(!duplicated(subject.ranges))]
  #pal.unique <- pal[which(!duplicated(subject.ranges))]
  pa.unique <- pa[which(!duplicated(subject.ranges))]

  if (length(subject.ranges.unique) > 1) {
    is.subset <- purrr::map_lgl(seq_along(subject.ranges.unique), function(i) {
      any(purrr::map_lgl(seq_along(subject.ranges.unique), function (j) {
        if (i == j) {
          return(F)
        } else {
          all(subject.ranges.unique[[i]] %in% subject.ranges.unique[[j]])
        }
      }))
    })
  } else {
    is.subset <- rep(F, length(subject.ranges.unique))
  }
  subject.ranges.unique <- subject.ranges.unique[which(!is.subset)]
  pa.unique <- pa.unique[which(!is.subset)]
  #pal.unique <- pal.unique[which(!is.subset)]

  # order alignment and subject ranges increasingly
  order_temp <- order(purrr::map_int(subject.ranges.unique, min))
  #pal.unique <- pal.unique[order_temp]
  pa.unique <- pa.unique[order_temp]
  subject.ranges.unique <- subject.ranges.unique[order_temp]

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
  gap.corr <- 0

  #pal <- pal[order(purrr::map_int(subject.ranges, min))]
  pa <- pa[order(purrr::map_int(subject.ranges, min))]

  # this is slow; how to speed up?
  '  for (x in 1:length(pa)) {
    print(gap.corr)
    temp <- data.frame(seq = strsplit(as.character(Biostrings::alignedPattern(pa[x])), ""),
                       position = (pa[x]@subject@range@start + gap.corr):(pa[x]@subject@range@start+nchar(as.character(Biostrings::alignedPattern(pa[x]))) - 1 + gap.corr))
    gap.corr <- gap.corr + sum(data.frame(pa[x]@subject@indel@unlistData)$width)
    df <- dplyr::left_join(df, temp, by = "position")
  }'

  #gap_corr <- purrr::accumulate(purrr::map_int(pal, function(x) sum(data.frame(x@subject@indel@unlistData)$width)), `+`)
  #which(Biostrings::nindel(pa)@deletion[,"WidthSum"] != 0) # what about deletion?

  gap_corr <- purrr::accumulate(Biostrings::nindel(pa)@insertion[,"WidthSum"], `+`)

  # outdated procedures
  '  dfs <- purrr::map2(pal, gap_corr, function(x, y) {
    alPa <- as.character(Biostrings::alignedPattern(x))
    start <- x@subject@range@start
    data.frame(seq = strsplit(alPa, ""), position = (start + y):(start+nchar(alPa) - 1 + y))
  })
'
'  dfs <- parallel::mcmapply(x = pal, y = gap_corr, FUN = function(x, y) {
    #alPa <- as.character(Biostrings::alignedPattern(x))
    # Biostrings::alignedPattern(x) returns leading NT that are aligned to gaps; whereas x@pattern does do it
    alPa <- stats::setNames(as.character(x@pattern), x@pattern@unaligned@ranges@NAMES)
    start <- x@subject@range@start
    data.frame(seq = strsplit(alPa, ""), position = (start + y):(start+nchar(alPa) - 1 + y))
  }, mc.cores = parallel::detectCores()-1, SIMPLIFY = F)'

  seq_vectorized <- Vectorize(seq.default, vectorize.args = c("from", "to"))

  start <- pa@subject@range@start
  alPa <- stats::setNames(as.character(pa@pattern), pa@pattern@unaligned@ranges@NAMES)


  seq <- stack(strsplit(alPa, ""))
  names(seq) <- c("seq", "pattern")

  # how to predefine output of seq_vectorized?
  position_var <- seq_vectorized(from = (start + gap_corr), to = (start+nchar(alPa) - 1 + gap_corr))
  if (methods::is(position_var, "matrix")) {
    position_var <- apply(position_var, 2, c, simplify = F)
  }
  position = stack(stats::setNames(position_var, names(alPa)))
  names(position) <- c("position", "pattern")

  dfs <- cbind(seq[,"seq",drop=F], position)
  dfs <- split(dfs, dfs$pattern)
  dfs <- purrr::map(dfs, function(x) {
    names(x)[1] <- unique(as.character(x[,"pattern",drop=T]))
    x <- x[-which(names(x) == "pattern")]
    return(x)
  })

  # do joining chunk-wise !!!
  # make that a separate function somewhen
  # dfs_join <- purrr::reduce(dfs, dplyr::full_join, by = "position")
  while (length(dfs) > 20) {
    dfs <- purrr::map(split(c(1:length(dfs)), ceiling(seq_along(c(1:length(dfs)))/10)), function(x) purrr::reduce(dfs[x], dplyr::full_join, by = c("position")))
  }
  # then left join with complete seq in index 1
  df <- purrr::reduce(c(list(df), dfs), dplyr::left_join, by = "position")

  df.match <- df
  for (x in names(patterns)) {
    df.match[,x] <- ifelse(df.match[,x] == df.match[,names(subject)], "match", ifelse(df.match[,x] == "-", "-", "mismatch"))
    df.match[,x] <- ifelse(df.match[,x] == "mismatch" & df.match[,names(subject)] == "-", "insertion", df.match[,x])
    df.match[,x] <- ifelse(df.match[,x] == "-" & df.match[,names(subject)] != "-", "gap", df.match[,x])
  }
  df.match[,names(subject)] <- ifelse(df.match[,names(subject)] == "-", "gap", df.match[,names(subject)])

  df <-
    df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(c(names(subject), names(patterns))), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(names(subject), names(patterns)[ifelse(rep(order_patterns, length(subject.ranges)), order(purrr::map_int(subject.ranges, min)), seq(1,length(subject.ranges)))])))

  df.match <-
    df.match %>%
    tidyr::pivot_longer(cols = dplyr::all_of(c(names(subject), names(patterns))), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(names(subject), names(patterns)[ifelse(rep(order_patterns, length(subject.ranges)), order(purrr::map_int(subject.ranges, min)), seq(1,length(subject.ranges)))])))

  g1 <- algnmt_plot(algnmt = df,
                    tile.border.color = NA,
                    font.family = "sans",
                    pattern.lim.size = 2,
                    pa = pa,
                    subject.lim.lines = F,
                    algnmt_type = seq_type)

  g2 <- algnmt_plot(algnmt = df.match,
                    tile.border.color = NA,
                    font.family = "sans",
                    pattern.lim.size = 2,
                    pa = pa,
                    subject.lim.lines = F,
                    algnmt_type = seq_type)

  return(list(base.plot = g1,
              match.plot = g2,
              algnmt.df.base = df,
              algnmt.df.match = df.match,
              min.max.subject.position = c(df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(subject.position),
                                           df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(subject.position)),
              pairwise_alignments = pa,
              #pairwise_alignment_list = pal,
              pattern = if(order_patterns) {patterns[order(purrr::map_int(subject.ranges, min))]} else {patterns},
              pattern_invalid = patterns_invalid,
              pattern_indel_inducing = pattern_indel_inducing,
              pattern_mismatching = pattern_mismatching_return))
}

