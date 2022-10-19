#' Align multiple DNA pattern sequences to one subject
#'
#' This function is useful plot the alignment position of multiple patterns in one subject.
#' It uses Biostrings::pairwiseAlignment(), obtains the individual alignment boundaries and converts
#' the results to a ggplot object with a few optional complementing information (see arguments).
#' The function will fail if a gap is induced in the subject and at least two pattern alignments overlap at this gap.
#'
#' @param subject a named character or named DNAStringSet of one subject (only the DNAStringSet but not DNAString can hold a name)
#' @param patterns a named character vector or named DNAStringSet of patterns to align to the subject sequence
#' @param type the type of alignment passed to Biostrings::pairwiseAlignment; not every type may work well with this function (if there are overlapping ranges of the alignments to the subject for example)
#' @param attach.nt add the length of the string to the name on the axis
#' @param order.patterns order pattern increasingly by alignment position (start)
#' @param perfect.matches.only filter patterns for those which match the subject without gaps, insertions or substitutions before pairwise alignment
#' @param fix_indels in case of overlapping indels and shared subject ranges, cut respective patterns to avoid indels
#' @param ... additional arguments to Biostrings::pairwiseAlignment apart from subject, pattern and type
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
#' als <- igsc::MultiplePairwiseAlignmentsToOneSubject(subject = s, patterns = p, tile.border.color = "black")
#' als_ordered <- igsc::MultiplePairwiseAlignmentsToOneSubject(subject = s, patterns = p, tile.border.color = "black", order.patterns = T)
#' }
MultiplePairwiseAlignmentsToOneSubject <- function(subject,
                                                   patterns,
                                                   type = "global-local",
                                                   perfect.matches.only = F,
                                                   order.patterns = F,
                                                   attach.nt = T,
                                                   fix_indels = F,
                                                   ...) {

  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }


  ## to do: also allow AAString, perform check
  if (!methods::is(patterns, "DNAStringSet")) {
    patterns <- Biostrings::DNAStringSet(patterns)
  }
  if (!methods::is(subject, "DNAStringSet")) {
    subject <- Biostrings::DNAStringSet(subject)
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

  if (attach.nt) {
    names(subject) <- paste0(names(subject), "_", nchar(as.character(subject)), "nt")
    names(patterns) <- paste0(names(patterns), "_", sapply(as.character(patterns), function(x) nchar(x)), "nt")
  }

  patterns.names <- names(patterns)
  subject.name <- names(subject)

  type <- match.arg(type, choices = c("global", "local", "overlap", "global-local", "local-global"))
  font.family <- match.arg(font.family, choices = c("sans", "mono", "serif"))

  if (perfect.matches.only) {
    perf <- Biostrings::vwhichPDict(subject = subject, pdict = Biostrings::PDict(patterns))[[1]]
    print(paste0(length(perf), " of ", length(patterns), " patterns found to perfectly match the subject."))
    patterns <- patterns[perf]
    if (length(patterns) == 0) {
      stop("No pattern with perfect match left.")
    }
  }


  # calculate all alignments
  pa <- Biostrings::pairwiseAlignment(subject = subject, pattern = patterns, type = type)

  # check for indels induced in the subject
  for (i in seq_along(pa)) {
    if (length(pa[i]@subject@indel@unlistData@start) > 0) {
      print(paste0(pa[i]@pattern@unaligned@ranges@NAMES, " caused ", length(pa[i]@subject@indel@unlistData@start), " indel(s) in the subject."))
    }
  }

  # find out if any pattern alignment overlap with gaps from another pattern alignment. this would cause problem in the alignment.
  ind <- as.data.frame(pa@subject@range)
  names(ind) <- c("al_start", "al_end", "al_width")
  ind$group <- 1:nrow(ind)
  ind <- dplyr::left_join(ind, as.data.frame(pa@subject@indel)[,-2], by = "group")
  ind$indel_start <- ind$start + ind$al_start - 1
  ind$indel_end <- ind$indel_start + ind$width - 1
  ind$corr_end <- NA

  als <- mapply(seq, ind$al_start, ind$al_end, SIMPLIFY = F)
  inds <- mapply(seq, ind$indel_start[!is.na(ind$indel_start)], ind$indel_end[!is.na(ind$indel_end)], SIMPLIFY = F)

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
        print(paste0(names(patterns)[k], " is cut at position ", min(ind[which(ind$group == k), "corr_end"], na.rm = T), " to avoid indel overlap with another's pattern range on the subject. Experimental, yet."))
        patterns[k] <- Biostrings::subseq(patterns[k], start = 1, end = min(ind[which(ind$group == k), "corr_end"], na.rm = T))
      }
    }
    pa <- Biostrings::pairwiseAlignment(subject = subject, pattern = patterns, type = type, ...)
  }

  # get ranges
  subject.ranges <- lapply(split(data.frame(pa@subject@range), seq(nrow(data.frame(pa@subject@range)))), function (x) {x$start:x$end})
  subject.ranges.unique <- subject.ranges[which(!duplicated(subject.ranges))]
  pa.unique <- pa[which(!duplicated(subject.ranges))]

  # maybe make it mapply
  if (length(subject.ranges.unique) > 1) {
    is.subset <- unlist(lapply(seq_along(subject.ranges.unique), function(i) {
      any(unlist(lapply(seq_along(subject.ranges.unique), function (j) {
        if (i == j) {
          return(F)
        } else {
          all(subject.ranges.unique[[i]] %in% subject.ranges.unique[[j]])
        }
      })))
    }))
  } else {
    is.subset <- rep(F, length(subject.ranges.unique))
  }
  subject.ranges.unique <- subject.ranges.unique[which(!is.subset)]
  pa.unique <- pa.unique[which(!is.subset)]

  # order alignment and subject ranges increasingly
  pa.unique <- pa.unique[order(sapply(subject.ranges.unique, function(x) min(x)))]
  subject.ranges.unique <- subject.ranges.unique[order(sapply(subject.ranges.unique, function(x) min(x)))]

  # paste together the complete subject
  total.subject.seq <- stringr::str_sub(subject, 1, (min(subject.ranges.unique[[1]]) - 1))
  for (i in seq_along(subject.ranges.unique)) {
    # test if there is a gap between the i-1th and the ith alignment; if so, fill with original sequence
    if (i != 1 && max(subject.ranges.unique[[i-1]])+1 < min(subject.ranges.unique[[i]])) {
      # if yes use original subject
      total.subject.seq <- paste0(total.subject.seq, substr(subject, max(subject.ranges.unique[[i-1]])+1, min(subject.ranges.unique[[i]])-1))
    }
    if (i < max(seq_along(subject.ranges.unique))) {
      r <- min(subject.ranges.unique[[i]])
      total.subject.seq <- paste0(total.subject.seq, substr(pa.unique[i]@subject, min(subject.ranges.unique[[i]])-r+1, min(subject.ranges.unique[[i+1]])-r))
    }
    if (i == max(seq_along(subject.ranges.unique))) {
      total.subject.seq <- paste0(total.subject.seq, pa.unique[i]@subject)
    }
    #print(total.subject.seq)
  }
  ## attach the remaining sequence from subject
  total.subject.seq <- paste0(total.subject.seq, substr(subject, max(subject.ranges.unique[[length(subject.ranges.unique)]])+1, nchar(as.character(subject))))


  ## create data frame for plotting
  df <- dplyr::mutate(data.frame(seq = stats::setNames(strsplit(total.subject.seq, ""), subject.name)), position = dplyr::row_number())
  df[df$seq != "-", "subject.position"] <- seq(1:nrow(df[df$seq != "-", ]))
  gap.corr <- 0
  pa <- pa[order(purrr::map_int(subject.ranges, min))]
  for (x in 1:length(pa)) {
    temp <- data.frame(seq = strsplit(as.character(Biostrings::alignedPattern(pa[x])), ""),
                       position = (pa[x]@subject@range@start + gap.corr):(pa[x]@subject@range@start+nchar(as.character(Biostrings::alignedPattern(pa[x]))) - 1 + gap.corr))
    gap.corr <- gap.corr + sum(data.frame(pa[x]@subject@indel@unlistData)$width)
    df <- dplyr::left_join(df, temp, by = "position")
  }

  df.match <- df
  for (x in patterns.names) {
    df.match[,x] <- ifelse(df.match[,x] == df.match[,subject.name], "match", ifelse(df.match[,x] == "-", "-", "mismatch"))
    df.match[,x] <- ifelse(df.match[,x] == "mismatch" & df.match[,subject.name] == "-", "insertion", df.match[,x])
    df.match[,x] <- ifelse(df.match[,x] == "-" & df.match[,subject.name] != "-", "gap", df.match[,x])
  }
  df.match[,subject.name] <- ifelse(df.match[,subject.name] == "-", "gap", df.match[,subject.name])

  acp1 <- acp[which(names(acp) %in% unique(df[,subject.name]))]
  df <-
    df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(c(subject.name, patterns.names)), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq = factor(seq, levels = names(acp1))) %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names[ifelse(rep(order.patterns, length(subject.ranges)), order(purrr::map_int(subject.ranges, min)), seq(1,length(subject.ranges)))])))

  acp2 <- acp[which(names(acp) %in% unique(unlist(df.match[,c(subject.name, patterns.names)])))]
  df.match <-
    df.match %>%
    tidyr::pivot_longer(cols = dplyr::all_of(c(subject.name, patterns.names)), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq = factor(seq, levels = names(acp2))) %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names[ifelse(rep(order.patterns, length(subject.ranges)), order(purrr::map_int(subject.ranges, min)), seq(1,length(subject.ranges)))])))

  g1 <- algnmt_plot(algnmt = df,
                    tile.border.color = NA,
                    font.family = "sans",
                    pattern.lim.size = 2,
                    pa = pa,
                    subject.lim.lines = F)

  g2 <- algnmt_plot(algnmt = df.match,
                    tile.border.color = NA,
                    font.family = "sans",
                    pattern.lim.size = 2,
                    pa = pa,
                    subject.lim.lines = F)

  return(list(base.plot = g1,
              match.plot = g2,
              base.df = df,
              match.df = df.match,
              min.max.subject.position = c(df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(subject.position),
                                           df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(subject.position))))
}

.integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    return(breaks[intersect(which(breaks > 0), which(breaks < max(x)))])
  }
  return(fxn)
}

algnmt_plot <- function(algnmt,
                        color_values = NULL,
                        tile.border.color = NA,
                        tile.border.on.NA = F,
                        text = F,
                        font.family = c("sans", "mono", "serif"),
                        theme = ggplot2::theme_bw(base_family = font.family),
                        legend.position = "none",
                        pattern.lim.size = 2,
                        pa = NULL,
                        subject.lim.lines = F,
                        ylab = "seq name",
                        legend.title = "",
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        coord_fixed_ratio = NULL,
                        x.breaks = NULL,
                        ...) {

  # color_values see igsc:::scheme_NT; igsc:::scheme_AA

  # coord_fixed_ratio
  # if !is.null(coord_fixed_ratio)
  # example: length(unique(algnmt[,pos_col,drop=T]))/length(unique(algnmt[,name_col,drop=T]))*0.5


  ## checks
  algnmt_type <- NULL
  if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet") || methods::is(algnmt, "AAStringSet")) {
    if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet")) {
      algnmt_type <- "NT"
    }
    if (methods::is(algnmt, "AAStringSet")) {
      algnmt_type <- "AA"
    }
    algnmt <- XStringSet_to_df(algnmt)
  }

  if (!all(c(pos_col, seq_col, name_col) %in% names(algnmt))) {
    stop("algnmt at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ".Alternative change function arguments.")
  }

  if (is.null(algnmt_type)) {
    inds <- intersect(which(!is.na(algnmt[,seq_col,drop=T])), which(!algnmt[,seq_col,drop=T] %in% c("match", "mismatch", "gap", "insertion", "ambiguous")))
    if (all(algnmt[,seq_col,drop=T][inds] %in% unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N")))) {
      algnmt_type <- "NT"
    } else if (all(algnmt[,seq_col,drop=T][inds] %in% c(Biostrings::AA_ALPHABET, "N"))) {
      algnmt_type <- "AA"
    }
  }

  # use preset colors
  if (is.null(color_values)) {
    if (algnmt_type == "NT") {
      color_values <- names(scheme_NT)[1]
    }
    if (algnmt_type == "AA") {
      color_values <- names(scheme_AA)[1]
    }
  }
  if (length(color_values) == 1) {
    if (algnmt_type == "NT") {
      color_values <- stats::setNames(scheme_NT[,match.arg(color_values, choices = names(scheme_NT)),drop=T], rownames(scheme_NT))
    } else if (algnmt_type == "AA") {
      color_values <- stats::setNames(scheme_AA[,match.arg(color_values, choices = names(scheme_AA)),drop=T], rownames(scheme_AA))
    } else {
      message("Type of alignment data (NT or AA) could not be determined. Choosing default ggplot colors.")
      color_values <- scales::hue_pal()(length(unique(as.character(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))))
    }
  } else {
    # provide own colors
    if (is.null(names(color_values))) {
      if (length(color_values) < length(unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))) {
        stop("Not enough color_values provided.")
      }
    } else {
      if (any(!unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]) %in% names(color_values))) {
        stop("Not all values in algnmt[,seq_col] found in names(color_values).")
      }
    }
  }

  # make sure it is not a factor
  algnmt[,pos_col,drop=T] <- as.numeric(as.character(algnmt[,pos_col,drop=T]))

  # https://stackoverflow.com/questions/45493163/ggplot-remove-na-factor-level-in-legend
  # --> na.translate = F

  #ex<-names(as.list(formals(ggplot2::theme)))
  if (is.null(x.breaks)) {
    x.breaks <- floor(base::pretty(algnmt[,pos_col,drop=T], n = 5))
    x.breaks <- x.breaks[which(x.breaks > 0)]
    x.breaks <- x.breaks[which(x.breaks < max(algnmt[,pos_col,drop=T]))]
  }

  plot <-
    ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(name_col))) +
    theme +
    ggplot2::theme(legend.position = legend.position,
                   panel.grid = ggplot2::element_blank(),
                   text = ggplot2::element_text(family = font.family),
                   ...) +
    ggplot2::labs(y = ylab, fill = legend.title) +
    ggplot2::scale_fill_manual(values = color_values, na.value = "white", na.translate = F) +
    ggplot2::scale_x_continuous(breaks = x.breaks)

  if (!is.na(tile.border.color)) {
    plot <- plot + ggplot2::geom_tile(data = if(tile.border.on.NA) {algnmt} else {algnmt[which(!is.na(algnmt[,seq_col,drop=T])), ]},
                                      color = tile.border.color,
                                      ggplot2::aes(fill = !!rlang::sym(seq_col)))
  } else {
    plot <- plot + ggplot2::geom_tile(ggplot2::aes(fill = !!rlang::sym(seq_col)))
  }

  if (text) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = !!rlang::sym(seq_col)), na.rm = T, family = font.family)
  }

  if (!is.null(coord_fixed_ratio)) {
    plot <- plot + ggplot2::coord_fixed(ratio = coord_fixed_ratio)
  }

  if (pattern.lim.size > 0 && !is.null(pa)) {
    pattern.ranges <- data.frame(pa@pattern@range, seq.name = "")
    for (i in 1:length(pa)) {pattern.ranges[i, "seq.name"] <- make.names(names(Biostrings::alignedPattern(pa[i])))}
    pattern.ranges <- tidyr::pivot_longer(pattern.ranges, cols = c(start, end), names_to = "pos", values_to = "values")
    pattern.ranges$position <- ifelse(pattern.ranges$pos == "start", -2, max(algnmt$position) + 2)
    plot <- plot + ggplot2::geom_text(data = pattern.ranges, ggplot2::aes(x = position, y = seq.name, label = values), size = pattern.lim.size, family = font.family, inherit.aes = F)
  }

  if (subject.lim.lines) {
    min.pos <- algnmt %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(position)
    max.pos <- algnmt %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_max(order_by = position, n = 1) %>% dplyr::pull(position)
    plot <- plot + ggplot2::geom_vline(xintercept = c(min.pos, max.pos), linetype = "dashed")
  }

  return(plot)
}

XStringSet_to_df <- function(xstringset) {
  out <- purrr::map(as.list(xstringset), as.character)
  out <- purrr::flatten(purrr::map(out, strsplit, split = ""))
  out <- purrr::map_dfr(out, function(x) utils::stack(stats::setNames(x, seq(1, length(x)))), .id = "seq.name")
  names(out)[c(2,3)] <- c("seq", "position")
  out$position <- as.numeric(as.character(out$position))
  return(out)
}

