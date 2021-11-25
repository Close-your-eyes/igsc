#' Align multiple DNA pattern sequences to one subject
#'
#' This function is useful plot the alignment position of multiple patterns in one subject.
#' It uses Biostrings::pairwiseAlignment(), obtains the individual alignment limits and converts that to a ggplot object with a few optional complementing information.
#' The function will not work if a gap is induced in the subject and at least two patterns overlap at this gap.
#'
#' @param subject a named character or named DNAStringSet of one subject (only the DNAStringSet but not DNAString can hold a name)
#' @param patterns a named character vector or named DNAStringSet of patterns to align to the subject sequence
#' @param type the type of alignment passed to Biostrings::pairwiseAlignment; not every type may work well with this function (if there are overlapping ranges of the alignments to the subject for example)
#' @param pattern.lim.size size of printed limits of aligned patterns (at which nt does the alignment to the subject starts and ends); set to 0 to avoid plotting
#' @param subject.lim.lines print vertical lines at the outermost subject-nts of all aligned patterns
#' @param attach.nt add the length of the string to the name on the axis
#' @param tile.border.color character; tiles from geom_tile are used to plot nts - should they have a border color, e.g. "black"; only useful for short alignment and only an aesthetic thing
#' @param order.patterns order pattern increasingly by alignment position (start)
#' @param perfect.matches.only filter patterns for those which match the subject without gaps, insertions or substitutions before pairwise alignment
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
                                                   pattern.lim.size = 2,
                                                   subject.lim.lines = F,
                                                   attach.nt = T,
                                                   tile.border.color = NA) {

  #alignment.color.palette
  acp <- c("A" = "#ffafaf",
           "T" = "#fed7af",
           "C" = "#afffaf",
           "G" = "#afffff",
           "-" = "#ffffff",
           "match" = "#999999", ##d9d9d9
           "mismatch" = "#a51515",
           "gap" = "#9248d4",
           "insertion" = "#000000",
           "ambiguous" = "#E69F00")

  if (class(patterns) != "DNAStringSet") {
    patterns <- Biostrings::DNAStringSet(patterns)
  }
  if (class(subject) != "DNAStringSet") {
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

  if (perfect.matches.only) {
    perf <- Biostrings::vwhichPDict(subject = subject, pdict = Biostrings::PDict(patterns))[[1]]
    print(paste0(length(perf), " of ", length(patterns), " patterns found to perfectly match the subject."))
    patterns <- patterns[perf]
  }

  # calculate all alignments
  pa <- Biostrings::pairwiseAlignment(subject = subject, pattern = patterns, type = type)

  # find out if any pattern alignment overlap with gaps from another pattern alignment. this would cause problem in the alignment.
  rr <-
    as.data.frame(pa@subject@range) %>%
    dplyr::mutate(alignment = seq(1, nrow(as.data.frame(pa@subject@range)), 1)) %>%
    dplyr::select(-c(width)) %>%
    dplyr::left_join(as.data.frame(pa@subject@indel) %>% dplyr::rename("alignment" = group, "indel.start" = start, "indel.end" = end, "indel.width" = width) %>% dplyr::select(-c(group_name)), by = "alignment") %>%
    tidyr::drop_na() %>%
    dplyr::mutate(indel.start.corr = indel.start + start - 1) %>%
    dplyr::mutate(indel.end.corr = indel.start.corr + indel.width - 1)

  combs <- expand.grid(unique(rr$alignment), unique(rr$alignment))
  combs <- as.matrix(combs[which(combs$Var1 != combs$Var2),])
  if (nrow(combs) > 0) {
    indel.al.overlaps <- unlist(lapply(split(combs, seq(nrow(combs))), function(x) {
      al.range <- unique(rr[which(rr$alignment == x[1]),"start"]):unique(rr[which(rr$alignment == x[1]),"end"])
      unlist(lapply(c(apply(as.matrix(rr[which(rr$alignment == x[2]),c("indel.start.corr", "indel.end.corr")]), 1, function(y) {seq(y[1], y[2], 1)})), function (z) {
        length(intersect(al.range, z))
      }))
    }))
    if (any(indel.al.overlaps > 0)) {
      stop("Overlapping indels and aligment ranges on the subject cannot be handled, yet.")
    }
  }

  # check for indels, just for information
  for (i in seq_along(pa)) {
    if (length(pa[i]@subject@indel@unlistData@start) > 0) {
      print(paste0(pa[i]@pattern@unaligned@ranges@NAMES, " caused ", length(pa[i]@subject@indel@unlistData@start), " indel(s) in the subject."))
    }
  }

  # get ranges
  subject.ranges <- lapply(split(data.frame(pa@subject@range), seq(nrow(data.frame(pa@subject@range)))), function (x) {x$start:x$end})
  subject.ranges.unique <- subject.ranges[which(!duplicated(subject.ranges))]
  pa.unique <- pa[which(!duplicated(subject.ranges))]

  if (length(subject.ranges.unique) > 1) {
    is.subset <- unlist(lapply(seq_along(subject.ranges.unique), function(i) {
      any(unlist(lapply(seq_along(subject.ranges.unique), function (j) {
        if (i == j) {
          return(F)
        } else {
          all(subject.ranges[[i]] %in% subject.ranges[[j]])
        }
      })))
    }))
  } else {
    is.subset <- rep(F, length(subject.ranges))
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

  df <- data.frame(seq = strsplit(total.subject.seq, "")[[1]],
                   position = seq(1:length(strsplit(total.subject.seq, "")[[1]])))

  df[df$seq != "-", "subject.position"] <- seq(1:nrow(df[df$seq != "-", ]))
  names(df)[1] <- subject.name

  pf <- list()
  gap.corr <- 0

  pa <- pa[order(sapply(subject.ranges, function(x) min(x)))]
  for (x in 1:length(pa)) {
    pf[[x]] <- data.frame(seq = strsplit(as.character(Biostrings::alignedPattern(pa[x])), ""),
                          position = (pa[x]@subject@range@start + gap.corr):(pa[x]@subject@range@start+nchar(as.character(Biostrings::alignedPattern(pa[x]))) - 1 + gap.corr))
    gap.corr <- gap.corr + sum(data.frame(pa[x]@subject@indel@unlistData)$width)
  }

  for (i in 1:length(pf)) {
    df <- df %>% dplyr::left_join(pf[[i]], by = "position")
  }

  df.match <- df
  for (x in patterns.names) {
    df.match[,x] <- ifelse(df.match[,x] == df.match[,subject.name], "match", ifelse(df.match[,x] == "-", "-", "mismatch"))
    df.match[,x] <- ifelse(df.match[,x] == "mismatch" & df.match[,subject.name] == "-", "insertion", df.match[,x])
    df.match[,x] <- ifelse(df.match[,x] == "-" & df.match[,subject.name] != "-", "gap", df.match[,x])
  }
  df.match[,subject.name] <- ifelse(df.match[,subject.name] == "-", "gap", df.match[,subject.name])

  df <-
    df %>%
    tidyr::pivot_longer(cols = all_of(c(subject.name, patterns.names)), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names)))
  acp1 <- acp[which(names(acp) %in% unique(df$seq))]
  df$seq <- factor(df$seq, levels = c(names(acp1)))
  df$seq.name <- factor(df$seq.name, levels = c(subject.name, patterns.names))
  if (order.patterns) {
    df$seq.name <- factor(df$seq.name, levels = c(subject.name, patterns.names[order(sapply(subject.ranges, function(x) min(x)))]))
  } else {
    df$seq.name <- factor(df$seq.name, levels = c(subject.name, patterns.names))
  }

  df.match <-
    df.match %>%
    tidyr::pivot_longer(cols = all_of(c(subject.name, patterns.names)), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names)))
  acp2 <- acp[which(names(acp) %in% unique(df.match$seq))]
  df.match$seq <- factor(df.match$seq, levels = c(names(acp2)))
  if (order.patterns) {
    df.match$seq.name <- factor(df.match$seq.name, levels = c(subject.name, patterns.names[order(sapply(subject.ranges, function(x) min(x)))]))
  } else {
    df.match$seq.name <- factor(df.match$seq.name, levels = c(subject.name, patterns.names))
  }

  g1 <- ggplot2::ggplot(df %>% dplyr::filter(!is.na(seq)), ggplot2::aes(x = position, y = seq.name, fill = seq)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::scale_fill_manual(values = acp1) +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::ylab("seq name")


  g2 <- ggplot2::ggplot(df.match %>% dplyr::filter(!is.na(seq)), ggplot2::aes(x = position, y = seq.name, fill = seq)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::scale_fill_manual(values = acp2) +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::ylab("seq name")

  if (!is.na(tile.border.color)) {
    g1 <- g1 + ggplot2::geom_tile(data = df[which(df$seq != "-"), ], color = tile.border.color)
    g2 <- g2 + ggplot2::geom_tile(data = df.match[which(df.match$seq != "-"), ], color = tile.border.color)
  }

  ### pull positions of patterns
  if (pattern.lim.size > 0) {
    pattern.ranges <- data.frame(pa@pattern@range, seq.name = "")
    for (i in 1:length(pa)) {pattern.ranges[i, "seq.name"] <- make.names(names(Biostrings::alignedPattern(pa[i])))}
    pattern.ranges <- pattern.ranges %>% tidyr::pivot_longer(cols = c(start, end), names_to = "pos", values_to = "values")
    pattern.ranges$position <- ifelse(pattern.ranges$pos == "start", -2, max(df.match$position) + 2)
    g2 <- g2 + ggplot2::geom_text(data = pattern.ranges, ggplot2::aes(x = position, y = seq.name, label = values), size = pattern.lim.size, inherit.aes = F)
  }

  min.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(position)
  min.subj.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(subject.position)

  max.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(position)
  max.subj.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(subject.position)

  if (subject.lim.lines) {
    g1 <- g1 + ggplot2::geom_vline(xintercept = min.pos, linetype = "dashed") + ggplot2::geom_vline(xintercept = max.pos, linetype = "dashed")
    g2 <- g2 + ggplot2::geom_vline(xintercept = min.pos, linetype = "dashed") + ggplot2::geom_vline(xintercept = max.pos, linetype = "dashed")
  }

  return(list(base.plot = g1,
              match.plot = g2,
              base.df = df,
              match.df = df.match,
              min.max.subject.position = c(min.subj.pos, max.subj.pos)))
}


integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


