#' Align multiple DNA pattern sequences to one subject
#'
#' This function is useful plot the alignment position of multiple patterns in one subject.
#' It uses Biostrings::pairwiseAlignment(), obtains the individual alignment limits and converts that to a ggplot object with a few optional complementing information.
#' The function will not work if a gap is induced in the subject and at least two patterns overlap at this gap.
#'
#' @param subject a character or DNAStringSet of one subject (to provide a name in final graphics, provide a named DNAStringSet)
#' @param patterns a character vector or DNAStringSet of patterns to align to the subject sequence
#' @param subject.name optional, provide a name for the subject; overwrites names provided with subject
#' @param patterns.names optional, provide names for the patterns; overwrites names provided with patterns
#' @param type the type of alignment passed to Biostrings::pairwiseAlignment; not every type may work well with this function (if there are overlapping ranges of the alignments to the subject for example)
#' @param subject.alignment.limits numeric vector of nt-positions (lower and upper) of the subject to manually limit (cut) the final graphic
#' @param print.pattern.positions check
#' @param pattern.positions.size check
#' @param print.subject.min.max check
#' @param attach.length.to.name add the length of the string to the name on the axis
#' @param tile.border.color character; tiles from geom_tile are used to plot nts - should they have a border color, e.g. "black"; only useful for short alignment and only an aesthetic thing
#' @param title print a title to the alignment
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
#' p <- stats::setNames(c("AAAA", "CCCC", "TTTT", "GGGG", "TTCC"), c("pat1", "pat2", "pat3", "pat4", "pat5"))
#' p <- Biostrings::DNAStringSet(p)
#' als <- MultiplePairwiseAlignmentsToOneSubject(subject = s, patterns = p, type = "local")
#'
#' }
MultiplePairwiseAlignmentsToOneSubject <- function(subject,
                                                   patterns,
                                                   subject.name,
                                                   patterns.names,
                                                   type = "global-local",
                                                   subject.alignment.limits,
                                                   print.pattern.positions = T,
                                                   pattern.positions.size = 2,
                                                   print.subject.min.max = F,
                                                   attach.length.to.name = T,
                                                   tile.border.color = NA,
                                                   title) {

  alignment.colour.palette <- c("A" = "#ffafaf",
                                "T" = "#fed7af",
                                "C" = "#afffaf",
                                "G" = "#afffff",
                                "-" = "#ffffff",
                                "match" = "#999999", ##d9d9d9
                                "mismatch" = "#a51515",
                                "gap" = "#9248d4",
                                "insertion" = "#000000",
                                "ambiguous" = "#E69F00")

  # fix missing
  if (missing(subject.name)) {
    if (is.null(names(subject))) {
      subject.name <- "subject"
      names(subject) <- subject.name
    } else {
      subject.name <- names(subject)
    }
  }

  if (missing(patterns.names)) {
    if (is.null(names(patterns))) {
      patterns.names <- paste0("pattern ", seq(1:length(patterns)))
      names(patterns) <- patterns.names
    } else {
      patterns.names <- names(patterns)
    }
  } else {
    if (length(patterns.names) != length(patterns)) {
      stop("patterns.names must have same length as patterns.")
    }
    if (length(unique(patterns.names)) != length(patterns.names)) {
      stop("patterns.names must be unique.")
    }
  }

  if (attach.length.to.name) {
    patterns.names <- paste0(patterns.names, "_", sapply(as.character(patterns), function(x) nchar(x)), "nt")
    subject.name <- paste0(subject.name, "_", nchar(as.character(subject)), "nt")
  }

  patterns.names <- make.names(patterns.names)
  subject.name <- make.names(subject.name)

  if (class(patterns) != "DNAStringSet") {
    patterns <- Biostrings::DNAStringSet(patterns)
  }
  if (class(subject) != "DNAStringSet") {
    subject <- Biostrings::DNAStringSet(subject)
  }

  names(patterns) <- patterns.names
  names(subject) <- subject.name

  type <- match.arg(type, choices = c("global", "local", "overlap", "global-local", "local-global"))

  # alignment limits
  if (!missing(subject.alignment.limits)) {
    if (length(subject.alignment.limits) != 2) {
      return("Please provide an upper and lower limit for the aligment.")
    }
    if (subject.alignment.limits[1] > subject.alignment.limits[2]) {
      subject.alignment.limits <- rev(subject.alignment.limits)
    }
    subject <- stringr::str_sub(subject, subject.alignment.limits[1], subject.alignment.limits[2])
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
  lapply(pa, function(x) {
    if (length(x@subject@indel@unlistData@start) > 0) {
      print(paste0(x@pattern@unaligned@ranges@NAMES, " caused ", length(x@subject@indel@unlistData@start), " indel(s) in the subject."))
    }
    return(NULL)
  })

  # get ranges
  subject.ranges <- lapply(split(data.frame(pa@subject@range), seq(nrow(data.frame(pa@subject@range)))), function (x) {x$start:x$end})

  ## error here, make new variable with
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
  for (i in rev(rev(seq_along(subject.ranges.unique))[-1])) {
    # test if there is a gap between the ith and the i+1th alignment; if so, fill with original sequence
    if (max(subject.ranges.unique[[i]]) < min(subject.ranges.unique[[i+1]])) {
      # if yes use original subject
      total.subject.seq <- paste0(total.subject.seq, substr(subject, max(subject.ranges.unique[[i]])+1, min(subject.ranges.unique[[i+1]])-1))
    } else {
      # if not use seq from alignment
      r <- min(subject.ranges.unique[[i]])
      total.subject.seq <- paste0(total.subject.seq, substr(pa.unique[i]@subject, min(subject.ranges.unique[[i]])-r+1, min(subject.ranges.unique[[i+1]])-r))
    }
    #print(total.subject.seq)
  }
  r <- min(subject.ranges.unique[[length(subject.ranges.unique)]])
  total.subject.seq <- paste0(total.subject.seq, substr(pa.unique[length(subject.ranges.unique)]@subject, 1, max(subject.ranges.unique[[length(subject.ranges.unique)]])+1-r))
  ## attach the remaining sequence from subject
  total.subject.seq <- paste0(total.subject.seq, substr(subject, max(subject.ranges.unique[[length(subject.ranges.unique)]])+1, nchar(as.character(subject))))

  df <- data.frame(seq = strsplit(total.subject.seq, "")[[1]],
                   position = seq(1:length(strsplit(total.subject.seq, "")[[1]])))

  df[df$seq != "-", "subject.position"] <- seq(1:nrow(df[df$seq != "-", ]))
  names(df)[1] <- subject.name

  pf <- list()
  gap.corr <- 0
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
    dplyr::mutate(seq = ifelse(is.na(seq), "-", seq)) %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names))) %>%
    dplyr::mutate(seq = factor(seq, levels = c(names(alignment.colour.palette)[1:4], names(alignment.colour.palette)[6:9], names(alignment.colour.palette)[5])))

  df.match <-
    df.match %>%
    tidyr::pivot_longer(cols = all_of(c(subject.name, patterns.names)), names_to = "seq.name", values_to = "seq") %>%
    dplyr::mutate(seq = ifelse(is.na(seq), "-", seq)) %>%
    dplyr::mutate(seq.name = factor(seq.name, levels = c(subject.name, patterns.names))) %>%
    dplyr::mutate(seq = factor(seq, levels = c(names(alignment.colour.palette)[1:4], names(alignment.colour.palette)[6:9], names(alignment.colour.palette)[5])))

  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = seq.name, fill = seq)) +
    ggplot2::geom_tile(color = tile.border.color) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = alignment.colour.palette)


  g2 <- ggplot2::ggplot(df.match, ggplot2::aes(x = position, y = seq.name, fill = seq)) +
    ggplot2::geom_tile(color = tile.border.color) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = alignment.colour.palette)

  ### pull positions of patterns
  if (print.pattern.positions) {
    pattern.ranges <- data.frame(pa@pattern@range, seq.name = "")
    for (i in 1:length(pa)) {pattern.ranges[i, "seq.name"] <- make.names(names(Biostrings::alignedPattern(pa[i])))}
    pattern.ranges <- pattern.ranges %>% tidyr::pivot_longer(cols = c(start, end), names_to = "pos", values_to = "values")
    pattern.ranges$position <- ifelse(pattern.ranges$pos == "start", -2, max(df.match$position) + 2)
    g2 <- g2 + ggplot2::geom_text(data = pattern.ranges, ggplot2::aes(x = position, y = seq.name, label = values), size = pattern.positions.size, inherit.aes = F)
  }

  min.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(position)
  min.subj.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = position, n = 1) %>% dplyr::pull(subject.position)

  max.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(position)
  max.subj.pos <- df %>% dplyr::filter(seq != "-") %>% dplyr::filter(seq.name != names(subject)) %>% dplyr::slice_min(order_by = -position, n = 1) %>% dplyr::pull(subject.position)

  if (print.subject.min.max) {
    g2 <- g2 + ggplot2::geom_vline(xintercept = min.pos, linetype = "dashed") + ggplot2::geom_vline(xintercept = max.pos, linetype = "dashed")
  }

  if (!missing(title)) {
    g1 <- g1 + ggplot2::ggtitle(title)
    g2 <- g2 + ggplot2::ggtitle(title)
  }

  return(list(base.plot = g1, match.plot = g2, base.df = df, match.df = df.match, min.max.subject.position = c(min.subj.pos, max.subj.pos)))
}

