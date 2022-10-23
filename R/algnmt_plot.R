#' Plot sequence alignments as ggplot object
#'
#' @param algnmt an alignment object returned from (i) igsc::MultiplePairwiseAlignmentsToOneSubject, (ii) DECIPHER::AlignSeqs or (iii) Biostrings::pairwiseAlignment;
#' in case (i) this is a data.frame, in case (ii) this is a XStringSet, in case (iii) this is a PairwiseAlignmentsSingleSubject;
#' alternatively provide a custom data frame which has at least to contain pos_col, seq_col, name_col (see other arguments)
#' @param color_values a color scale for NTs or AAs (depending on type of algnmt); provide a named vector of colors where names are NTs or AAs;
#' or choose a name from igsc:::scheme_NT or igsc:::scheme_AA; or leave NULL to have a scheme choosen automatically
#' @param tile.border.color a color to draw tile borders with; e.g. "black" or a hex code;
#' leave NA to have no borders (quicker plotting and advised for long alignments)
#' @param tile.border.on.NA logical whether to draw borders on NA positions
#' @param text logical whether to draw text in tiles (NT or AA identifier);
#' only advisable for short alignments
#' @param font.family which font to use for plotting
#' @param theme ggplot theme to use as basis
#' @param pattern.lim.size numeric; plot annotation of pattern alignment limits; value indicates size;
#' set to 0 to omit plotting; only applicable if pa is provided
#' @param pa list of pairwise alignment; rather for internal use by igsc::MultiplePairwiseAlignmentsToOneSubject;
#' from pa pattern alignment limits are derived
#' @param subject.lim.lines logical whether to plot vertical lines of subject alignment limits
#' @param pos_col name of position column in algnmt (applicable if algnmt is a data.frame)
#' @param seq_col name of sequence column in algnmt (applicable if algnmt is a data.frame)
#' @param name_col name of column which holds sequence names in algnmt (applicable if algnmt is a data.frame)
#' @param coord_fixed_ratio numeric; fixed aspect ratio of plot; leave NULL to not force a ratio
#' @param x.breaks numeric vector; manually provide breaks (ticks) on x-axis; leave NULL to have breaks
#' picked automatically based on data provided
#' @param algnmt_type 'NT' or 'AA'; only required if algnmt is a data.frame and only if
#' NT or AA cannot be guessed; leave NULL to have it guessed based on data
#' @return ggplot2 object of alignment
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
algnmt_plot <- function(algnmt,
                        color_values = NULL,
                        tile.border.color = NA,
                        tile.border.on.NA = F,
                        text = F,
                        font.family = c("sans", "mono", "serif"),
                        theme = ggplot2::theme_bw(base_family = font.family),
                        pattern.lim.size = 2,
                        pa = NULL,
                        subject.lim.lines = F,
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        coord_fixed_ratio = NULL,
                        x.breaks = NULL,
                        algnmt_type = NULL) {

  # example: length(unique(algnmt[,pos_col,drop=T]))/length(unique(algnmt[,name_col,drop=T]))*0.5

  font.family <- match.arg(font.family, choices = c("sans", "mono", "serif"))

  ## checks
  if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet") || methods::is(algnmt, "AAStringSet")) {
    # from DECIPHER::AlignSeqs
    if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet")) {
      algnmt_type <- "NT"
    }
    if (methods::is(algnmt, "AAStringSet")) {
      algnmt_type <- "AA"
    }
    algnmt <- XStringSet_to_df(algnmt)
  } else if (methods::is(algnmt, "PairwiseAlignmentsSingleSubject")) {
    # from Biostrings::pairwiseAlignment
    if (methods::is(pa@pattern@unaligned, "QualityScaledDNAStringSet") || methods::is(pa@pattern@unaligned, "QualityScaledRNAStringSet")) {
      algnmt_type <- "NT"
    }
    if (methods::is(pa@pattern@unaligned, "QualityScaledAAStringSet")) {
      algnmt_type <- "AA"
    }
    algnmt <- pa_to_df(algnmt)
  }
  # else: data frame from igsc::MultiplePairwiseAlignmentsToOneSubject

  if (!all(c(pos_col, seq_col, name_col) %in% names(algnmt))) {
    stop("algnmt at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ".Alternative change function arguments.")
  }

  # if seq names are numeric; order them increasingly
  if (!anyNA(suppressWarnings(as.numeric(algnmt$seq.name)))) {
    algnmt$seq.name <- factor(algnmt$seq.name, levels = as.character(unique(as.numeric(as.character(algnmt$seq.name)))))
  } else {
    algnmt$seq.name <- as.factor(algnmt$seq.name)
  }


  if (is.null(algnmt_type)) {
    inds <- intersect(which(!is.na(algnmt[,seq_col,drop=T])), which(!algnmt[,seq_col,drop=T] %in% c("match", "mismatch", "gap", "insertion", "ambiguous")))
    algnmt_type <- guess_type(seq_vector = algnmt[,seq_col,drop=T][inds])
  }

  # use preset colors
  if (is.null(color_values)) {
    if (algnmt_type == "NT") {
      color_values <- colnames(scheme_NT)[1]
    }
    if (algnmt_type == "AA") {
      color_values <- colnames(scheme_AA)[1]
    }
  }
  if (length(color_values) == 1) {
    if (algnmt_type == "NT") {
      color_values <- scheme_NT[,match.arg(color_values, choices = colnames(scheme_NT)),drop=T]
    } else if (algnmt_type == "AA") {
      if (color_values == "Chemistry_AA") {
        # special case; change legend to chemical property
        col_col <- "chem prop"
        chem_col <- utils::stack(aa_info[["aa_main_prop"]])
        names(chem_col) <- c(col_col, "seq")
        chem_col$seq <- as.character(chem_col$seq)
        algnmt <- dplyr::left_join(algnmt, chem_col, by = "seq")
        algnmt[,col_col] <- ifelse(is.na(algnmt[,col_col,drop=T]), algnmt$seq, algnmt[,col_col,drop=T])
        color_values <- scheme_AA[,match.arg(color_values, choices = colnames(scheme_AA)),drop=T]
        color_values <- color_values[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
        for (i in 1:length(color_values)) {
          if (names(color_values)[i] %in% names(aa_info[["aa_main_prop"]])) {
            names(color_values)[i] <- aa_info[["aa_main_prop"]][names(color_values)[i]]
          }
        }
      } else {
        col_col <- seq_col
        color_values <- scheme_AA[,match.arg(color_values, choices = colnames(scheme_AA)),drop=T]
      }
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

  if (!is.null(names(color_values)) && col_col == seq_col) {
    color_values <- color_values[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
  }

  # make sure it is not a factor
  if (!is.numeric(algnmt[,pos_col,drop=T])) {
    algnmt[,pos_col] <- as.numeric(as.character(algnmt[,pos_col,drop=T]))
  }
  algnmt[,seq_col] <- as.character(algnmt[,seq_col,drop=T])

  # https://stackoverflow.com/questions/45493163/ggplot-remove-na-factor-level-in-legend
  # --> na.translate = F

  if (is.null(x.breaks)) {
    x.breaks <- floor(base::pretty(algnmt[,pos_col,drop=T], n = 5))
    x.breaks <- x.breaks[which(x.breaks > 0)]
    x.breaks <- x.breaks[which(x.breaks < max(algnmt[,pos_col,drop=T]))]
  }

  plot <-
    ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(name_col))) +
    theme +
    ggplot2::theme(panel.grid = element_blank(),
                   text = ggplot2::element_text(family = font.family)) +
    ggplot2::scale_fill_manual(values = color_values, na.value = "white", na.translate = F) +
    ggplot2::scale_x_continuous(breaks = x.breaks)

  if (!is.na(tile.border.color)) {
    plot <- plot + ggplot2::geom_tile(data = if(tile.border.on.NA) {algnmt} else {algnmt[which(!is.na(algnmt[,seq_col,drop=T])), ]},
                                      color = tile.border.color,
                                      ggplot2::aes(fill = !!rlang::sym(col_col)))
  } else {
    plot <- plot + ggplot2::geom_tile(ggplot2::aes(fill = !!rlang::sym(col_col)))
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
    pattern.ranges$position <- ifelse(pattern.ranges$pos == "start", -2, max(algnmt[,pos_col,drop=T]) + 2)
    plot <- plot + ggplot2::geom_text(data = pattern.ranges, ggplot2::aes(x = position, y = seq.name, label = values), size = pattern.lim.size, family = font.family, inherit.aes = F)
  }

  if (subject.lim.lines) {
    min.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!!rlang::sym(name_col) != names(subject)) %>% dplyr::slice_min(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
    max.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!!rlang::sym(name_col) != names(subject)) %>% dplyr::slice_max(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
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

  # do this in parent fun
  #if (!anyNA(suppressWarnings(as.numeric(out$seq.name)))) {out$seq.name <- factor(out$seq.name, levels = as.character(unique(as.numeric(as.character(out$seq.name)))))}
  return(out)
}

pa_to_df <- function(pa) {
  if (length(pa) > 1) {
    stop("Please provide one pairwiseAlignment only (= one pattern only). length(algnmt) should be 1.")
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledDNAStringSet")) {
    xstringfun <- Biostrings::DNAStringSet
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledRNAStringSet")) {
    xstringfun <- Biostrings::RNAStringSet
  }
  if (methods::is(pa@pattern@unaligned, "QualityScaledAAStringSet")) {
    xstringfun <- Biostrings::AAStringSet
  }

  if (any(c(is.null(pa@pattern@unaligned@ranges@NAMES), is.null(pa@subject@unaligned@ranges@NAMES)))) {
    message("No names provided for pattern and/or subject. If desired provide named XStringSets, each of length 1, to Biostrings::pairwiseAlignment.")
  }
  pattern_name <- ifelse(is.null(pa@pattern@unaligned@ranges@NAMES), "pattern", pa@pattern@unaligned@ranges@NAMES)
  subject_name <- ifelse(is.null(pa@subject@unaligned@ranges@NAMES), "subject", pa@subject@unaligned@ranges@NAMES)
  return(xstringfun(stats::setNames(c(as.character(pa@pattern), as.character(pa@subject)),
                                    nm = c(pattern_name, subject_name))))
}

.integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    return(breaks[intersect(which(breaks > 0), which(breaks < max(x)))])
  }
  return(fxn)
}

guess_type <- function(seq_vector) {
  # current not so strict because even c("A", "T", "G", "C") could not be guessed
  if (all(unique(seq_vector) %in% unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N")))) {
    # && !all(seq_vector %in% c(Biostrings::AA_ALPHABET, "N"))
    return("NT")
  } else if (all(seq_vector %in% c(Biostrings::AA_ALPHABET, "N"))) {
    # !all(unique(seq_vector) %in% unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N"))) &&
    return("AA")
  } else {
    stop("Alignment type could not be determined unambigously. Please define it: 'NT' or 'AA'.")
  }
}
