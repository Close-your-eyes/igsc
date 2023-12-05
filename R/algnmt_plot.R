#' Plot sequence alignments as ggplot object
#'
#' @param algnmt an alignment object returned from (i) igsc::MultiplePairwiseAlignmentsToOneSubject, (ii) DECIPHER::AlignSeqs or (iii) Biostrings::pairwiseAlignment;
#' in case (i) this is a data.frame, in case (ii) this is a XStringSet, in case (iii) this is a PairwiseAlignmentsSingleSubject;
#' alternatively provide a custom data frame which has at least to contain pos_col, seq_col, name_col (see other arguments)
#' @param color_values a color scale for NTs or AAs (depending on type of algnmt); provide a named vector of colors where names are NTs or AAs;
#' or choose a name from igsc:::scheme_NT or igsc:::scheme_AA; or choose one from names(purrr::flatten(Peptides:::AAdata)); or leave NULL for the default scheme
#' @param tile.border.color a color to draw tile borders with; e.g. "black" or a hex code;
#' leave NA to have no borders (quicker plotting and advised for long alignments)
#' @param tile.border.on.NA logical whether to draw borders on NA positions
#' @param text logical whether to draw text in tiles (NT or AA identifier); or numeric indicating text size for geom_text()
#' advisable only for short alignments
#' @param font.family which font to use for plotting
#' @param theme ggplot theme to use as basis
#' @param pattern.lim.size numeric; plot annotation of pattern alignment limits; value indicates size;
#' set to 0 to omit plotting; only applicable if pa is provided
#' @param subject.lim.lines logical whether to plot vertical lines of subject alignment limits
#' @param pos_col name of position column in algnmt (applicable if algnmt is a data.frame)
#' @param seq_col name of sequence column in algnmt (applicable if algnmt is a data.frame)
#' @param name_col name of column which holds sequence names in algnmt (applicable if algnmt is a data.frame)
#' @param coord_fixed_ratio numeric; fixed aspect ratio of plot; leave NULL to not force a ratio
#' @param x.breaks numeric vector; manually provide breaks (ticks) on x-axis; leave NULL to have breaks
#' picked automatically based on data provided
#' @param algnmt_type 'NT' or 'AA'; only required if algnmt is a data.frame and only if
#' NT or AA cannot be guessed; leave NULL to have it guessed based on data
#' @param ref a reference sequence to compare all other sequnces to; e.g. a consensus sequence
#' @param subject_name name of the subject sequence in algnmt; only needed if subject.lim.lines = TRUE
#' @param pattern.lim.pos
#' @param plot.pattern.names
#' @param pairwiseAlignment
#' @param y_group_col
#' @param plot.pattern.names_fun
#' @param add_length_suffix
#' @param group_on_yaxis
#' @param min_gap
#' made with DECIPHER::ConsensusSequence; if provided matching residues are replaced by a dot (.)
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
                        pattern.lim.pos = "inner",
                        plot.pattern.names = F,
                        plot.pattern.names_fun = ggrepel::geom_text_repel, # ggplot2::geom_text
                        pairwiseAlignment = NULL,
                        subject.lim.lines = F,
                        subject_name = NULL,
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        y_group_col = NULL, # set groups of patterns to plot on same y-axis level
                        coord_fixed_ratio = NULL,
                        x.breaks = NULL,
                        algnmt_type = NULL,
                        add_length_suffix = F,
                        group_on_yaxis = F,
                        min_gap = 100,
                        ref = NULL) {
  ## check duplicate names
  ## y_group_col actually makes only sense when a algnmt is a dataframe, e.g. from multi pairwisealignment
  ## otherwise y_group_col is set to NULL

if (!requireNamespace("Peptides", quietly = T)){
  utils::install.packages("Peptides")
}
if (!requireNamespace("farver", quietly = T)){
  utils::install.packages("farver")
}
if (!requireNamespace("viridisLite", quietly = T)){
  utils::install.packages("viridisLite")
}

if (subject.lim.lines && is.null(subject_name)) {
  message("When subject.lim.lines = T, subject_name cannot be NULL. subject.lim.lines set to FALSE")
  subject.lim.lines <- F
}

if (subject.lim.lines && !subject_name %in% algnmt[,name_col,drop=T]) {
  stop("subject_name not found in ", name_col, " of algnmt.")
}

font.family <- match.arg(font.family, choices = c("sans", "mono", "serif"))
pattern.lim.pos <- match.arg(pattern.lim.pos, c("inner", "outer"))
plot.pattern.names_fun <- match.fun(plot.pattern.names_fun)

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
  if (!is.null(y_group_col)) {
    message("y_group_col is set to NULL.")
    y_group_col <- NULL
  }
} else if (methods::is(algnmt, "PairwiseAlignmentsSingleSubject")) {
  # from Biostrings::pairwiseAlignment
  if (methods::is(pairwiseAlignment, "PairwiseAlignmentsSingleSubject")) {
    temp <- pairwiseAlignment[1]@pattern@unaligned
  } else if (methods::is(pairwiseAlignment, "list")) {
    temp <- pairwiseAlignment[[1]]@pattern@unaligned
  }
  if (methods::is(temp, "QualityScaledDNAStringSet") || methods::is(temp, "QualityScaledRNAStringSet")) {
    algnmt_type <- "NT"
  }
  if (methods::is(temp, "QualityScaledAAStringSet")) {
    algnmt_type <- "AA"
  }

  algnmt <- pa_to_df(pa = algnmt)
  if (!is.null(y_group_col)) {
    message("y_group_col is set to NULL.")
    y_group_col <- NULL
  }
}
# else: data frame from igsc::MultiplePairwiseAlignmentsToOneSubject

if (group_on_yaxis && is.null(y_group_col)) {
  # if pairwise alignment is provided, use this one to derived ranges: [TODO]
  if (!is.null(pairwiseAlignment)) {
    subject.ranges <- seq2(pairwiseAlignment@subject@range@start, pairwiseAlignment@subject@range@start+pairwiseAlignment@subject@range@width-1)
    names(subject.ranges) <- pairwiseAlignment@pattern@unaligned@ranges@NAMES
    subject_name <- setdiff(unique(algnmt[,name_col,drop=T]), names(subject.ranges))
  } else {
    # here the subject has to be separated from patterns
    subject.ranges <- seq2(algnmt %>% tidyr::drop_na() %>% dplyr::group_by(!!rlang::sym(name_col)) %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col)),
                           algnmt %>% tidyr::drop_na() %>% dplyr::group_by(!!rlang::sym(name_col)) %>% dplyr::slice_max(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col)))
    names(subject.ranges) <- as.character(algnmt %>% tidyr::drop_na() %>% dplyr::group_by(!!rlang::sym(name_col)) %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(name_col)))
    if (is.null(subject_name)) {
      max_len <- max(lengths(subject.ranges))
      if (length(which(lengths(subject.ranges) == max_len)) > 1) {
        message("Subject could not be identified as there are min. 2 sequences which have the max length. Cannot order patterns on y-axis.")
        ## TODO: set variable for ordering pattern on y-axis to FALSE here
      } else {
        subject_name <- names(which(lengths(subject.ranges) == max_len))
      }
    }
    subject.ranges <- subject.ranges[which(names(subject.ranges) != subject_name)]
  }

  subject.ranges_ordered <- subject.ranges[order(sapply(subject.ranges, min))]
  rows <- numeric(length(subject.ranges_ordered))
  rows[1] <- 1
  # priority to row 1, or lowest row in general
  for (i in seq_along(subject.ranges_ordered)[-1]) {
    row_set <- 1
    j <- max(which(rows[1:(i-1)] == row_set))
    overlap <- T
    while (overlap) {
      # loop through all pattern of that row and check if there a larger overlap than allowed (or gap smaller than allowed)
      # what if min_gap are negative values? and what if min_gap is larger than the current range - select the larger of range or min_gap
      # as subject.ranges are ordered: check for every pattern or just the most recent one? - max(j) does that
      if (any(sapply(j, function(x) length(intersect(c(subject.ranges_ordered[[x]], (max(subject.ranges_ordered[[x]])+1):(max(subject.ranges_ordered[[x]])+1+min_gap)), subject.ranges_ordered[[i]])) != 0))) {
        row_set <- row_set + 1
        j <- which(rows == row_set) # which patterns are in next row to consider
        if (length(j) == 0) {
          # first pattern in that row - no further checking needed
          overlap <- F
        } else {
          j <- max(j)
        }
      } else {
        overlap <- F
      }
    }
    rows[i] <- row_set
  }
  rows <- paste0("Group_", rows)
  names(rows) <- names(subject.ranges_ordered)
  rows <- c(stats::setNames(subject_name, subject_name), rows)
  y_group_col <- "group"
  algnmt$group <- rows[as.character(unname(algnmt[,name_col,drop=T]))]
  algnmt$group <- factor(algnmt$group, levels = c(subject_name, unique(rows[-1])))

} else if (group_on_yaxis && !is.null(y_group_col)) {
  message("y_group_col is not NULL. Using this one. Ignoring group_on_yaxis.")
}


if (!is.null(ref)) {
  if (!ref %in% unique(algnmt[,name_col,drop=T])) {
    stop("ref not found in names of algnmt.")
  }
  seq_original <- "seq_original"
  algnmt <- compare_seqs_df(df = algnmt,
                            ref = ref,
                            pos_col = pos_col,
                            seq_col = seq_col,
                            name_col = name_col,
                            seq_original = seq_original)

  # seq_original will cause that tiles are filled according to chemical property even if they are replaced by a dots when
  # sequence is equal to the reference sequence (ref)
} else {
  seq_original <- seq_col
}

if (!all(c(pos_col, seq_col, name_col, y_group_col) %in% names(algnmt))) {
  stop("algnmt at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ". Alternatively change function arguments.")
}

# if seq names are numeric; order them increasingly
if (!anyNA(suppressWarnings(as.numeric(as.character(algnmt[,name_col,drop=T]))))) {
  algnmt[,name_col] <- factor(algnmt[,name_col,drop=T], levels = as.character(unique(as.numeric(as.character(algnmt[,name_col,drop=T])))))
} else {
  algnmt[,name_col] <- as.factor(algnmt[,name_col,drop=T])
}

if (is.null(algnmt_type)) {
  inds <- intersect(which(!is.na(algnmt[,seq_col,drop=T])), which(!algnmt[,seq_col,drop=T] %in% c("match", "mismatch", "gap", "insertion", "ambiguous")))
  algnmt_type <- guess_type(seq_vector = algnmt[,seq_col,drop=T][inds])
}


# use preset colors
if (is.null(color_values)) {
  if (algnmt_type == "NT") {
    #tile_color <- colnames(igsc:::scheme_NT)[1]
    color_values <- colnames(igsc:::scheme_NT)[1]
  }
  if (algnmt_type == "AA") {
    color_values <- colnames(igsc:::scheme_AA)[1]
    # set color_values to Chemistry_AA by default which will cause falling into the special case below
    # to avoid that, change this here to tile_color (like in case of algnmt_type == "NT") and connect
  }
}

#default; may change when algnmt_type == "AA" && color_values == "Chemistry_AA"
col_col <- seq_col

if (length(color_values) == 1 && color_values %in% c(colnames(igsc:::scheme_AA),
                                                     colnames(igsc:::scheme_NT),
                                                     names(purrr::flatten(Peptides:::AAdata)))) {
  if (algnmt_type == "NT") {
    tile_color <- igsc:::scheme_NT[,match.arg(color_values, choices = colnames(igsc:::scheme_NT)),drop=T]
  } else if (algnmt_type == "AA") {
    if (color_values == "Chemistry_AA") {
      # special case; change legend to chemical property
      col_col <- color_values
      chem_col <- utils::stack(igsc:::aa_info[["aa_main_prop"]])
      names(chem_col) <- c(col_col, seq_original)
      chem_col[,seq_original] <- as.character(chem_col[,seq_original,drop=T])
      algnmt <- dplyr::left_join(algnmt, chem_col, by = seq_original)
      algnmt[,col_col] <- ifelse(is.na(algnmt[,col_col,drop=T]), algnmt[,seq_original,drop=T], algnmt[,col_col,drop=T])
      tile_color <- igsc:::scheme_AA[,match.arg(color_values, choices = colnames(igsc:::scheme_AA)),drop=T]
      tile_color <- tile_color[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
      for (i in 1:length(tile_color)) {
        if (names(tile_color)[i] %in% names(igsc:::aa_info[["aa_main_prop"]])) {
          names(tile_color)[i] <- igsc:::aa_info[["aa_main_prop"]][names(tile_color)[i]]
        }
      }
      tile_color <- tile_color[unique(names(tile_color))]
    } else if (color_values %in% names(purrr::flatten(Peptides:::AAdata))) {
      col_col <- color_values
      algnmt[,col_col] <- purrr::flatten(Peptides:::AAdata)[[color_values]][algnmt[,seq_original,drop=T]]
      tile_color <- viridisLite::viridis(100)[cut(unique(algnmt[,col_col,drop=T])[which(!is.na(unique(algnmt[,col_col,drop=T])))], 100)]
      names(tile_color) <- unique(algnmt[,col_col,drop=T])[which(!is.na(unique(algnmt[,col_col,drop=T])))]
    } else {
      tile_color <- igsc:::scheme_AA[,match.arg(color_values, choices = colnames(igsc:::scheme_AA)),drop=T]
    }
  } else {
    message("Type of alignment data (NT or AA) could not be determined. Choosing default ggplot colors.")
    tile_color <- scales::hue_pal()(length(unique(as.character(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))))
  }
} else {
  # provide own colors
  if (is.null(names(color_values))) {
    if (length(color_values) < length(unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))) {
      warning("Less colors provided than entities in the alignment.")
    }
  } else {
    if (any(!unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]) %in% names(color_values))) {
      warning("Not all values in algnmt[,seq_col] found in names(color_values).")
    }
  }
  tile_color <- color_values
}

if (!is.null(names(tile_color)) && col_col == seq_col) {
  tile_color <- tile_color[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
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

if (algnmt_type == "AA" && length(color_values) == 1 && color_values == "Chemistry_AA") {
  col_breaks = unique(igsc:::aa_info[["aa_main_prop"]])[which(unique(igsc:::aa_info[["aa_main_prop"]]) != "stop")]
} else {
  col_breaks = ggplot2::waiver()
}


if (add_length_suffix) {
  nts1 <- stats::setNames(c(pairwiseAlignment@subject@unaligned@ranges@width,
                            pairwiseAlignment@pattern@unaligned@ranges@width),
                          c(pairwiseAlignment@subject@unaligned@ranges@NAMES,
                            pairwiseAlignment@pattern@unaligned@ranges@NAMES))

  nts2 <- paste0(names(nts1), "\n", nts1, " ", tolower(algnmt_type))
  names(nts2) <- names(nts1)
}

yaxis <- ifelse(is.null(y_group_col), rlang::sym(name_col), rlang::sym(y_group_col))
plot <-
  ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!yaxis)) + # name_col
  theme +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), text = ggplot2::element_text(family = font.family)) +
  ggplot2::scale_x_continuous(breaks = x.breaks)

if (add_length_suffix && is.null(y_group_col)) {
  # nt can only be added to y-axis text when sequences are not grouped
  plot <- plot + ggplot2::scale_y_discrete(labels = nts2[levels(algnmt$seq.name)])
}


if (is.numeric(algnmt[,col_col,drop=T])) {
  plot <- plot + ggplot2::scale_fill_viridis_c(na.value = "white")
} else {
  plot <- plot + ggplot2::scale_fill_manual(values = tile_color, breaks = col_breaks, na.value = "white", na.translate = F)
}

if (!is.na(tile.border.color)) {
  plot <- plot + ggplot2::geom_tile(data = if(tile.border.on.NA) {algnmt} else {algnmt[which(!is.na(algnmt[,seq_col,drop=T])),]},
                                    color = tile.border.color,
                                    ggplot2::aes(fill = !!rlang::sym(col_col)))
  # geom_raster seems not take color?!
} else {
  plot <- plot + ggplot2::geom_tile(data = algnmt[which(!is.na(algnmt[,seq_col,drop=T])),],
                                    ggplot2::aes(fill = !!rlang::sym(col_col)))
  # geom_raster # geom_tile
}

if ((is.logical(text) && text) || (is.numeric(text) && text > 0)) {
  if (is.logical(text)) {
    text_size <- 4
  } else {
    text_size <- text
  }

  bckgr_colors <- farver::decode_colour(tile_color[as.character(algnmt[,col_col,drop=T])], to = "hcl")
  algnmt$text_colors <- ifelse(bckgr_colors[, "l"] > 50, "black", "white")

  plot <-
    plot +
    ggplot2::geom_text(data = algnmt, ggplot2::aes(label = !!rlang::sym(seq_col), color = I(text_colors)),
                       na.rm = T, family = font.family, size = text_size)

  #plot <- plot + ggplot2::geom_text(ggplot2::aes(label = !!rlang::sym(seq_col), color = farver::decode_colour(tile_color[as.character(algnmt[,col_col,drop=T])], to = "hcl")[,"l"] > 50), na.rm = T, family = font.family, size = text_size) + ggplot2::scale_color_manual(guide = "none", values = c("black", "white"))
}

if (!is.null(coord_fixed_ratio)) {
  plot <- plot + ggplot2::coord_fixed(ratio = coord_fixed_ratio)
}

if (pattern.lim.size > 0 && !is.null(pairwiseAlignment)) {
  subject.ranges <- as.data.frame(pairwiseAlignment@subject@range)[-3]
  names(subject.ranges)[1:2] <- paste0("subject_", names(subject.ranges)[1:2])
  pattern.ranges <- data.frame(pairwiseAlignment@pattern@range,
                               seq.name = ifelse(rep(is.null(pairwiseAlignment@pattern@unaligned@ranges@NAMES), length(pairwiseAlignment)),
                                                 paste0("pattern_", seq(1,length(pairwiseAlignment))),
                                                 pairwiseAlignment@pattern@unaligned@ranges@NAMES))
  names(pattern.ranges)[1:2] <- paste0("pattern_", names(pattern.ranges)[1:2])
  ranges <- cbind(pattern.ranges, subject.ranges)
  ranges$subject_mid <- ranges$subject_start + (c(ranges$subject_end - ranges$subject_start)/2)
  ranges <- tidyr::pivot_longer(ranges, cols = c(subject_start, subject_end), names_to = "pos", values_to = "inner_position")
  ranges$outer_position <- ifelse(ranges$pos == "subject_start", -2, max(algnmt[,pos_col,drop=T]) + 2)
  ranges$label <- ifelse(ranges$pos == "subject_start", ranges$pattern_start, ranges$pattern_end)
  # maybe adjust style of labels
  label_plot_pos <- ifelse(pattern.lim.pos == "inner", rlang::sym("inner_position"), rlang::sym("outer_position"))

  if (!is.null(y_group_col)) {
    unique_group_seq <- unique(algnmt[,c(y_group_col, "seq.name")])
    unique_group_seq <- stats::setNames(unique_group_seq[,y_group_col,drop=T], unique_group_seq$seq.name)
    ranges[,y_group_col] <- unique_group_seq[ranges$seq.name]
  }
  # nudge labels alternatingly up and down?
  plot <- plot + ggplot2::geom_text(data = ranges,
                                    ggplot2::aes(x = !!label_plot_pos, y = !!yaxis, label = label),
                                    size = pattern.lim.size,
                                    family = font.family,
                                    inherit.aes = F)
  if (plot.pattern.names) {
    ## option to plot this independent of pattern.lims
    ## ggrepel considers the pattern.lim labels - great.
    ranges_sub <- dplyr::distinct(ranges, !!yaxis, subject_mid, seq.name)
    if (add_length_suffix) {
      ranges_sub$label <- nts2[ranges_sub$seq.name]
    } else {
      ranges_sub$label <- ranges_sub$seq.name
    }

    plot <- plot + plot.pattern.names_fun(data = ranges_sub,
                                          ggplot2::aes(x = subject_mid, y = !!yaxis, label = label),
                                          size = pattern.lim.size,
                                          family = font.family,
                                          inherit.aes = F)
  }
}

if (subject.lim.lines) {
  min.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_min(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
  max.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_max(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
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
  # maintain the original order of sequences
  out$seq.name <- factor(out$seq.name, levels = names(xstringset))

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
  seq_vector <- unique(toupper(seq_vector))
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

compare_seqs_df <- function(df,
                            ref,
                            pos_col = "position",
                            seq_col = "seq",
                            name_col = "seq.name",
                            seq_original = "seq_original") {

  seq.name.order <- levels(df[,name_col,drop=T])
  df <-
    df %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(!names(.) %in% c(ref, pos_col))])) %>%
    dplyr::mutate(value = paste0(value, "_", ifelse(value == !!rlang::sym(ref), ".", value))) %>%
    dplyr::mutate({{ref}} := paste0(!!rlang::sym(ref), "_", !!rlang::sym(ref))) %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(!names(.) %in% pos_col)]), names_to = name_col) %>%
    tidyr::separate(value, into = c(seq_original, seq_col), sep = "_") %>%
    dplyr::mutate({{name_col}} := factor(!!rlang::sym(name_col), levels = seq.name.order))
  return(df)
}

compare_seqs <- function(subject,
                         patterns,
                         match_dots = "pattern") {

  # adjusted from printPairwiseAlignment
  subject_name <- names(subject)
  subject <- strsplit(subject, "")[[1]]
  algnmt_chr <- purrr::map_chr(patterns, function(x) {
    pattern <- strsplit(x, "")[[1]]
    for (j in 1:length(pattern)) {
      if (pattern[j] == subject[j]) {
        if (match_dots == "subject") {
          subject[j] <- "."
        }
        if (match_dots == "pattern") {
          pattern[j] <- "."
        }
      }
    }
    return(paste(pattern, collapse = ""))
  })

  return(Biostrings::AAStringSet(c(algnmt_chr, stats::setNames(paste(subject, collapse = ""), subject_name))))
}
