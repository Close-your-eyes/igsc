#' Plot sequence alignments as ggplot object
#'
#' @param aln an alignment object returned from (i) igsc::pwalign_multi, (ii) DECIPHER::AlignSeqs or (iii) pwalign::pairwise_alignment;
#' in case (i) this is a data.frame, in case (ii) this is a XStringSet, in case (iii) this is a pairwise_alignmentsSingleSubject;
#' alternatively provide a custom data frame which has at least to contain pos_col, seq_col, name_col (see other arguments)
#' @param tile_fill a color scale for NTs or AAs (depending on type of aln); provide a named vector of colors where names are NTs or AAs;
#' or choose a name from igsc:::scheme_NT or igsc:::scheme_AA; or choose one from names(purrr::flatten(Peptides:::AAdata)); or leave NULL for the default scheme
#' @param tile_color a color to draw tile borders with; e.g. "black" or a hex code;
#' leave NA to have no borders (quicker plotting and advised for long alignments);
#' rather suited for short detailed alignments
#' @param tile_color_NA logical whether to draw borders on NA positions
#' @param tile_text logical whether to draw text in tiles (NT or AA identifier); or numeric indicating text size for geom_text()
#' advisable only for short alignments
#' @param base_theme ggplot theme to use as basis
#' @param pattern_lim_size numeric; plot annotation of pattern alignment limits; value indicates size;
#' set to 0 to omit plotting; only applicable if pa is provided
#' @param subject_lim_lines logical whether to plot vertical lines of subject alignment limits;
#' requires subject_name
#' @param pos_col name of position column in aln (applicable if aln is a data.frame)
#' @param seq_col name of sequence column in aln (applicable if aln is a data.frame)
#' @param name_col name of column which holds sequence names in aln (applicable if aln is a data.frame)
#' @param coord_fixed_ratio numeric; fixed aspect ratio of plot; leave NULL to not force a ratio
#' @param x_breaks numeric vector; manually provide breaks (ticks) on x-axis; set NULL to have breaks
#' picked automatically
#' @param aln_type typically 'NT' or 'AA' to influence the tile_fill automatically,
#' or any other string like "other" to have other colors;
#' only required if aln is a data.frame and only if
#' NT or AA cannot be guessed; leave NULL to have it guessed based on data
#' @param ref a reference sequence to compare all other sequences to;
#' e.g. a consensus sequence or a gene sequence to align reads to;
#' should be a name that appears in name_col; passed to compare_seq_df_long
#' @param subject_name name of the subject sequence in aln; subject will be plotted
#' first (at bottom of plot)
#' @param pattern_lim_pos where to plot the pattern_lims; only if pattern_lim_size > 0
#' @param pattern_names whether to plot pattern names within the plot
#' @param pairwise_alignment provide the pairwise_alignment that aln is based on;
#' this will enable the use of other function arguments
#' @param y_group_col name of column that contains information on which patterns
#' to plot in one row (as one group so to say); this argument competes with
#' group_on_yaxis (applicable if aln is a data.frame)
#' @param pattern_names_fun which function to use for plotting pattern_names;
#' e.g., ggrepel::geom_text_repel or ggplot2::geom_text
#' @param add_length_suffix whether to add the pattern length as suffix to
#' pattern_names
#' @param group_on_yaxis have an algorithm decide which patterns to plot in one
#' row without overlaps; this will use the plot area most efficiently;
#' this argument competes with y_group_col; use min_gap to define the minimal
#' gap size to other patterns
#' @param min_gap minimal gap size between patterns to trigger plotting in
#' separate rows; only applies if group_on_yaxis is TRUE
#' @param line_args arguments for line drawing; passed to geom_segment
#' @param theme_args theme arguments for ggplot
#' @param start_end_col name of column with start and end position of patterns;
#' not mandatory, can be derived from position column; but this may be wrong in
#' case when the first exon (or first part of a pattern in general) has a later
#' position as subsequent ones (think of circular reference sequences);
#' (applicable if aln is a data.frame)
#' @param pos_shift how many positions to shift the whole x-axis; e.g. in order
#' to make sure that the first part of a pattern comes most left;
#' can be a fixed position which becomes the start or a relative shift in form
#' of +100 or +2500 or so
#' @param pos_shift_adjust_axis adjust axis labels to maintain true positions
#' @param verbose print messages or not
#' @param focus focus on aligned patterns by cutting the subject range;
#' a positive integer limiting the range to n positions before the first pattern
#' position and n positions after the last pattern position
#' @param tile_line
#' @param base_theme_args
#' @param pattern_names_fun_args
#' @param y_order
#' @param y_breaks
#' @param order_numeric_seq_names
#'
#' @return ggplot2 object of alignment
#' @export
#'
#' @importFrom zeallot "%<-%"
#'
#' @examples
#' granzymes <- c("GZMA","GZMB","GZMH","GZMK","GZMM")
#' out <- get_sequences_from_biomart(granzymes)
#'
#' gzma <- dplyr::filter(out, hgnc_symbol == "GZMA")
#' gzmb <- dplyr::filter(out, hgnc_symbol == "GZMB")
#'
#' ### pairwise alignment
#' # just character vectors
#' aln1 <- pwalign::pairwise_alignment(subject = gzma$transcript_exon_intron,
#'                                    pattern = gzma$coding,
#'                                    type = "local-global")
#' aln_plot(aln1)
#'
#' # provide names via DNAStringSets
#' GZMA_RNA <- stats::setNames(gzma$transcript_exon_intron, "GZMA_premRNA")
#' GZMA_CDS <- stats::setNames(gzma$coding, "GZMA_CDS")
#' aln2 <- pwalign::pairwise_alignment(subject = Biostrings::DNAStringSet(GZMA_RNA),
#'                                    pattern = Biostrings::DNAStringSet(GZMA_CDS),
#'                                    type = "local-global")
#' aln_plot(aln2)
#'
#' ### multiple pairwise alignment
#' # GZMA is on plus strand
#' GZMA_exons <- gzma$gene_exon[[1]]$gene_exon
#' GZMA_gene <- gzma$gene_exon_intron
#' aln3 <- pwalign_multi(subject = GZMA_gene, patterns = GZMA_exons)
#' aln3[["plot"]]
#'
#' # GZMB is on minus strand, but also works
#' GZMB_exons <- gzmb$gene_exon[[1]]$gene_exon
#' GZMB_gene <- gzmb$gene_exon_intron
#' # no argument changed: focus on potential seq mismatches
#' aln4 <- pwalign_multi(subject = GZMB_gene, patterns = GZMB_exons)
#' aln4[["plot"]]
#' # different coloring and pattern in order of alignment position
#' aln5 <- pwalign_multi(subject = GZMB_gene, patterns = GZMB_exons,
#'                                                order_patterns = T,
#'                                                compare_seq_df_long_args = list(
#'                                                  change_ref = F,
#'                                                  change_nonref = T))
#' aln5[["plot"]]
#'
#' GZMB_CDS <- stats::setNames(gzmb$coding, "GZMB_CDS")
#' aln6 <- pwalign::pairwise_alignment(subject = Biostrings::DNAStringSet(GZMB_CDS),
#'                                    pattern = Biostrings::DNAStringSet(GZMA_CDS),
#'                                    type = "global")
#' # no modification of plotting colors
#' aln_plot(aln6)
#'
#' # alter matches/mismatches with compare_seq_df_long
#' padf <- pwalign_to_df(aln6)
#' padf2 <- compare_seq_df_long(padf,
#'                              ref = "GZMA_CDS",
#'                              change_nonref = T,
#'                              change_ref = T)
#' aln_plot(padf2)
#'
#' # change color of pwalign gaps
#' fillcol <- igsc:::scheme_NT[["Chemistry_NT"]
#' fillcol[which(names(fillcol) == "-")] <- NA # or "white"
#' aln_plot(padf2, tile_fill = fillcol)
aln_plot <- function(aln,
                     tile_fill = NULL,
                     tile_color = NA,
                     tile_color_NA = F,
                     tile_text = F,
                     tile_line = F,
                     line_args = list(linewidth = 0.1, color = "black"),
                     base_theme = colrr::theme_material,
                     base_theme_args = list(white = T),
                     theme_args = list(panel.grid = ggplot2::element_blank()),
                     pattern_lim_size = 0,
                     pattern_lim_pos = c("inner", "outer"),
                     pattern_names = 0, # make numeric, zero = disabled
                     pattern_names_fun = ggrepel::geom_text_repel,
                     pattern_names_fun_args = list(),
                     pairwise_alignment = NULL,
                     subject_lim_lines = F,
                     subject_name = NULL,
                     pos_col = "position",
                     seq_col = "seq",
                     name_col = "seq.name",
                     start_end_col = "start_end",
                     y_group_col = NULL,
                     coord_fixed_ratio = NULL,
                     x_breaks = NULL,
                     y_breaks = "..auto..",
                     aln_type = NULL,
                     add_length_suffix = F,
                     group_on_yaxis = F,
                     min_gap = 10,
                     ref = NULL,
                     pos_shift = NULL,
                     pos_shift_adjust_axis = T,
                     verbose = T,
                     focus = NULL,
                     y_order = c("as_is", "increasing", "decreasing"),
                     order_numeric_seq_names = F) {

  # document and clean up

  ## check duplicate names

  pattern_lim_pos <- rlang::arg_match(pattern_lim_pos)
  y_order <- rlang::arg_match(y_order)
  pattern_names_fun <- match.fun(pattern_names_fun)
  base_theme <- match.fun(base_theme)

  c(aln, aln_type, y_group_col) %<-% convert_aln_and_get_type(aln = aln,
                                                              aln_type = aln_type,
                                                              y_group_col = y_group_col,
                                                              verbose = verbose,
                                                              name_col = name_col,
                                                              seq_col = seq_col,
                                                              pos_col = pos_col,
                                                              order_numeric_seq_names = order_numeric_seq_names)


  if (!is.null(pairwise_alignment) && is.null(subject_name)) {
    if (!is.null(pairwise_alignment@subject@unaligned@ranges@NAMES)) {
      subject_name <- pairwise_alignment@subject@unaligned@ranges@NAMES
    }
  }

  c(subj_rngs, subject_name) %<-% infer_subject_name(aln = aln,
                                                     seq_col = seq_col,
                                                     name_col = name_col,
                                                     pos_col = pos_col,
                                                     subject_name = subject_name,
                                                     subject_name_infer = subject_name_infer)

  if (subject_lim_lines && is.null(subject_name)) {
    if (verbose) {
      message("subject_name is NULL or could not be inferred subject_lim_lines set to FALSE.")
    }
    subject_lim_lines <- F
  }

  if (subject_lim_lines && !is.null(subject_name) && !subject_name %in% aln[[name_col]]) {
    message("subject_name not found in ", name_col, " of aln. Can't draw subject_lim_lines.")
    subject_lim_lines <- F
  }

  aln <- change_aln_focus(aln = aln,
                          focus = focus,
                          subject_name = subject_name,
                          name_col = name_col,
                          pos_col = pos_col)

  c(aln, x_breaks) %<-% shift_aln(aln = aln,
                                  pos_shift = pos_shift,
                                  pos_col = pos_col,
                                  verbose = verbose,
                                  x_breaks = x_breaks)

  c(aln, yaxis, y_group_col) %<-% make_yaxis(aln = aln,
                                             group_on_yaxis = group_on_yaxis,
                                             y_group_col = y_group_col,
                                             y_order = y_order,
                                             subject_name = subject_name,
                                             subj_rngs = subj_rngs,
                                             min_gap = min_gap,
                                             name_col = name_col,
                                             pos_col = pos_col,
                                             seq_col = seq_col,
                                             verbose = verbose)

  c(aln, seq_original) %<-% compare_pattern_to_ref(aln = aln,
                                                   ref = ref,
                                                   pos_col = pos_col,
                                                   seq_col = seq_col,
                                                   name_col = name_col,
                                                   aln_type = aln_type)

  c(aln, tile_fill_internal, col_breaks, col_col) %<-% make_color(aln = aln,
                                                                  aln_type = aln_type,
                                                                  tile_fill = tile_fill,
                                                                  seq_col = seq_col,
                                                                  seq_original = seq_original,
                                                                  verbose = verbose)

  aln_summary <- make_aln_summary(aln = aln,
                                  seq_col = seq_col,
                                  pos_col = pos_col,
                                  name_col = name_col,
                                  start_end_col = start_end_col,
                                  y_group_col = y_group_col,
                                  yaxis = yaxis,
                                  add_length_suffix = add_length_suffix,
                                  pairwise_alignment = pairwise_alignment,
                                  verbose = verbose,
                                  aln_type = aln_type)


  if (y_breaks == "..auto.." && length(unique(aln[[yaxis]])) > 100) {
    if (verbose) {
      message("removed y axis breaks due to y_breaks=..auto.. and more than 100 breaks")
    }
    theme_args <- c(ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                   axis.ticks.y = ggplot2::element_blank()), theme_args)
  }

  plot <-
    ggplot2::ggplot(aln, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(yaxis))) + # name_col
    Gmisc::fastDoCall(base_theme, args = base_theme_args) +
    Gmisc::fastDoCall(ggplot2::theme, args = theme_args) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::coord_cartesian(expand = FALSE)

  if (add_length_suffix && is.null(y_group_col)) {
    # nt can only be added to y-axis text when sequences are not grouped
    plot <- plot + ggplot2::scale_y_discrete(labels = aln_summary[["label"]][levels(aln[[name_col]])])
  } else if (add_length_suffix && !is.null(y_group_col)) {
    message("length suffices not plotable when y-axis is grouped.")
  }

  if (is.numeric(aln[[col_col]])) {
    plot <- plot + ggplot2::scale_fill_viridis_c(na.value = "white")
  } else {
    plot <- plot + ggplot2::scale_fill_manual(values = tile_fill_internal,
                                              breaks = col_breaks,
                                              na.value = "white",
                                              na.translate = F)
  }

  plot <- add_tile_line(aln = aln,
                        tile_line = tile_line,
                        start_end_col = start_end_col,
                        verbose = verbose,
                        aln_summary = aln_summary,
                        line_args = line_args,
                        pos_col = pos_col,
                        yaxis = yaxis,
                        plot = plot)

  if (!is.na(tile_color)) {
    plot <- plot + ggplot2::geom_tile(data = if(tile_color_NA) {aln} else {aln[which(!is.na(aln[[seq_col]])),]},
                                      color = tile_color,
                                      ggplot2::aes(fill = !!rlang::sym(col_col)))
    # geom_raster seems not take color?!
  } else {
    plot <- plot + ggplot2::geom_tile(data = aln[which(!is.na(aln[[seq_col]])),],
                                      ggplot2::aes(fill = !!rlang::sym(col_col)))
    # geom_raster # geom_tile
  }

  if ((is.logical(tile_text) && tile_text) || (is.numeric(tile_text) && tile_text > 0)) {
    text_size <- ifelse(is.logical(tile_text), 4, tile_text)
    aln$text_colors <- brathering:::bw_txt(tile_fill_internal[as.character(aln[[col_col]])])
    plot <- plot +
      ggplot2::geom_text(data = aln, ggplot2::aes(label = !!rlang::sym(seq_col), color = I(text_colors)),
                         na.rm = T, size = text_size)
  }

  if (!is.null(coord_fixed_ratio)) {
    plot <- plot + ggplot2::coord_fixed(ratio = coord_fixed_ratio)
  }

  if (!is.null(pos_shift) && pos_shift_adjust_axis) {
    plot_x_breaks <- ggplot2::ggplot_build(plot)[["layout"]][["panel_params"]][[1]][["x"]][["breaks"]]
    plot_x_breaks <- plot_x_breaks[which(plot_x_breaks <= pos_shift)]
    plot <- suppressMessages(plot + ggplot2::scale_x_continuous(breaks = max(aln[[pos_col]]) - pos_shift + 1 + c(1, plot_x_breaks),
                                                                labels = c(1, plot_x_breaks)))
  }

  plot <- add_pattern_lims(plot = plot,
                           pattern_lim_size = pattern_lim_size,
                           pairwise_alignment = pairwise_alignment,
                           aln_summary = aln_summary,
                           name_col = name_col,
                           pos_col = pos_col,
                           pattern_lim_pos = pattern_lim_pos,
                           verbose = verbose)


  if (pattern_names > 0) {
    ## ggrepel considers the pattern.lim labels - great.
    plot <- plot + Gmisc::fastDoCall(pattern_names_fun,
                                     args = c(list(data = if(is.null(subject_name)) {aln_summary} else {dplyr::filter(aln_summary, !!rlang::sym(name_col) != subject_name)},
                                                   mapping = ggplot2::aes(x = mid_pos, y = !!rlang::sym(yaxis), label = label),
                                                   size = pattern_names,
                                                   inherit.aes = F),
                                              pattern_names_fun_args))
  }

  if (subject_lim_lines) {
    #minmaxpos <- c(aln_summary |> dplyr::filter(!!rlang::sym(name_col) != subject_name) |> dplyr::pull(min_pos), aln_summary |> dplyr::filter(!!rlang::sym(name_col) != subject_name) |> dplyr::pull(max_pos))
    min.pos <- aln |> dplyr::filter(!!rlang::sym(seq_col) != "-") |> dplyr::filter(!is.na(!!rlang::sym(seq_col))) |> dplyr::filter(!!rlang::sym(name_col) != subject_name) |> dplyr::slice_min(order_by = !!rlang::sym(pos_col), n = 1) |> dplyr::pull(!!rlang::sym(pos_col))
    max.pos <- aln |> dplyr::filter(!!rlang::sym(seq_col) != "-") |> dplyr::filter(!is.na(!!rlang::sym(seq_col))) |> dplyr::filter(!!rlang::sym(name_col) != subject_name) |> dplyr::slice_max(order_by = !!rlang::sym(pos_col), n = 1) |> dplyr::pull(!!rlang::sym(pos_col))
    plot <- plot + ggplot2::geom_vline(xintercept = c(min.pos, max.pos), linetype = "dashed")
  }


  return(plot)
}


convert_aln_and_get_type <- function(aln,
                                     aln_type,
                                     y_group_col,
                                     verbose,
                                     name_col,
                                     seq_col,
                                     pos_col,
                                     order_numeric_seq_names) {

  if (methods::is(aln, "DNAStringSet") || methods::is(aln, "RNAStringSet") || methods::is(aln, "AAStringSet")) {
    # from DECIPHER::AlignSeqs
    if (methods::is(aln, "DNAStringSet") || methods::is(aln, "RNAStringSet")) {
      aln_type <- "NT"
    }
    if (methods::is(aln, "AAStringSet")) {
      aln_type <- "AA"
    }
    aln <- xstringset_to_df(xstringset = aln,
                            name_col = name_col,
                            seq_col = seq_col,
                            pos_col = pos_col)
    if (!is.null(y_group_col)) {
      ## y_group_col makes only sense when a aln is a dataframe, e.g. from multi pairwise_alignment
      if (verbose) {
        message("y_group_col is set to NULL.")
      }
      y_group_col <- NULL
    }
  } else if (methods::is(aln, "pairwise_alignmentsSingleSubject") || methods::is(aln, "list")) {
    # from pwalign::pairwise_alignment
    if (methods::is(aln, "pairwise_alignmentsSingleSubject")) {
      aln_type <- ifelse(guess_type2(aln) == "AA", "AA", "NT")
    } else if (methods::is(aln, "list")) {
      aln_type <- ifelse(guess_type2(aln[[1]]) == "AA", "AA", "NT")
    }

    aln <- pwalign_to_df(
      pa = aln,
      verbose = verbose,
      subject_width = "whole"
    )
    if (!is.null(y_group_col)) {
      if (verbose) {
        message("y_group_col is set to NULL.")
      }
      y_group_col <- NULL
    }
  }

  if (!name_col %in% names(aln)) {
    stop("name_col not found in aln.")
  }
  if (!seq_col %in% names(aln)) {
    stop("seq_col not found in aln.")
  }
  if (!pos_col %in% names(aln)) {
    stop("pos_col not found in aln.")
  }
  if (!is.null(y_group_col) && !y_group_col %in% names(aln)) {
    stop("y_group_col not found in aln.")
  }

  if (is.null(aln_type)) {
    inds <- intersect(which(!is.na(aln[[seq_col]])), which(!aln[[seq_col]] %in% c(".", "x", "match", "mismatch", "gap", "insertion", "ambiguous")))
    aln_type <- guess_type(seq_vector = aln[[seq_col]][inds])
  }

  # make sure it is not a factor
  if (!is.numeric(aln[[pos_col]])) {
    aln[[pos_col]] <- as.numeric(as.character(aln[[pos_col]]))
  }
  aln[[seq_col]] <- as.character(aln[[seq_col]])

  # if seq names are numeric; order them increasingly
  if (order_numeric_seq_names && !anyNA(suppressWarnings(as.numeric(as.character(aln[[name_col]]))))) {
    aln[[name_col]] <- factor(aln[[name_col]], levels = as.character(unique(as.numeric(as.character(aln[[name_col]])))))
  } else if (!is.factor(aln[[name_col]])) {
    aln[[name_col]] <- as.factor(aln[[name_col]])
  }

  return(list(aln = aln, aln_type = aln_type, y_group_col = y_group_col))
}


change_aln_focus <- function(aln,
                             focus,
                             subject_name,
                             name_col,
                             pos_col) {
  if (!is.null(focus) && !is.null(subject_name)) {
    if (!is.numeric(focus)) {
      stop("focus should be a positive integer.")
    }
    pos_col_rng <- range(aln[which(aln[[name_col]] != subject_name), pos_col])
    aln <- dplyr::filter(aln, dplyr::between(!!rlang::sym(pos_col), pos_col_rng[1]-focus, pos_col_rng[2]+focus))
  } else if (!is.null(focus) && is.null(subject_name)) {
    message("focus requires subject_name which is NULL though.")
  } else if (!is.null(focus) && !is.null(subject_name) && !subject_name %in% aln[[name_col]]) {
    message("subject name not found in ", name_col, ". Cannot use focus.")
  }

  return(aln)
}

shift_aln <- function(aln,
                      pos_shift,
                      pos_col,
                      verbose,
                      x_breaks) {

  if (!is.null(pos_shift)) {
    if (!is.numeric(pos_shift) && !grepl("^\\+", pos_shift)) {
      stop("pos_shift has to be numeric giving and absolute position as start or start with a '+' giving a relative shift.")
    }
    if (grepl("^\\+", pos_shift)) {
      pos_shift <- max(aln[[pos_col]]) - as.numeric(gsub("\\+", "", pos_shift)) + 2
      if (pos_shift < 0) {
        stop("pos_shift cannot be larger than the largest alignment position.")
      }
    }
    conv <- stats::setNames(shifted_pos(
      x = unique(aln[[pos_col]]),
      start_pos = pos_shift,
      verbose = verbose
    ), unique(aln[[pos_col]]))
    aln[[pos_col]] <- conv[aln[[pos_col]]]
  }

  # https://stackoverflow.com/questions/45493163/ggplot-remove-na-factor-level-in-legend
  # --> na.translate = F
  if (is.null(x_breaks)) {
    ## fix with expand = F in coord_cartesian
    x_breaks <- floor(base::pretty(unique(aln[[pos_col]]), n = 5))
    if (!any(x_breaks < 0)) {
      x_breaks <- x_breaks[which(x_breaks > 0)]
    }
    #x_breaks <- x_breaks[which(x_breaks < max(aln[[pos_col]]))]
  }

  return(list(aln = aln, x_breaks = x_breaks))
}


make_yaxis <- function(aln,
                       group_on_yaxis,
                       y_group_col,
                       y_order,
                       subject_name,
                       subj_rngs,
                       min_gap,
                       name_col,
                       pos_col,
                       seq_col,
                       verbose) {

  if (group_on_yaxis && is.null(y_group_col) && length(unique(aln$seq.name))>2) {
    # if pairwise alignment is provided, use this one to derived ranges - but this may not be compatible with shifting?!
    if (is.null(subject_name)) {
      message("Cannot group on y-axis without subject name.")
    } else {
      subj_rngs <- subj_rngs[names(subj_rngs) != subject_name]
      if (length(subj_rngs) == 0) {
        rows <- character(0)
      } else {
        rng_df <- data.frame(
          name = names(subj_rngs),
          start = vapply(subj_rngs, min, numeric(1)),
          end   = vapply(subj_rngs, max, numeric(1)),
          stringsAsFactors = FALSE
        )
        rng_df <- rng_df[order(rng_df$start, rng_df$end), ]
        assigned_rows <- integer(nrow(rng_df))
        assigned_rows[1] <- 1

        for (i in 2:nrow(rng_df)) {
          current_start <- rng_df$start[i]
          current_end   <- rng_df$end[i]
          candidate_row <- 1
          repeat {
            existing <- which(assigned_rows == candidate_row)
            if (length(existing) == 0) {
              break
            }
            too_close <- any(vapply(existing, function(j) {
              dist <- range_distance(
                rng_df$start[j],
                rng_df$end[j],
                current_start,
                current_end
              )
              dist < min_gap
            }, logical(1)))
            if (!too_close) {
              break
            }
            candidate_row <- candidate_row + 1
          }
          assigned_rows[i] <- candidate_row
        }
        rows <- paste0("Group_", assigned_rows)
        names(rows) <- rng_df$name
      }

      rows <- c(stats::setNames(subject_name, subject_name), rows)
      y_group_col <- "group"
      aln[[y_group_col]] <- rows[as.character(aln[[name_col]])]
      aln[[y_group_col]] <- factor(
        aln[[y_group_col]],
        levels = c(subject_name, unique(rows[-1]))
      )
    }

  } else if (group_on_yaxis && !is.null(y_group_col) && length(unique(aln$seq.name))>2) {
    if (verbose) {
      message("y_group_col is not NULL. Using this one. Ignoring group_on_yaxis.")
    }
  }

  if (y_order != "as_is" && length(unique(aln$seq.name))>2) {
    # if (is.null(y_group_col) && is.factor(aln[[name_col]])) {
    #   message("name_col is a factor. Will stick to this order.")
    #   order <- F
    # } else if (!is.null(y_group_col) && is.factor(aln[[y_group_col]])) {
    #   message("y_group_col is a factor. Will stick to this order.")
    #   order <- F
    # }
    if (is.null(y_group_col)) {
      # second level of ordering: length: shortest first
      yorder <- data.frame(name = names(subj_rngs),
                           minpos = sapply(subj_rngs, min, simplify = T),
                           len = sapply(subj_rngs, length, simplify = T)) |>
        dplyr::arrange(minpos, len) |>
        dplyr::pull(name)
      if (y_order == "decreasing") {
        yorder <- rev(yorder)
      }
      # if (!is.null(subject_name)) {
      #   yorder <- unique(c(subject_name, yorder))
      # }
      yorder <- unique(c(subject_name, yorder))
      aln[[name_col]] <- factor(aln[[name_col]], levels = yorder)
    } else {
      yorder <- dplyr::summarise(aln, min = min(!!rlang::sym(pos_col)), .by = !!rlang::sym(y_group_col))
      yorder <- yorder[which(yorder[[y_group_col]] != subject_name),]
      yorder <- names(sort(stats::setNames(yorder$min, yorder[,1,drop=T]), decreasing = y_order == "decreasing"))
      yorder <- c(subject_name, yorder)
      # yorder <- yorder |>
      #   dplyr::arrange(!!rlang::sym(pos_col))) |>
      #   dplyr::distinct(!!rlang::sym(y_group_col))
      # if (!is.null(subject_name)) {
      #   subject_name2 <- aln |>
      #     dplyr::filter(!!rlang::sym(name_col) == subject_name) |>
      #     dplyr::distinct(!!rlang::sym(y_group_col)) |>
      #     dplyr::pull(!!rlang::sym(y_group_col))
      #   # yorder <- dplyr::filter(yorder, !!rlang::sym(y_group_col) != subject_name2) # unique does the job
      #   yorder <- unique(c(subject_name2, names(sort(stats::setNames(yorder$min, yorder[[1,drop=T]), decreasing = y_order == "decreasing"))))
      # } else {
      #   yorder <- names(sort(stats::setNames(yorder$min, yorder[[1,drop=T]), decreasing = y_order == "decreasing"))
      # }
      aln[[y_group_col]] <- factor(aln[[y_group_col]], levels = yorder)
    }
  }

  if (!is.null(y_group_col)) {
    test_dups <- aln |>
      tidyr::drop_na(!!rlang::sym(seq_col)) |> # y_group_col
      dplyr::summarise(n = dplyr::n(), .by = c(!!rlang::sym(y_group_col), !!rlang::sym(pos_col)))

    if (any(test_dups$n > 1)) {
      message("Your y_group_col causes multiple sequnces (subject/patterns) to occupy the same position in the plot. You may want to check on that.")
      message("dplyr::summarise(tidyr::drop_na(aln, !!rlang::sym(y_group_col)), n = dplyr::n(), .by = c(!!rlang::sym(y_group_col), !!rlang::sym(pos_col)))")
    }
  }

  #yaxis <- ifelse(is.null(y_group_col), rlang::sym(name_col), rlang::sym(y_group_col))
  yaxis <- ifelse(is.null(y_group_col), name_col, y_group_col)

  return(list(aln = aln, yaxis = yaxis, y_group_col = y_group_col))
}

range_distance <- function(a_start, a_end, b_start, b_end) {

  # overlapping intervals
  if (a_end >= b_start && b_end >= a_start) {
    return(-min(a_end, b_end) + max(a_start, b_start) - 1)
  }

  # a before b
  if (a_end < b_start) {
    return(b_start - a_end)
  }

  # b before a
  return(a_start - b_end)
}

compare_pattern_to_ref <- function(aln,
                                   ref,
                                   pos_col,
                                   seq_col,
                                   name_col,
                                   aln_type) {

  if (!is.null(ref)) {
    if (!ref %in% unique(aln[[name_col]])) {
      message("ref not found in names of aln.")
      ref <- NULL
    }
    if (aln_type == "AA") {
      seq_original <- "seq_original"
    } else {
      seq_original <- NULL
    }

    aln <- compare_seq_df_long(df = aln,
                               ref = ref,
                               pos_col = pos_col,
                               seq_col = seq_col,
                               name_col = name_col,
                               seq_original = seq_original,
                               match_symbol = ".",
                               mismatch_symbol = "x",
                               change_nonref = T,
                               nonref_mismatch_as = "mismatch_symbol",
                               change_ref = F)

    # seq_original will cause that tiles are filled according to chemical property even if they are replaced by a dots when
    # sequence is equal to the reference sequence (ref)
  } else {
    seq_original <- seq_col
  }

  if (!all(c(pos_col, seq_col, name_col) %in% names(aln))) {
    stop("after comparison to ref: aln at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ". Alternatively change function arguments.")
  }

  return(list(aln = aln, seq_original = seq_original))
}


make_color <- function(aln,
                       aln_type,
                       tile_fill,
                       seq_col,
                       seq_original,
                       verbose) {

  NT <- igsc:::scheme_NT
  AA <- igsc:::scheme_AA

  # use preset colors
  if (is.null(tile_fill)) {
    if (aln_type == "NT") {
      tile_fill <- colnames(NT)[1]
    }
    if (aln_type == "AA") {
      tile_fill <- colnames(AA)[1]
      # set tile_fill to Chemistry_AA by default which will cause falling into the special case below
      # to avoid that, change this here to tile_color (like in case of aln_type == "NT") and connect
    }
  }

  #default; may change when aln_type == "AA" && tile_fill == "Chemistry_AA"
  col_col <- seq_col
  if (length(tile_fill) == 1 && tile_fill %in% c(colnames(AA),
                                                 colnames(NT),
                                                 names(purrr::flatten(Peptides:::AAdata)))) {

    if (aln_type == "NT") {
      tile_fill_internal <- NT[,match.arg(tile_fill, choices = colnames(NT))]
    } else if (aln_type == "AA") {
      if (tile_fill == "Chemistry_AA") {
        # special case; change legend to chemical property
        col_col <- tile_fill
        chem_col <- utils::stack(igsc:::aa_info[["aa_main_prop"]])
        names(chem_col) <- c(col_col, seq_original)
        chem_col[[seq_original]] <- as.character(chem_col[[seq_original]])
        aln <- dplyr::left_join(aln, chem_col, by = seq_original)
        aln[[col_col]] <- ifelse(is.na(aln[[col_col]]), aln[[seq_original]], aln[[col_col]])
        tile_fill_internal <- AA[match.arg(tile_fill, choices = colnames(AA)),]
        tile_fill_internal <- tile_fill_internal[unique(aln[[seq_col]][which(!is.na(aln[[seq_col]]))])]
        for (i in 1:length(tile_fill_internal)) {
          if (names(tile_fill_internal)[i] %in% names(igsc:::aa_info[["aa_main_prop"]])) {
            names(tile_fill_internal)[i] <- igsc:::aa_info[["aa_main_prop"]][names(tile_fill_internal)[i]]
          }
        }
        tile_fill_internal <- tile_fill_internal[unique(names(tile_fill_internal))]
      } else if (tile_fill %in% names(purrr::flatten(Peptides:::AAdata))) {
        col_col <- tile_fill
        aln[[col_col]] <- purrr::flatten(Peptides:::AAdata)[[tile_fill]][aln[[seq_original]]]
        tile_fill_internal <- colrr::col_pal("viridis", return = "c", n = 100)[cut(unique(aln[[col_col]])[which(!is.na(unique(aln[[col_col]])))], 100)]
        names(tile_fill_internal) <- unique(aln[[col_col]])[which(!is.na(unique(aln[[col_col]])))]
      } else {
        tile_fill_internal <- AA[match.arg(tile_fill, choices = colnames(AA)),]
      }
    } else {
      if (verbose) {
        message("Type of alignment data (NT or AA) could not be determined. Choosing default ggplot colors.")
      }
      tile_fill_internal <- scales::hue_pal()(length(unique(as.character(aln[[seq_col]][which(!is.na(aln[[seq_col]]))]))))
    }
  } else {
    # provide own colors
    if (is.null(tile_fill)) {
      # happens when algmnt_type is not NT or AA
      # redundant to above when NT or AA cannot be guessed
      tile_fill <- scales::hue_pal()(length(unique(as.character(aln[[seq_col]][which(!is.na(aln[[seq_col]]))]))))
    } else {
      if (is.null(names(tile_fill))) {
        if (length(tile_fill) < length(unique(aln[[seq_col]][which(!is.na(aln[[seq_col]]))]))) {
          warning("Less colors provided than entities in the alignment.")
        }
      } else {
        if (any(!unique(aln[[seq_col]][which(!is.na(aln[[seq_col]]))]) %in% names(tile_fill))) {
          warning("Not all values in aln[[seq_col]] found in names(tile_fill).")
        }
      }
    }
    tile_fill_internal <- tile_fill
  }

  if (!is.null(names(tile_fill_internal)) && col_col == seq_col) {
    tile_fill_internal <- tile_fill_internal[unique(aln[[seq_col]][which(!is.na(aln[[seq_col]]))])]
  }

  if (aln_type == "AA" && length(tile_fill) == 1 && tile_fill == "Chemistry_AA") {
    col_breaks = unique(igsc:::aa_info[["aa_main_prop"]])[which(unique(igsc:::aa_info[["aa_main_prop"]]) != "stop")]
  } else {
    col_breaks = ggplot2::waiver()
  }

  return(list(aln = aln, tile_fill_internal = tile_fill_internal, col_breaks = col_breaks, col_col = col_col))
}


make_aln_summary <- function(aln,
                             seq_col,
                             pos_col,
                             name_col,
                             start_end_col,
                             y_group_col,
                             yaxis,
                             add_length_suffix,
                             pairwise_alignment,
                             verbose,
                             aln_type) {

  ## aln summary prep
  ## exclude subject ??
  aln_summary <- aln |>
    dplyr::filter(!is.na(!!rlang::sym(seq_col))) |>
    dplyr::summarise(n_not_NA = sum(!is.na(!!rlang::sym(seq_col))),
                     min_pos = min(!!rlang::sym(pos_col)),
                     max_pos = max(!!rlang::sym(pos_col)),
                     .by = !!rlang::sym(name_col)) |>
    dplyr::mutate(mid_pos = max_pos - (max_pos-min_pos)/2)


  if (start_end_col %in% names(aln)) {
    #seq_original may come from compare_seqs_df, when !is.null(ref)
    aln_summary <- aln_summary |>
      dplyr::left_join(aln |>
                         dplyr::filter(!is.na(!!rlang::sym(start_end_col))) |>
                         dplyr::select(!!!rlang::syms(names(aln)[which(!names(aln) %in% c(seq_col, "seq_original"))])) |>
                         tidyr::pivot_wider(names_from = !!rlang::sym(start_end_col), values_from = !!rlang::sym(pos_col)),
                       by = name_col)
  }

  if (!is.null(y_group_col)) {
    y_group_conv <- dplyr::distinct(tidyr::drop_na(aln, !!rlang::sym(seq_col)), !!rlang::sym(y_group_col), !!rlang::sym(name_col))
    y_group_conv <- stats::setNames(as.character(y_group_conv[[y_group_col]]), y_group_conv[[name_col]])
    aln_summary[[yaxis]] <- y_group_conv[as.character(aln_summary[[name_col]])]
    if (is.factor(aln[[yaxis]])) {
      aln_summary[[yaxis]] <- factor(aln_summary[[yaxis]], levels = levels(aln[[yaxis]]))
    } else {
      aln_summary[[yaxis]] <- factor(aln_summary[[yaxis]], levels = unique(aln_summary[[yaxis]]))
    }
  }

  if (add_length_suffix) {
    # keep this, with pairwise_alignment
    if (is.null(pairwise_alignment)) {
      if (verbose) {
        message("pairwise_alignment is not provided and orignal sequences are unknown. Length suffix is inferred from number of non-NA elements in seq_col.")
      }
      seq_lengths <- stats::setNames(aln_summary[["n_not_NA"]],
                                     aln_summary[[name_col]])
    } else {
      if (verbose) {
        message("Length suffix is inferred from pairwise_alignment.")
      }
      seq_lengths <- stats::setNames(c(pairwise_alignment@subject@unaligned@ranges@width,
                                       pairwise_alignment@pattern@unaligned@ranges@width),
                                     c(pairwise_alignment@subject@unaligned@ranges@NAMES,
                                       pairwise_alignment@pattern@unaligned@ranges@NAMES))
    }
    if (aln_type %in% c("NT", "AA")) {
      seq_lengths <- stats::setNames(paste0(names(seq_lengths), "\n", seq_lengths, " ", tolower(aln_type)),
                                     nm = names(seq_lengths))
    } else {
      seq_lengths <- stats::setNames(paste0(names(seq_lengths), "\n", seq_lengths),
                                     nm = names(seq_lengths))
    }

    aln_summary$label <- seq_lengths[aln_summary[[name_col]]]
  } else {
    aln_summary$label <- stats::setNames(as.character(aln_summary[[name_col]]),
                                         as.character(aln_summary[[name_col]]))
  }

  return(aln_summary)
}


add_pattern_lims <- function(plot,
                             pattern_lim_size,
                             pairwise_alignment,
                             aln_summary,
                             name_col,
                             pos_col,
                             pattern_lim_pos,
                             verbose) {

  if (pattern_lim_size > 0 && !is.null(pairwise_alignment)) {
    pattern.ranges <- data.frame(pairwise_alignment@pattern@range,
                                 seq.name = ifelse(rep(is.null(pairwise_alignment@pattern@unaligned@ranges@NAMES), length(pairwise_alignment)),
                                                   paste0("pattern_", seq(1,length(pairwise_alignment))),
                                                   pairwise_alignment@pattern@unaligned@ranges@NAMES))
    names(pattern.ranges)[1:2] <- paste0("pattern_", names(pattern.ranges)[1:2])
    # stats::setNames(as.data.frame(pairwise_alignment@subject@range)[1:2], nm = paste0("subject_", names(as.data.frame(pairwise_alignment@subject@range)[1:2])))
    # subject position not from pairwise_alignment but from aln_summary as the latter can account for pos_shift
    ranges <- dplyr::left_join(pattern.ranges, aln_summary, by = dplyr::join_by(!!rlang::sym(name_col)))
    ranges <- tidyr::pivot_longer(ranges, cols = c(min_pos, max_pos), names_to = "pos", values_to = "inner_position")
    ranges$outer_position <- ifelse(ranges$pos == "min_pos", -2, max(aln[[pos_col]]) + 2)
    ranges$label <- ifelse(ranges$pos == "min_pos", ranges$pattern_start, ranges$pattern_end)
    # maybe adjust style of labels
    label_plot_pos <- ifelse(pattern_lim_pos == "inner", rlang::sym("inner_position"), rlang::sym("outer_position"))
    # nudge labels alternating up and down?
    plot <- plot + ggplot2::geom_text(data = ranges,
                                      ggplot2::aes(x = !!label_plot_pos, y = !!rlang::sym(yaxis), label = label),
                                      size = pattern_lim_size,
                                      inherit.aes = F)
  } else if (pattern_lim_size > 0) {
    if (verbose) {
      message("pattern limits can only be plotted when pairwise_alignment is provided.")
    }
  }

  return(plot)

}

add_tile_line <- function(aln,
                          tile_line,
                          start_end_col,
                          verbose,
                          aln_summary,
                          line_args,
                          pos_col,
                          yaxis,
                          plot) {

  if (tile_line) {
    if (!start_end_col %in% names(aln)) {
      if (verbose) {
        message("start_end_col not found in aln data frame. Using min and max position to draw line. But this could be wrong if, e.g. the first exon is at later position as the first one.")
      }
      xmin <- "min_pos"
      xmax <- "max_pos"
    } else {
      xmin <- "start"
      xmax <- "end"
    }

    if (start_end_col %in% names(aln) && any(aln_summary$start > aln_summary$end)) {
      # split the plotting of line_segment
      # start to very last position
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(aln_summary, start > end),
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = max(aln[[pos_col]]), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # very first position to very end
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(aln_summary, start > end),
                                                                            ggplot2::aes(x = 1, xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # plot those where start < end
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(aln_summary, start < end),
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # reordering of axis needed as step-wise plotting of lines corrupts the original order
      plot <- suppressMessages(plot + ggplot2::scale_y_discrete(limits = levels(aln[[yaxis]])))
    } else {
      # all starts are smaller than ends
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = aln_summary,
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
    }
  }

  return(plot)
}


#' Convert XStringSet to data frame
#'
#' @param xstringset
#' @param name_col column with seq names
#' @param seq_col column with seq data (nt or aa)
#' @param pos_col position column
#' @param subject_name name of sequence that is subject or reference;
#' is passed as attribute to return df
#'
#' @returns long data frame
#' @export
#'
#' @examples
xstringset_to_df <- function(xstringset,
                             name_col = "seq.name",
                             seq_col = "seq",
                             pos_col = "position",
                             subject_name = NULL) {
  terminal_gap_sym <- "&"
  out <- purrr::map(as.list(xstringset), as.character)
  out <- purrr::map(out, replace_terminal_dashes, replacement = terminal_gap_sym)
  out <- purrr::flatten(purrr::map(out, strsplit, split = ""))
  out <- purrr::map(out, ~gsub(terminal_gap_sym, NA, .x))

  out <- purrr::map_dfr(out, function(x) utils::stack(stats::setNames(x, seq(1, length(x)))), .id = name_col)
  names(out)[c(2,3)] <- c(seq_col, pos_col)
  out[[pos_col]] <- as.numeric(as.character(out[[pos_col]]))
  # maintain the original order of sequences
  out[[name_col]] <- factor(out[[name_col]], levels = names(xstringset))

  if (!is.null(subject_name)) {
    attr(out, "subject_name") <- subject_name
  }

  # do this in parent fun
  #if (!anyNA(suppressWarnings(as.numeric(out$seq.name)))) {out$seq.name <- factor(out$seq.name, levels = as.character(unique(as.numeric(as.character(out$seq.name)))))}
  return(out)
}

replace_terminal_dashes <- function(x, replacement = "*") {
  stopifnot(is.character(x), is.character(replacement), nchar(replacement) == 1L)

  # replace start-of-string runs
  m1 <- regexpr("^-+", x, perl = TRUE)
  repl1 <- regmatches(x, m1)
  len1  <- attr(m1, "match.length")
  hit1  <- which(m1 != -1L)
  for (i in hit1) repl1[[i]] <- strrep(replacement, len1[i])
  regmatches(x, m1) <- repl1

  # replace end-of-string runs
  m2 <- regexpr("-+$", x, perl = TRUE)
  repl2 <- regmatches(x, m2)
  len2  <- attr(m2, "match.length")
  hit2  <- which(m2 != -1L)
  for (i in hit2) repl2[[i]] <- strrep(replacement, len2[i])
  regmatches(x, m2) <- repl2

  return(x)
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

compare_seqs <- function(subject,
                         patterns,
                         match_dots = "pattern") {

  # adjusted from pwalign_print
  subject_name <- names(subject)
  subject <- strsplit(subject, "")[[1]]
  aln_chr <- purrr::map_chr(patterns, function(x) {
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

  return(Biostrings::AAStringSet(c(aln_chr, stats::setNames(paste(subject, collapse = ""), subject_name))))
}

## work with start pos - then derive shift
shifted_pos <- function(x,
                        n = NULL,
                        start_pos = NULL,
                        verbose) {


  if (!is.null(n) && is.na(n)) {
    stop("n is NA, not numeric.")
  }
  if (!is.null(start_pos) && is.na(start_pos)) {
    stop("start_pos is NA, not numeric.")
  }

  if (!is.null(start_pos)) {
    if (!is.null(n)) {
      if (verbose) {
        message("n will be ignored as start_pos is provided.")
      }
    }
    n <- which(start_pos == x) - 1
    if (length(n) == 0) {
      stop("start_pos not found in x.")
    }
  }
  x2 <- c(x[(length(x)-n+1):length(x)], dplyr::lag(x,n)[(n+1):length(x)])
  return(x2)
}

infer_subject_name <- function(aln,
                               seq_col,
                               name_col,
                               pos_col,
                               subject_name,
                               subject_name_infer = T,
                               verbose = T) {

  subj_rngs <- aln |>
    tidyr::drop_na(!!rlang::sym(seq_col)) |>
    dplyr::group_by(!!rlang::sym(name_col))
  subj_rngs <- stats::setNames(igsc:::seq2(subj_rngs |>
                                             dplyr::slice_min(!!rlang::sym(pos_col)) |>
                                             dplyr::pull(!!rlang::sym(pos_col)),
                                           subj_rngs |>
                                             dplyr::slice_max(!!rlang::sym(pos_col)) |>
                                             dplyr::pull(!!rlang::sym(pos_col))),
                               nm = as.character(subj_rngs |>
                                                   dplyr::slice_min(!!rlang::sym(pos_col)) |>
                                                   dplyr::pull(!!rlang::sym(name_col))))
  #assign("subj_rngs", subj_rngs, envir = parent.frame())


  if (is.null(subject_name)) {
    # just try, possible from pa_to_df
    # remains NULL if does not exist
    subject_name <- attr(aln, "subject_name")
  }

  if (is.null(subject_name) && subject_name_infer) {
    max_len <- max(lengths(subj_rngs))
    if (length(which(lengths(subj_rngs) == max_len)) > 1) {
      stop("Subject could not be identified.")
      ## TODO: set variable for ordering pattern on y-axis to FALSE here
    } else {
      subject_name <- names(which(lengths(subj_rngs) == max_len))
    }
    # assigns in parent environment (https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r?rq=1)
    #assign("subject_name", subject_name, envir = parent.frame())
    if (verbose) {
      message("subject sequence inferred: ", subject_name)
    }
  }

  return(list(subj_rngs = subj_rngs, subject_name = subject_name))
}
