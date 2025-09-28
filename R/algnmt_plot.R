#' Plot sequence alignments as ggplot object
#'
#' @param algnmt an alignment object returned from (i) igsc::pwalign_multi, (ii) DECIPHER::AlignSeqs or (iii) pwalign::pairwiseAlignment;
#' in case (i) this is a data.frame, in case (ii) this is a XStringSet, in case (iii) this is a PairwiseAlignmentsSingleSubject;
#' alternatively provide a custom data frame which has at least to contain pos_col, seq_col, name_col (see other arguments)
#' @param tile_fill a color scale for NTs or AAs (depending on type of algnmt); provide a named vector of colors where names are NTs or AAs;
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
#' @param pos_col name of position column in algnmt (applicable if algnmt is a data.frame)
#' @param seq_col name of sequence column in algnmt (applicable if algnmt is a data.frame)
#' @param name_col name of column which holds sequence names in algnmt (applicable if algnmt is a data.frame)
#' @param coord_fixed_ratio numeric; fixed aspect ratio of plot; leave NULL to not force a ratio
#' @param x_breaks numeric vector; manually provide breaks (ticks) on x-axis; set NULL to have breaks
#' picked automatically
#' @param algnmt_type typically 'NT' or 'AA' to influence the tile_fill automatically,
#' or any other string like "other" to have other colors;
#' only required if algnmt is a data.frame and only if
#' NT or AA cannot be guessed; leave NULL to have it guessed based on data
#' @param ref a reference sequence to compare all other sequences to;
#' e.g. a consensus sequence or a gene sequence to align reads to;
#' should be a name that appears in name_col; passed to compare_seq_df_long
#' @param subject_name name of the subject sequence in algnmt; subject will be plotted
#' first (at bottom of plot)
#' @param pattern_lim_pos where to plot the pattern_lims; only if pattern_lim_size > 0
#' @param pattern_names whether to plot pattern names within the plot
#' @param pairwiseAlignment provide the pairwiseAlignment that algnmt is based on;
#' this will enable the use of other function arguments
#' @param y_group_col name of column that contains information on which patterns
#' to plot in one row (as one group so to say); this argument competes with
#' group_on_yaxis (applicable if algnmt is a data.frame)
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
#' (applicable if algnmt is a data.frame)
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
#' @param subject_name_infer
#' @param y_order
#'
#' @return ggplot2 object of alignment
#' @export
#'
#' @importFrom magrittr "%>%"
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
#' aln1 <- pwalign::pairwiseAlignment(subject = gzma$transcript_exon_intron,
#'                                    pattern = gzma$coding,
#'                                    type = "local-global")
#' algnmt_plot(aln1)
#'
#' # provide names via DNAStringSets
#' GZMA_RNA <- stats::setNames(gzma$transcript_exon_intron, "GZMA_premRNA")
#' GZMA_CDS <- stats::setNames(gzma$coding, "GZMA_CDS")
#' aln2 <- pwalign::pairwiseAlignment(subject = Biostrings::DNAStringSet(GZMA_RNA),
#'                                    pattern = Biostrings::DNAStringSet(GZMA_CDS),
#'                                    type = "local-global")
#' algnmt_plot(aln2)
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
#' aln6 <- pwalign::pairwiseAlignment(subject = Biostrings::DNAStringSet(GZMB_CDS),
#'                                    pattern = Biostrings::DNAStringSet(GZMA_CDS),
#'                                    type = "global")
#' # no modification of plotting colors
#' algnmt_plot(aln6)
#'
#' # alter matches/mismatches with compare_seq_df_long
#' padf <- pwalign_to_df(aln6)
#' padf2 <- compare_seq_df_long(padf,
#'                              ref = "GZMA_CDS",
#'                              change_nonref = T,
#'                              change_ref = T)
#' algnmt_plot(padf2)
#'
#' # change color of pwalign gaps
#' fillcol <- igsc:::scheme_NT[,"Chemistry_NT"]
#' fillcol[which(names(fillcol) == "-")] <- NA # or "white"
#' algnmt_plot(padf2, tile_fill = fillcol)
algnmt_plot <- function(algnmt,
                        tile_fill = NULL,
                        tile_color = NA,
                        tile_color_NA = F,
                        tile_text = F,
                        tile_line = F,
                        line_args = list(linewidth = 0.1, color = "black"),
                        base_theme = colrr::theme_material,
                        base_theme_args = list(),
                        theme_args = list(panel.grid = ggplot2::element_blank()),
                        pattern_lim_size = 0,
                        pattern_lim_pos = c("inner", "outer"),
                        pattern_names = 0, # make numeric, zero = disabled
                        pattern_names_fun = ggrepel::geom_text_repel,
                        pattern_names_fun_args = list(),
                        pairwiseAlignment = NULL,
                        subject_lim_lines = F,
                        subject_name = NULL,
                        subject_name_infer = T,
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        start_end_col = "start_end",
                        y_group_col = NULL,
                        coord_fixed_ratio = NULL,
                        x_breaks = NULL,
                        algnmt_type = NULL,
                        add_length_suffix = F,
                        group_on_yaxis = F,
                        min_gap = 10,
                        ref = NULL,
                        pos_shift = NULL,
                        pos_shift_adjust_axis = T,
                        verbose = T,
                        focus = NULL,
                        y_order = c("as_is", "increasing", "decreasing")) {

  # document and clean up

  ## check duplicate names

  if (!requireNamespace("Peptides", quietly = T)){
    utils::install.packages("Peptides")
  }
  if (!requireNamespace("farver", quietly = T)){
    utils::install.packages("farver")
  }
  if (!requireNamespace("viridisLite", quietly = T)){
    utils::install.packages("viridisLite")
  }

  pattern_lim_pos <- match.arg(pattern_lim_pos, c("inner", "outer"))
  pattern_names_fun <- match.fun(pattern_names_fun)
  y_order <- match.arg(y_order, c("as_is", "increasing", "decreasing"))
  base_theme <- match.fun(base_theme)

  ## checks
  if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet") || methods::is(algnmt, "AAStringSet")) {
    # from DECIPHER::AlignSeqs
    if (methods::is(algnmt, "DNAStringSet") || methods::is(algnmt, "RNAStringSet")) {
      algnmt_type <- "NT"
    }
    if (methods::is(algnmt, "AAStringSet")) {
      algnmt_type <- "AA"
    }
    algnmt <- xstringset_to_df(xstringset = algnmt,
                               name_col = name_col,
                               seq_col = seq_col,
                               pos_col = pos_col)
    if (!is.null(y_group_col)) {
      ## y_group_col makes only sense when a algnmt is a dataframe, e.g. from multi pairwisealignment
      if (verbose) {
        message("y_group_col is set to NULL.")
      }
      y_group_col <- NULL
    }
  } else if (methods::is(algnmt, "PairwiseAlignmentsSingleSubject") || methods::is(algnmt, "list")) {
    # from pwalign::pairwiseAlignment
    if (methods::is(algnmt, "PairwiseAlignmentsSingleSubject")) {
      algnmt_type <- ifelse(guess_type2(algnmt) == "AA", "AA", "NT")
    } else if (methods::is(algnmt, "list")) {
      algnmt_type <- ifelse(guess_type2(algnmt[[1]]) == "AA", "AA", "NT")
    }

    algnmt <- pwalign_to_df(
      pa = algnmt,
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
  # else: data frame from igsc::pwalign_multi

  if (!name_col %in% names(algnmt)) {
    stop("name_col not found in algnmt.")
  }
  if (!seq_col %in% names(algnmt)) {
    stop("seq_col not found in algnmt.")
  }
  if (!pos_col %in% names(algnmt)) {
    stop("pos_col not found in algnmt.")
  }
  if (!is.null(y_group_col) && !y_group_col %in% names(algnmt)) {
    stop("y_group_col not found in algnmt.")
  }

  if (!is.null(pairwiseAlignment) && is.null(subject_name)) {
    if (!is.null(pairwiseAlignment@subject@unaligned@ranges@NAMES)) {
      subject_name <- pairwiseAlignment@subject@unaligned@ranges@NAMES
    }
  }

  infer_subject_name(algnmt = algnmt,
                     seq_col = seq_col,
                     name_col = name_col,
                     pos_col = pos_col,
                     subject_name = subject_name,
                     subject_name_infer = subject_name_infer)

  if (subject_lim_lines && !is.null(subject_name) && !subject_name %in% algnmt[,name_col,drop=T]) {
    stop("subject_name not found in ", name_col, " of algnmt. Can't draw subject_lim_lines.")
  }

  if (!is.null(focus) && !is.null(subject_name)) {
    if (!is.numeric(focus)) {
      stop("focus should be a positive integer.")
    }
    pos_col_rng <- range(algnmt[which(algnmt[,name_col,drop=T] != subject_name), pos_col])
    algnmt <-
      algnmt %>%
      dplyr::filter(dplyr::between(!!rlang::sym(pos_col), pos_col_rng[1]-focus, pos_col_rng[2]+focus))
  } else if (!is.null(focus) && is.null(subject_name)) {
    message("focus requires subject_name which is NULL though.")
  } else if (!is.null(focus) && !is.null(subject_name) && !subject_name %in% algnmt[,name_col,drop=T]) {
    message("subject name not found in ", name_col, ". Cannot use focus.")
  }

  if (!is.null(pos_shift)) {
    if (!is.numeric(pos_shift) && !grepl("^\\+", pos_shift)) {
      stop("pos_shift has to be numeric giving and absolute position as start or start with a '+' giving a relative shift.")
    }
    if (grepl("^\\+", pos_shift)) {
      pos_shift <- max(algnmt[,pos_col,drop=T]) - as.numeric(gsub("\\+", "", pos_shift)) + 2
      if (pos_shift < 0) {
        stop("pos_shift cannot be larger than the largest alignment position.")
      }
    }
    conv <- stats::setNames(shifted_pos(x = unique(algnmt[,pos_col,drop=T]), start_pos = pos_shift, verbose = verbose), unique(algnmt[,pos_col,drop=T]))
    algnmt[,pos_col] <- conv[algnmt[,pos_col,drop=T]]
  }

  if (group_on_yaxis && is.null(y_group_col)) {
    # if pairwise alignment is provided, use this one to derived ranges - but this may not be compatible with shifting?!
    if (is.null(subject_name)) {
      message("Cannot group on y-axis without subject name.")
    } else {
      subject.ranges <- subject.ranges[which(names(subject.ranges) != subject_name)]
      subject.ranges <- subject.ranges[order(sapply(subject.ranges, min))]
      rows <- numeric(length(subject.ranges))
      rows[1] <- 1
      # priority to row 1, or lowest row in general
      for (i in seq_along(subject.ranges)[-1]) {
        row_set <- 1
        j <- which(rows == row_set) # consider gapped patterns; j are the indices of subject.ranges to consider
        overlap <- T
        while (overlap) {
          # loop through all pattern of that row and check if there a larger overlap than allowed (or gap smaller than allowed)
          # what if min_gap are negative values? and what if min_gap is larger than the current range - select the larger of range or min_gap
          # as subject.ranges are ordered: check for every pattern or just the most recent one? - max(j) does that
          # how to account for gapped alignment - one pattern being within a gap of another pattern - like below.
          if (any(sapply(j, function(x) length(intersect(subject.ranges[[x]], subject.ranges[[i]])) > min_gap))) {
            row_set <- row_set + 1
            j <- which(rows == row_set) # which patterns are in next row to consider
            if (length(j) == 0) {
              # first pattern in that row - no further checking needed
              overlap <- F
            }
          } else {
            overlap <- F
          }
        }
        rows[i] <- row_set
      }

      rows <- paste0("Group_", rows)
      names(rows) <- names(subject.ranges)
      rows <- c(stats::setNames(subject_name, subject_name), rows)
      y_group_col <- "group"
      algnmt[,y_group_col] <- rows[as.character(unname(algnmt[,name_col,drop=T]))]
      algnmt[,y_group_col] <- factor(algnmt[,y_group_col,drop=T], levels = c(subject_name, unique(rows[-1])))
    }

  } else if (group_on_yaxis && !is.null(y_group_col)) {
    if (verbose) {
      message("y_group_col is not NULL. Using this one. Ignoring group_on_yaxis.")
    }
  }

  if (y_order != "as_is") {
    # if (is.null(y_group_col) && is.factor(algnmt[[name_col]])) {
    #   message("name_col is a factor. Will stick to this order.")
    #   order <- F
    # } else if (!is.null(y_group_col) && is.factor(algnmt[[y_group_col]])) {
    #   message("y_group_col is a factor. Will stick to this order.")
    #   order <- F
    # }
    if (is.null(y_group_col)) {
      # second level of ordering: length: shortest first
      yorder <- data.frame(name = names(subject.ranges),
                           minpos = sapply(subject.ranges, min, simplify = T),
                           len = sapply(subject.ranges, length, simplify = T)) %>%
        dplyr::arrange(minpos, len) %>%
        dplyr::pull(name)
      if (y_order == "decreasing") {
        yorder <- rev(yorder)
      }
      if (!is.null(subject_name)) {
        yorder <- unique(c(subject_name, yorder))
      }
      algnmt[[name_col]] <- factor(algnmt[[name_col]], levels = yorder)
    } else {
      yorder <-
        algnmt %>%
        dplyr::group_by(!!rlang::sym(y_group_col)) %>%
        dplyr::summarise(min = min(!!rlang::sym(pos_col)))
      if (!is.null(subject_name)) {
        subject_name2 <-
          algnmt %>%
          dplyr::filter(!!rlang::sym(name_col) == subject_name) %>%
          dplyr::distinct(!!rlang::sym(y_group_col)) %>%
          dplyr::pull(!!rlang::sym(y_group_col))
        # yorder <- dplyr::filter(yorder, !!rlang::sym(y_group_col) != subject_name2) # unique does the job
        yorder <- unique(c(subject_name2, names(sort(stats::setNames(yorder$min, yorder[,1,drop=T]), decreasing = y_order == "decreasing"))))
      } else {
        yorder <- names(sort(stats::setNames(yorder$min, yorder[,1,drop=T]), decreasing = y_order == "decreasing"))
      }
      algnmt[[y_group_col]] <- factor(algnmt[[y_group_col]], levels = yorder)
    }
  }

  if (!is.null(ref)) {
    if (!ref %in% unique(algnmt[,name_col,drop=T])) {
      stop("ref not found in names of algnmt.")
    }
    seq_original <- "seq_original"
    algnmt <- compare_seq_df_long(df = algnmt,
                                  ref = ref,
                                  pos_col = pos_col,
                                  seq_col = seq_col,
                                  name_col = name_col,
                                  seq_original = seq_original,
                                  match_symbol = ".",
                                  mismatch_symbol = "x",
                                  change_pattern = T,
                                  pattern_mismatch_as = "mismatch_symbol",
                                  change_ref = F)

    # seq_original will cause that tiles are filled according to chemical property even if they are replaced by a dots when
    # sequence is equal to the reference sequence (ref)
  } else {
    seq_original <- seq_col
  }

  if (!all(c(pos_col, seq_col, name_col) %in% names(algnmt))) {
    stop("algnmt at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ". Alternatively change function arguments.")
  }

  # if seq names are numeric; order them increasingly
  if (!anyNA(suppressWarnings(as.numeric(as.character(algnmt[,name_col,drop=T]))))) {
    algnmt[,name_col] <- factor(algnmt[,name_col,drop=T], levels = as.character(unique(as.numeric(as.character(algnmt[,name_col,drop=T])))))
  } else if (!is.factor(algnmt[,name_col,drop=T])) {
    algnmt[,name_col] <- as.factor(algnmt[,name_col,drop=T])
  }

  if (is.null(algnmt_type)) {
    inds <- intersect(which(!is.na(algnmt[,seq_col,drop=T])), which(!algnmt[,seq_col,drop=T] %in% c(".", "x", "match", "mismatch", "gap", "insertion", "ambiguous")))
    algnmt_type <- guess_type(seq_vector = algnmt[,seq_col,drop=T][inds])
  }

  # use preset colors
  if (is.null(tile_fill)) {
    if (algnmt_type == "NT") {
      tile_fill <- colnames(igsc:::scheme_NT)[1]
    }
    if (algnmt_type == "AA") {
      tile_fill <- colnames(igsc:::scheme_AA)[1]
      # set tile_fill to Chemistry_AA by default which will cause falling into the special case below
      # to avoid that, change this here to tile_color (like in case of algnmt_type == "NT") and connect
    }
  }

  #default; may change when algnmt_type == "AA" && tile_fill == "Chemistry_AA"
  col_col <- seq_col
  if (length(tile_fill) == 1 && tile_fill %in% c(colnames(igsc:::scheme_AA),
                                                 colnames(igsc:::scheme_NT),
                                                 names(purrr::flatten(Peptides:::AAdata)))) {
    if (algnmt_type == "NT") {
      tile_fill_internal <- igsc:::scheme_NT[,match.arg(tile_fill, choices = colnames(igsc:::scheme_NT)),drop=T]
    } else if (algnmt_type == "AA") {
      if (tile_fill == "Chemistry_AA") {
        # special case; change legend to chemical property
        col_col <- tile_fill
        chem_col <- utils::stack(igsc:::aa_info[["aa_main_prop"]])
        names(chem_col) <- c(col_col, seq_original)
        chem_col[,seq_original] <- as.character(chem_col[,seq_original,drop=T])
        algnmt <- dplyr::left_join(algnmt, chem_col, by = seq_original)
        algnmt[,col_col] <- ifelse(is.na(algnmt[,col_col,drop=T]), algnmt[,seq_original,drop=T], algnmt[,col_col,drop=T])
        tile_fill_internal <- igsc:::scheme_AA[,match.arg(tile_fill, choices = colnames(igsc:::scheme_AA)),drop=T]
        tile_fill_internal <- tile_fill_internal[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
        for (i in 1:length(tile_fill_internal)) {
          if (names(tile_fill_internal)[i] %in% names(igsc:::aa_info[["aa_main_prop"]])) {
            names(tile_fill_internal)[i] <- igsc:::aa_info[["aa_main_prop"]][names(tile_fill_internal)[i]]
          }
        }
        tile_fill_internal <- tile_fill_internal[unique(names(tile_fill_internal))]
      } else if (tile_fill %in% names(purrr::flatten(Peptides:::AAdata))) {
        col_col <- tile_fill
        algnmt[,col_col] <- purrr::flatten(Peptides:::AAdata)[[tile_fill]][algnmt[,seq_original,drop=T]]
        tile_fill_internal <- viridisLite::viridis(100)[cut(unique(algnmt[,col_col,drop=T])[which(!is.na(unique(algnmt[,col_col,drop=T])))], 100)]
        names(tile_fill_internal) <- unique(algnmt[,col_col,drop=T])[which(!is.na(unique(algnmt[,col_col,drop=T])))]
      } else {
        tile_fill_internal <- igsc:::scheme_AA[,match.arg(tile_fill, choices = colnames(igsc:::scheme_AA)),drop=T]
      }
    } else {
      if (verbose) {
        message("Type of alignment data (NT or AA) could not be determined. Choosing default ggplot colors.")
      }
      tile_fill_internal <- scales::hue_pal()(length(unique(as.character(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))))
    }
  } else {
    # provide own colors
    if (is.null(tile_fill)) {
      # happens when algmnt_type is not NT or AA
      # redundant to above when NT or AA cannot be guessed
      tile_fill <- scales::hue_pal()(length(unique(as.character(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))))
    } else {
      if (is.null(names(tile_fill))) {
        if (length(tile_fill) < length(unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]))) {
          warning("Less colors provided than entities in the alignment.")
        }
      } else {
        if (any(!unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))]) %in% names(tile_fill))) {
          warning("Not all values in algnmt[,seq_col] found in names(tile_fill).")
        }
      }
    }
    tile_fill_internal <- tile_fill
  }

  if (!is.null(names(tile_fill_internal)) && col_col == seq_col) {
    tile_fill_internal <- tile_fill_internal[unique(algnmt[,seq_col,drop=T][which(!is.na(algnmt[,seq_col,drop=T]))])]
  }

  # make sure it is not a factor
  if (!is.numeric(algnmt[,pos_col,drop=T])) {
    algnmt[,pos_col] <- as.numeric(as.character(algnmt[,pos_col,drop=T]))
  }
  algnmt[,seq_col] <- as.character(algnmt[,seq_col,drop=T])

  # https://stackoverflow.com/questions/45493163/ggplot-remove-na-factor-level-in-legend
  # --> na.translate = F
  if (is.null(x_breaks)) {
    ## fix with expand = F in coord_cartesian
    x_breaks <- floor(base::pretty(unique(algnmt[[pos_col]]), n = 5))
    if (!any(x_breaks < 0)) {
      x_breaks <- x_breaks[which(x_breaks > 0)]
    }
    #x_breaks <- x_breaks[which(x_breaks < max(algnmt[,pos_col,drop=T]))]
  }

  if (algnmt_type == "AA" && length(tile_fill) == 1 && tile_fill == "Chemistry_AA") {
    col_breaks = unique(igsc:::aa_info[["aa_main_prop"]])[which(unique(igsc:::aa_info[["aa_main_prop"]]) != "stop")]
  } else {
    col_breaks = ggplot2::waiver()
  }
  yaxis <- ifelse(is.null(y_group_col), rlang::sym(name_col), rlang::sym(y_group_col))


  ## algnmt summary prep
  ## exclude subject ??
  algnmt_summary <-
    algnmt %>%
    dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>%
    dplyr::group_by(!!rlang::sym(name_col)) %>%
    dplyr::summarise(n_not_NA = sum(!is.na(!!rlang::sym(seq_col))), min_pos = min(!!rlang::sym(pos_col)), max_pos = max(!!rlang::sym(pos_col))) %>%
    dplyr::mutate(mid_pos = max_pos - (max_pos-min_pos)/2)

  if (subject_lim_lines && is.null(subject_name)) {
    if (verbose) {
      message("subject_name is NULL or could not be inferred subject_lim_lines set to FALSE.")
    }
    subject_lim_lines <- F
  }

  if (start_end_col %in% names(algnmt)) {
    #seq_original may come from compare_seqs_df, when !is.null(ref)
    algnmt_summary <-
      algnmt_summary %>%
      dplyr::left_join(algnmt %>%
                         dplyr::filter(!is.na(!!rlang::sym(start_end_col))) %>%
                         dplyr::select(!!!rlang::syms(names(algnmt)[which(!names(algnmt) %in% c(seq_col, "seq_original"))])) %>%
                         tidyr::pivot_wider(names_from = !!rlang::sym(start_end_col), values_from = !!rlang::sym(pos_col)),
                       by = name_col)
  }

  if (!is.null(y_group_col)) {
    test_dups <-
      algnmt %>%
      tidyr::drop_na(!!rlang::sym(seq_col)) %>% # y_group_col
      dplyr::group_by(!!rlang::sym(y_group_col), !!rlang::sym(pos_col)) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    if (any(test_dups$n > 1)) {
      message("Your y_group_col causes multiple sequnces (subject/patterns) to occupy the same position in the plot. You may want to check on that.")
      message("algnmt %>%
  tidyr::drop_na(!!rlang::sym(y_group_col)) %>%
  dplyr::group_by(!!rlang::sym(y_group_col), !!rlang::sym(pos_col)) %>%
  dplyr::summarise(n = dplyr::n())")
    }
    y_group_conv <- dplyr::distinct(tidyr::drop_na(algnmt, !!rlang::sym(seq_col)), !!rlang::sym(y_group_col), !!rlang::sym(name_col))
    y_group_conv <- stats::setNames(as.character(y_group_conv[,y_group_col,drop=T]), y_group_conv[,name_col,drop=T])
    algnmt_summary[[yaxis]] <- y_group_conv[as.character(algnmt_summary[,name_col,drop=T])]
    if (is.factor(algnmt[,yaxis,drop=T])) {
      algnmt_summary[,yaxis] <- factor(algnmt_summary[,yaxis,drop=T], levels = levels(algnmt[,yaxis,drop=T]))
    } else {
      algnmt_summary[,yaxis] <- factor(algnmt_summary[,yaxis,drop=T], levels = unique(algnmt_summary[,yaxis,drop=T]))
    }
  }

  if (add_length_suffix) {
    # keep this, with pairwiseAlignment
    if (is.null(pairwiseAlignment)) {
      if (verbose) {
        message("pairwiseAlignment is not provided and orignal sequences are unknown. Length suffix is inferred from number of non-NA elements in seq_col.")
      }
      seq_lengths <- stats::setNames(algnmt_summary[,"n_not_NA",drop=T], algnmt_summary[,name_col,drop=T])
    } else {
      if (verbose) {
        message("Length suffix is inferred from pairwiseAlignment.")
      }
      seq_lengths <- stats::setNames(c(pairwiseAlignment@subject@unaligned@ranges@width,
                                       pairwiseAlignment@pattern@unaligned@ranges@width),
                                     c(pairwiseAlignment@subject@unaligned@ranges@NAMES,
                                       pairwiseAlignment@pattern@unaligned@ranges@NAMES))
    }
    if (algnmt_type %in% c("NT", "AA")) {
      seq_lengths <- stats::setNames(paste0(names(seq_lengths), "\n", seq_lengths, " ", tolower(algnmt_type)), nm = names(seq_lengths))
    } else {
      seq_lengths <- stats::setNames(paste0(names(seq_lengths), "\n", seq_lengths), nm = names(seq_lengths))
    }

    algnmt_summary$label <- seq_lengths[algnmt_summary[,name_col,drop=T]]
  } else {
    algnmt_summary$label <- stats::setNames(as.character(algnmt_summary[, name_col,drop=T]), as.character(algnmt_summary[, name_col,drop=T]))
  }

  plot <-
    ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(yaxis))) + # name_col
    Gmisc::fastDoCall(base_theme, args = base_theme_args) +
    Gmisc::fastDoCall(ggplot2::theme, args = theme_args) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::coord_cartesian(expand = FALSE)

  if (add_length_suffix && is.null(y_group_col)) {
    # nt can only be added to y-axis text when sequences are not grouped
    plot <- plot + ggplot2::scale_y_discrete(labels = algnmt_summary[,"label",drop=T][levels(algnmt[,name_col,drop=T])])
  } else if (add_length_suffix && !is.null(y_group_col)) {
    message("length suffices not plotable when y-axis is grouped.")
  }

  if (is.numeric(algnmt[,col_col,drop=T])) {
    plot <- plot + ggplot2::scale_fill_viridis_c(na.value = "white")
  } else {
    plot <- plot + ggplot2::scale_fill_manual(values = tile_fill_internal, breaks = col_breaks, na.value = "white", na.translate = F)
  }

  if (tile_line) {
    if (!start_end_col %in% names(algnmt)) {
      if (verbose) {
        message("start_end_col not found in algnmt data frame. Using min and max position to draw line. But this could be wrong if, e.g. the first exon is at later position as the first one.")
      }
      xmin <- "min_pos"
      xmax <- "max_pos"
    } else {
      xmin <- "start"
      xmax <- "end"
    }

    if (start_end_col %in% names(algnmt) && any(algnmt_summary$start > algnmt_summary$end)) {
      # split the plotting of line_segment
      # start to very last position
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(algnmt_summary, start > end),
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = max(algnmt[[pos_col]]), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # very first position to very end
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(algnmt_summary, start > end),
                                                                            ggplot2::aes(x = 1, xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # plot those where start < end
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = dplyr::filter(algnmt_summary, start < end),
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
      # reordering of axis needed as step-wise plotting of lines corrupts the original order
      plot <- suppressMessages(plot + ggplot2::scale_y_discrete(limits = levels(algnmt[[yaxis]])))
    } else {
      # all starts are smaller than ends
      plot <- plot + Gmisc::fastDoCall(ggplot2::geom_segment, args = c(line_args,
                                                                       list(data = algnmt_summary,
                                                                            ggplot2::aes(x = !!rlang::sym(xmin), xend = !!rlang::sym(xmax), y = !!rlang::sym(yaxis), yend = !!rlang::sym(yaxis)))))
    }
  }


  if (!is.na(tile_color)) {
    plot <- plot + ggplot2::geom_tile(data = if(tile_color_NA) {algnmt} else {algnmt[which(!is.na(algnmt[,seq_col,drop=T])),]},
                                      color = tile_color,
                                      ggplot2::aes(fill = !!rlang::sym(col_col)))
    # geom_raster seems not take color?!
  } else {
    plot <- plot + ggplot2::geom_tile(data = algnmt[which(!is.na(algnmt[,seq_col,drop=T])),],
                                      ggplot2::aes(fill = !!rlang::sym(col_col)))
    # geom_raster # geom_tile
  }

  if ((is.logical(tile_text) && tile_text) || (is.numeric(tile_text) && tile_text > 0)) {
    if (is.logical(tile_text)) {
      text_size <- 4
    } else {
      text_size <- tile_text
    }

    bckgr_colors <- farver::decode_colour(tile_fill_internal[as.character(algnmt[,col_col,drop=T])], to = "hcl")
    algnmt$text_colors <- ifelse(bckgr_colors[, "l"] > 50, "black", "white")

    plot <-
      plot +
      ggplot2::geom_text(data = algnmt, ggplot2::aes(label = !!rlang::sym(seq_col), color = I(text_colors)),
                         na.rm = T, size = text_size)
  }

  if (!is.null(coord_fixed_ratio)) {
    plot <- plot + ggplot2::coord_fixed(ratio = coord_fixed_ratio)
  }

  if (!is.null(pos_shift) && pos_shift_adjust_axis) {
    plot_x_breaks <- ggplot2::ggplot_build(plot)[["layout"]][["panel_params"]][[1]][["x"]][["breaks"]]
    plot_x_breaks <- plot_x_breaks[which(plot_x_breaks <= pos_shift)]
    plot <- suppressMessages(plot + ggplot2::scale_x_continuous(breaks = max(algnmt[[pos_col]]) - pos_shift + 1 + c(1, plot_x_breaks),
                                                                labels = c(1, plot_x_breaks)))
  }

  if (pattern_lim_size > 0 && !is.null(pairwiseAlignment)) {
    pattern.ranges <- data.frame(pairwiseAlignment@pattern@range,
                                 seq.name = ifelse(rep(is.null(pairwiseAlignment@pattern@unaligned@ranges@NAMES), length(pairwiseAlignment)),
                                                   paste0("pattern_", seq(1,length(pairwiseAlignment))),
                                                   pairwiseAlignment@pattern@unaligned@ranges@NAMES))
    names(pattern.ranges)[1:2] <- paste0("pattern_", names(pattern.ranges)[1:2])
    # stats::setNames(as.data.frame(pairwiseAlignment@subject@range)[1:2], nm = paste0("subject_", names(as.data.frame(pairwiseAlignment@subject@range)[1:2])))
    # subject position not from pairwiseAlignment but from algnmt_summary as the latter can account for pos_shift
    ranges <- dplyr::left_join(pattern.ranges, algnmt_summary, by = dplyr::join_by(!!rlang::sym(name_col)))
    ranges <- tidyr::pivot_longer(ranges, cols = c(min_pos, max_pos), names_to = "pos", values_to = "inner_position")
    ranges$outer_position <- ifelse(ranges$pos == "min_pos", -2, max(algnmt[,pos_col,drop=T]) + 2)
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
      message("pattern limits can only be plotted if pairwiseAlignment is provided.")
    }
  }


  if (pattern_names > 0) {
    ## ggrepel considers the pattern.lim labels - great.
    plot <- plot + Gmisc::fastDoCall(pattern_names_fun, args = c(list(data = if(is.null(subject_name)) {algnmt_summary} else {algnmt_summary %>% dplyr::filter(!!rlang::sym(name_col) != subject_name)},
                                                                      mapping = ggplot2::aes(x = mid_pos, y = !!rlang::sym(yaxis), label = label),
                                                                      size = pattern_names,
                                                                      inherit.aes = F),
                                                                 pattern_names_fun_args))
  }

  if (subject_lim_lines) {
    #minmaxpos <- c(algnmt_summary %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::pull(min_pos), algnmt_summary %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::pull(max_pos))
    min.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_min(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
    max.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_max(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
    plot <- plot + ggplot2::geom_vline(xintercept = c(min.pos, max.pos), linetype = "dashed")
  }


  return(plot)
}


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

infer_subject_name <- function(algnmt,
                               seq_col,
                               name_col,
                               pos_col,
                               subject_name,
                               subject_name_infer,
                               verbose = T) {

  subject.ranges <- algnmt %>% tidyr::drop_na(!!rlang::sym(seq_col)) %>% dplyr::group_by(!!rlang::sym(name_col))
  subject.ranges <- stats::setNames(igsc:::seq2(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col)),
                                                subject.ranges %>% dplyr::slice_max(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col))),
                                    nm = as.character(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(name_col))))
  assign("subject.ranges", subject.ranges, envir = parent.frame())


  if (is.null(subject_name)) {
    # just try, possible from pa_to_df
    # remains NULL if does not exist
    subject_name <- attr(algnmt, "subject_name")
  }

  if (is.null(subject_name) && subject_name_infer) {
    max_len <- max(lengths(subject.ranges))
    if (length(which(lengths(subject.ranges) == max_len)) > 1) {
      if (verbose) {
        message("Subject could not be identified.")
      }
      ## TODO: set variable for ordering pattern on y-axis to FALSE here
    } else {
      subject_name <- names(which(lengths(subject.ranges) == max_len))
    }
    # assigns in parent environment (https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r?rq=1)
    assign("subject_name", subject_name, envir = parent.frame())
    message("subject sequence inferred: ", subject_name)
  }
}
