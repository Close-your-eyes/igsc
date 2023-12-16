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
#' @param theme ggplot theme to use as basis
#' @param pattern_lim_size numeric; plot annotation of pattern alignment limits; value indicates size;
#' set to 0 to omit plotting; only applicable if pa is provided
#' @param subject_lim_lines logical whether to plot vertical lines of subject alignment limits
#' @param pos_col name of position column in algnmt (applicable if algnmt is a data.frame)
#' @param seq_col name of sequence column in algnmt (applicable if algnmt is a data.frame)
#' @param name_col name of column which holds sequence names in algnmt (applicable if algnmt is a data.frame)
#' @param coord_fixed_ratio numeric; fixed aspect ratio of plot; leave NULL to not force a ratio
#' @param x_breaks numeric vector; manually provide breaks (ticks) on x-axis; set NULL to have breaks
#' picked automatically
#' @param algnmt_type 'NT' or 'AA'; only required if algnmt is a data.frame and only if
#' NT or AA cannot be guessed; leave NULL to have it guessed based on data
#' @param ref a reference sequence to compare all other sequnces to; e.g. a consensus sequence
#' @param subject_name name of the subject sequence in algnmt; only needed if subject_lim_lines = TRUE
#' @param pattern_lim_pos
#' @param pattern_names
#' @param pairwiseAlignment
#' @param y_group_col
#' @param pattern_names_fun
#' @param add_length_suffix
#' @param group_on_yaxis
#' @param min_gap
#' @param line
#' @param line_args
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
                        line = F,
                        line_args = list(linewidth = 0.2, color = "black"),
                        theme = ggplot2::theme_bw(),
                        theme_args = list(panel.grid = ggplot2::element_blank()),
                        pattern_lim_size = 0,
                        pattern_lim_pos = "inner",
                        pattern_names = F,
                        pattern_names_fun = ggrepel::geom_text_repel, # ggplot2::geom_text
                        pairwiseAlignment = NULL,
                        subject_lim_lines = F,
                        subject_name = NULL,
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        start_end_col = "start_end",
                        y_group_col = NULL, # set groups of patterns to plot on same y-axis level
                        coord_fixed_ratio = NULL,
                        x_breaks = NULL,
                        algnmt_type = NULL,
                        add_length_suffix = F,
                        group_on_yaxis = F,
                        min_gap = 10,
                        ref = NULL,
                        pos_shift = NULL, # fixed position like: 100000 or relative: "+100"
                        pos_shift_adjust_axis = T) {
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

  if (subject_lim_lines && !is.null(subject_name) && !subject_name %in% algnmt[,name_col,drop=T]) {
    stop("subject_name not found in ", name_col, " of algnmt.")
  }

  if (!is.null(pos_shift)) {
    if (!is.numeric(pos_shift) && !grepl("^\\+", pos_shift)) {
      stop("pos_shift has to be numeric giving and absolute position as start or start with a '+' giving a relative shift.")
    }
    if (grepl("^\\+", pos_shift)) {
      pos_shift <- max(algnmt$position)-as.numeric(gsub("\\+", "", pos_shift)) + 2
      if (pos_shift < 0) {
        stop("pos_shift cannot be larger than the largest alignment position.")
      }
    }
    conv <- setNames(shifted_pos(x = unique(algnmt$position), start_pos = pos_shift), unique(algnmt$position))
    algnmt$position <- conv[algnmt$position]
  }

  pattern_lim_pos <- match.arg(pattern_lim_pos, c("inner", "outer"))
  pattern_names_fun <- match.fun(pattern_names_fun)

  if (!is.null(pairwiseAlignment) && is.null(subject_name)) {
    if (!is.null(pairwiseAlignment@subject@unaligned@ranges@NAMES)) {
      subject_name <- pairwiseAlignment@subject@unaligned@ranges@NAMES
    }
  }

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
    # if pairwise alignment is provided, use this one to derived ranges - but this may not be compatible with shifting?!
    # remove this procedure - too complicated
   ' if (!is.null(pairwiseAlignment)) {
      subject.ranges <- seq2(pairwiseAlignment@subject@range@start, pairwiseAlignment@subject@range@start+pairwiseAlignment@subject@range@width-1)
      names(subject.ranges) <- pairwiseAlignment@pattern@unaligned@ranges@NAMES
      subject_name <- setdiff(unique(algnmt[,name_col,drop=T]), names(subject.ranges))
    } else {
      # here the subject has to be separated from patterns
      subject.ranges <- algnmt %>% tidyr::drop_na(!!rlang::sym(seq_col)) %>% dplyr::group_by(!!rlang::sym(name_col))
      subject.ranges <- stats::setNames(seq2(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col)),
                                             subject.ranges %>% dplyr::slice_max(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col))),
                                        nm = as.character(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(name_col))))
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
    }'
    # here the subject has to be separated from patterns
    subject.ranges <- algnmt %>% tidyr::drop_na(!!rlang::sym(seq_col)) %>% dplyr::group_by(!!rlang::sym(name_col))
    subject.ranges <- stats::setNames(seq2(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col)),
                                           subject.ranges %>% dplyr::slice_max(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(pos_col))),
                                      nm = as.character(subject.ranges %>% dplyr::slice_min(!!rlang::sym(pos_col)) %>% dplyr::pull(!!rlang::sym(name_col))))
    if (is.null(subject_name)) {
      max_len <- max(lengths(subject.ranges))
      if (length(which(lengths(subject.ranges) == max_len)) > 1) {
        message("Subject could not be identified as there are min. 2 sequences which have the max length. Cannot order patterns on y-axis. Provide is as argument 'subject_name'.")
        ## TODO: set variable for ordering pattern on y-axis to FALSE here
      } else {
        subject_name <- names(which(lengths(subject.ranges) == max_len))
      }
    }
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
    stop("algnmt at least has to contain columns named: ", pos_col, ", ", seq_col, ", ", name_col, ", ", y_group_col, ". Alternatively change function arguments.")
  }

  # if seq names are numeric; order them increasingly
  if (!anyNA(suppressWarnings(as.numeric(as.character(algnmt[,name_col,drop=T]))))) {
    algnmt[,name_col] <- factor(algnmt[,name_col,drop=T], levels = as.character(unique(as.numeric(as.character(algnmt[,name_col,drop=T])))))
  } else if (!is.factor(algnmt[,name_col,drop=T])) {
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
  if (is.null(x_breaks)) {
    x_breaks <- floor(base::pretty(algnmt[,pos_col,drop=T], n = 5))
    x_breaks <- x_breaks[which(x_breaks > 0)]
    x_breaks <- x_breaks[which(x_breaks < max(algnmt[,pos_col,drop=T]))]
  }

  if (algnmt_type == "AA" && length(color_values) == 1 && color_values == "Chemistry_AA") {
    col_breaks = unique(igsc:::aa_info[["aa_main_prop"]])[which(unique(igsc:::aa_info[["aa_main_prop"]]) != "stop")]
  } else {
    col_breaks = ggplot2::waiver()
  }
  yaxis <- ifelse(is.null(y_group_col), rlang::sym(name_col), rlang::sym(y_group_col))


  ## algnmt summary prep
  algnmt_summary <-
    algnmt %>%
    dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>%
    dplyr::group_by(!!rlang::sym(name_col)) %>%
    dplyr::summarise(n_not_NA = sum(!is.na(!!rlang::sym(seq_col))), min_pos = min(!!rlang::sym(pos_col)), max_pos = max(!!rlang::sym(pos_col))) %>%
    dplyr::mutate(mid_pos = max_pos - (max_pos-min_pos)/2)

  if (subject_lim_lines && is.null(subject_name)) {
    subject_guess <- algnmt_summary %>% dplyr::slice_max(n_not_NA, n = 1, with_ties = T)
    if (nrow(subject_guess) > 1) {
      message("subject_name is NULL or could not be guessed. subject_lim_lines set to FALSE.")
      subject_lim_lines <- F
    } else {
      subject_name <- as.character(subject_guess[,name_col,drop=T])
      message("subject guessed: ", subject_name)
    }
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
    y_group_conv <- dplyr::distinct(algnmt, !!rlang::sym(y_group_col), !!rlang::sym(name_col))
    y_group_conv <- stats::setNames(as.character(y_group_conv[,y_group_col,drop=T]), y_group_conv[,name_col,drop=T])
    algnmt_summary[[yaxis]] <- y_group_conv[as.character(algnmt_summary[,name_col,drop=T])]
    algnmt_summary[,y_group_col] <- factor(algnmt_summary[,y_group_col,drop=T], levels = levels(algnmt[,y_group_col,drop=T]))
  }

  if (add_length_suffix) {
    # keep this, with pairwiseAlignment
    if (is.null(pairwiseAlignment)) {
      message("pairwiseAlignment is not provided and orignal sequences are unknown. Length suffix is inferred from number of non-NA elements in seq_col.")
      seq_lengths <- stats::setNames(algnmt_summary[,"n_not_NA",drop=T], algnmt_summary[,name_col,drop=T])
    } else {
      message("Length suffix is inferred from pairwiseAlignment.")
      seq_lengths <- stats::setNames(c(pairwiseAlignment@subject@unaligned@ranges@width,
                                       pairwiseAlignment@pattern@unaligned@ranges@width),
                                     c(pairwiseAlignment@subject@unaligned@ranges@NAMES,
                                       pairwiseAlignment@pattern@unaligned@ranges@NAMES))
    }
    seq_lengths <- stats::setNames(paste0(names(seq_lengths), "\n", seq_lengths, " ", tolower(algnmt_type)), nm = names(seq_lengths))
    algnmt_summary$label <- seq_lengths[algnmt_summary[,name_col,drop=T]]
  } else {
    algnmt_summary$label <- stats::setNames(as.character(algnmt_summary[, name_col,drop=T]), as.character(algnmt_summary[, name_col,drop=T]))
  }

  plot <-
    ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(yaxis))) + # name_col
    theme +
    Gmisc::fastDoCall(ggplot2::theme, args = theme_args) +
    ggplot2::scale_x_continuous(breaks = x_breaks)

  if (add_length_suffix && is.null(y_group_col)) {
    # nt can only be added to y-axis text when sequences are not grouped
    plot <- plot + ggplot2::scale_y_discrete(labels = algnmt_summary[,"label",drop=T][levels(algnmt[,name_col,drop=T])])
  }


  if (is.numeric(algnmt[,col_col,drop=T])) {
    plot <- plot + ggplot2::scale_fill_viridis_c(na.value = "white")
  } else {
    plot <- plot + ggplot2::scale_fill_manual(values = tile_color, breaks = col_breaks, na.value = "white", na.translate = F)
  }

  if (line) {
    if (!start_end_col %in% names(algnmt)) {
      message("start_end_col not found in algnmt data frame. Using min and max position to draw line. But this could be wrong if, e.g. the first exon is at later position as the first one.")
      xmin <- "min_pos"
      xmax <- "max_pos"
    } else {
      xmin <- "start"
      xmax <- "end"
    }

    if (any(algnmt_summary$start > algnmt_summary$end)) {
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
                         na.rm = T, size = text_size)
  }

  if (!is.null(coord_fixed_ratio)) {
    plot <- plot + ggplot2::coord_fixed(ratio = coord_fixed_ratio)
  }

  if (!is.null(pos_shift) && pos_shift_adjust_axis) {
    plot_x_breaks <- ggplot2::ggplot_build(plot)[["layout"]][["panel_params"]][[1]][["x"]][["breaks"]]
    plot_x_breaks <- plot_x_breaks[which(plot_x_breaks <= pos_shift)]
    plot <- suppressMessages(plot + ggplot2::scale_x_continuous(breaks = max(algnmt$position) - pos_shift + 1 + c(1, plot_x_breaks),
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
    message("pattern limits can only be plotted if pairwiseAlignment is provided.")
  }

  if (pattern_names) {
    ## ggrepel considers the pattern.lim labels - great.
    plot <- plot + pattern_names_fun(data = algnmt_summary,
                                     ggplot2::aes(x = mid_pos, y = !!rlang::sym(yaxis), label = label),
                                     size = pattern_lim_size,
                                     inherit.aes = F)
  }

  if (subject_lim_lines) {
    #minmaxpos <- c(algnmt_summary %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::pull(min_pos), algnmt_summary %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::pull(max_pos))
    min.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_min(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
    max.pos <- algnmt %>% dplyr::filter(!!rlang::sym(seq_col) != "-") %>% dplyr::filter(!is.na(!!rlang::sym(seq_col))) %>% dplyr::filter(!!rlang::sym(name_col) != subject_name) %>% dplyr::slice_max(order_by = !!rlang::sym(pos_col), n = 1) %>% dplyr::pull(!!rlang::sym(pos_col))
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

  temp_order <- levels(df[,name_col,drop=T])
  temp_cols <- dplyr::select(df, -!!rlang::sym(seq_col)) # save cols for joining below; if other cols are not excluded, pivotting does not work
  df <-
    df %>%
    dplyr::select(!!rlang::sym(pos_col), !!rlang::sym(seq_col), !!rlang::sym(name_col)) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(name_col), values_from = !!rlang::sym(seq_col)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(!names(.) %in% c(ref, pos_col))])) %>%
    dplyr::mutate(value = paste0(value, "_", ifelse(value == !!rlang::sym(ref), ".", value))) %>%
    dplyr::mutate({{ref}} := paste0(!!rlang::sym(ref), "_", !!rlang::sym(ref))) %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(.)[which(names(.) != pos_col)]), names_to = name_col) %>%
    tidyr::separate(value, into = c(seq_original, seq_col), sep = "_") %>%
    dplyr::mutate({{name_col}} := factor(!!rlang::sym(name_col), levels = temp_order))
  df <- dplyr::left_join(df, temp_cols, by = dplyr::join_by(!!rlang::sym(pos_col), !!rlang::sym(name_col)))
  df[,seq_col] <- ifelse(df[,seq_col,drop=T] == "NA", NA, df[,seq_col,drop=T])
  df[,seq_original] <- ifelse(df[,seq_original,drop=T] == "NA", NA, df[,seq_original,drop=T])
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

## work with start pos - then derive shift
shifted_pos <- function(x,
                        n = NULL,
                        start_pos = NULL) {


  if (!is.null(n) && is.na(n)) {
    stop("n is NA, not numeric.")
  }
  if (!is.null(start_pos) && is.na(start_pos)) {
    stop("start_pos is NA, not numeric.")
  }

  if (!is.null(start_pos)) {
    if (!is.null(n)) {
      message("n will be ignored as start_pos is provided.")
    }
    n <- which(start_pos == x) - 1
    if (length(n) == 0) {
      stop("start_pos not found in x.")
    }
  }
  x2 <- c(x[(length(x)-n+1):length(x)], dplyr::lag(x,n)[(n+1):length(x)])
  return(x2)
}
