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
#' @param pa list of pairwise alignments in form of PairwiseAlignmentsSingleSubject
#' or as actual list of single pairwise alignments; the latter will sometime make the alignment quicker
#' rather for internal use by igsc::MultiplePairwiseAlignmentsToOneSubject;
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
#' @param ref a reference sequence to compare all other sequnces to; e.g. a consensus sequence
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
                        pa = NULL,
                        subject.lim.lines = F,
                        pos_col = "position",
                        seq_col = "seq",
                        name_col = "seq.name",
                        coord_fixed_ratio = NULL,
                        x.breaks = NULL,
                        algnmt_type = NULL,
                        ref = NULL) {
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
    if (methods::is(pa, "PairwiseAlignmentsSingleSubject")) {
      temp <- pa[1]@pattern@unaligned
    } else if (methods::is(pa, "list")) {
      temp <- pa[[1]]@pattern@unaligned
    }
    if (methods::is(temp, "QualityScaledDNAStringSet") || methods::is(temp, "QualityScaledRNAStringSet")) {
      algnmt_type <- "NT"
    }
    if (methods::is(temp, "QualityScaledAAStringSet")) {
      algnmt_type <- "AA"
    }

    algnmt <- pa_to_df(algnmt)
  }
  # else: data frame from igsc::MultiplePairwiseAlignmentsToOneSubject

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


  if (!all(c(pos_col, seq_col, name_col) %in% names(algnmt))) {
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

  plot <-
    ggplot2::ggplot(algnmt, ggplot2::aes(x = !!rlang::sym(pos_col), y = !!rlang::sym(name_col))) +
    theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   text = ggplot2::element_text(family = font.family)) +
    ggplot2::scale_x_continuous(breaks = x.breaks)

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
    plot <- plot + ggplot2::geom_raster(data = algnmt[which(!is.na(algnmt[,seq_col,drop=T])),],
                                        ggplot2::aes(fill = !!rlang::sym(col_col)))
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

  if (pattern.lim.size > 0 && !is.null(pa)) {
    if (!methods::is(pa, "list")) {
      pa <- as.list(pa)
    }
    pattern.ranges <- purrr::map_dfr(pa, function(x) data.frame(x@pattern@range, seq.name = ""))
    if (is.null(names(pa))) {
      pattern.ranges$seq.name <- make.names(purrr::map(pa, function(x) names(Biostrings::alignedPattern(x))))
    } else {
      pattern.ranges$seq.name <- names(pa)
    }
    #for (i in 1:length(pa)) {pattern.ranges[i, "seq.name"] <- make.names(names(Biostrings::alignedPattern(pa[i])))}
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
