#' Print pairwise alignments of DNA or AA to console or txt file
#'
#' @param pa pwalign::pariwiseAlignment or XStringSet
#' @param linewidth numeric, how many letters to print per line
#' @param match_dots should matching position be printed as dots?;
#' 'subject' to have dots in the subject, 'pattern' to
#' have them in the pattern
#' @param print_pos print the positions of pattern and subject above and below
#' respectively
#' @param print_pos_end print positions on the end of each line
#' @param use_align_starts if TRUE, the first positions are not 1 but refer to
#' position of alignment within the provided sequences
#' @param out_file path to a file where to print results to; if NULL results are
#' printed in console
#' @param col_out color the printed alignment?
#' @param extend_subject extend subject further than alignment limits (-,+)
#'
#' @return alignment in printed format in console or file
#' @export
#'
#' @importFrom rlang %||%
#' @importFrom zeallot %<-%
#'
#' @examples
#' granzymes <- c("GZMA","GZMB","GZMH","GZMK","GZMM")
#' out <- get_sequences_from_biomart(granzymes)
#'
#' gzma <- dplyr::filter(out, hgnc_symbol == "GZMA")
#' gzmb <- dplyr::filter(out, hgnc_symbol == "GZMB")
#' GZMA_CDS <- stats::setNames(gzma$coding, "GZMA_CDS")
#' GZMB_CDS <- stats::setNames(gzmb$coding, "GZMB_CDS")
#' aln6 <- pwalign::pairwiseAlignment(subject = Biostrings::DNAStringSet(GZMB_CDS),
#'                                    pattern = Biostrings::DNAStringSet(GZMA_CDS),
#'                                    type = "global")
#' pwalign_print(pwalign_to_xstringset(aln6), match_dots = "subject", col_out = T)
#' pwalign_print(aln6, match_dots = "subject", col_out = T)
pwalign_print <- function(pa,
                          linewidth = 100,
                          match_dots = NULL,
                          print_pos = T,
                          print_pos_end = F,
                          use_align_starts = T,
                          col_out = T,
                          extend_subject = c(0,0),
                          out_file = NULL) {

  # pa formats
  #http://emboss.sourceforge.net/docs/themes/AlignFormats.html

  if (!is.null(out_file)) {
    if (!grepl("\\.txt", out_file)) {
      print("You may want to save the output to a .txt file. If so, have a file with .txt in the out_file path.")
    }
    sink(out_file)
  }

  if (methods::is(pa, "PairwiseAlignmentsSingleSubject")) {
    if (is.null(pa@subject@unaligned@ranges@NAMES) || is.null(pa@pattern@unaligned@ranges@NAMES)) {
      message("In order to pass names to a pairwiseAligment object, create named XStringSets to pass to pwalign::pairwiseAlignment")
    }

    p_name <- pa@pattern@unaligned@ranges@NAMES %||% "pattern"
    s_name <- pa@subject@unaligned@ranges@NAMES %||% "subject"

    ## extend subject?? - handle gaps - no not necessary, no gaps anyway
    #pattern <- as.character(pa@pattern) #Biostrings::pattern(pa)
    #subject <- as.character(pa@subject) #Biostrings::subject(pa)
    pattern <- as.character(pwalign::alignedPattern(pa))
    subject <- as.character(pwalign::alignedSubject(pa))

    c(subject, pattern) %<-% extend_subject_fun(extend_subject, pa, subject, pattern)

  } else if (methods::is(pa, "DNAStringSet") || methods::is(pa, "RNAStringSet") || methods::is(pa, "AAStringSet")) {
    if (length(pa) > 2) {
      message("pa has more than 2 sequences. Using first two.")
    } else if (length(pa) == 1) {
      stop("only 1 sequence in pa.")
    }
    pa <- as.character(pa)
    subject <- pa[2]
    pattern <- pa[1]
    p_name <- names(pattern) %||% "pattern"
    s_name <- names(subject) %||% "subject"
  } else {
    stop("pa must be pairwiseAlignment or XStringSet.")
  }


  if (nchar(subject) != nchar(pattern)) {
    stop("pattern and subject differ in length which cannot be handled.")
  }

  # make seq names equal in length
  c(p_name, s_name) %<-% pad_strings(c(p_name, s_name))
  linewidth <- min(c(linewidth, nchar(pattern)))

  if (use_align_starts && methods::is(pa, "PairwiseAlignmentsSingleSubject")) {
    p_start <- pa@pattern@range@start
    s_start <- pa@subject@range@start - extend_subject[1]
  } else {
    p_start <- 1
    s_start <- 1
  }

  ## split seqs to assign_dots_two_seq and get_symbols
  c(subject, pattern) %<-% strsplit(c(subject, pattern), "")

  c(subject, pattern) %<-% assign_dots_two_seq(
    match_dots = match_dots,
    subject = subject,
    pattern = pattern
  )
  # symbol lines above/below sequences
  p_sym <- paste(get_symbols(s = pattern, start.pos = p_start), collapse = "")
  s_sym <- paste(get_symbols(s = subject, start.pos = s_start), collapse = "")

  ## then re-collapse
  pattern <- paste(pattern, collapse = "")
  subject <- paste(subject, collapse = "")

  starts <- seq(1, nchar(pattern), by = linewidth)

  p_res <- s_res <- 0
  p_pos <- p_start - 1
  s_pos <- s_start - 1

  ## this is done to make the previous function work which expects start at 0
  while (p_pos %% 5 != 0) {
    p_pos <- p_pos - 1
  }
  while (s_pos %% 5 != 0) {
    s_pos <- s_pos - 1
  }

  # vectorized string split for every row
  c(s_chunk,
    p_chunk,
    p_sym_chunk,
    s_sym_chunk) %<-% stringr::str_sub_all(
      c(subject, pattern, p_sym, s_sym),
      starts,
      starts+linewidth-1
    )

  # end pos after every line, only needed when print_pos_end
  p_res <- p_res + cumsum(linewidth - stringr::str_count(p_chunk, "-"))
  s_res <- s_res + cumsum(linewidth - stringr::str_count(s_chunk, "-"))

  # add name spacing before sym_line on every line
  p_sym_chunk2 <- paste0(paste(rep(" ", nchar(p_name)), collapse = ""), p_sym_chunk)
  s_sym_chunk2 <- paste0(paste(rep(" ", nchar(s_name)), collapse = ""), s_sym_chunk)

  if (col_out) {
    p_chunk <- purrr::map_chr(p_chunk, col_letters)
    s_chunk <- purrr::map_chr(s_chunk, col_letters)
  }
  # add name before seq lines
  p_chunk <- paste0(p_name, p_chunk, "  ")
  s_chunk <- paste0(s_name, s_chunk, "  ")
  # add position at end
  if (print_pos_end) {
    p_chunk <- paste0(p_chunk, p_res + p_start - 1)
    s_chunk <- paste0(s_line, s_res + s_start - 1)
  }

  for (i in seq_along(p_sym_chunk)) {
    # did not know how to vectorize, so: loop
    c(p_sym_chunk[i], p_pos) %<-% .get_pos_line(name = p_name, chunk = p_sym_chunk[i], pos = p_pos)
    c(s_sym_chunk[i], s_pos) %<-% .get_pos_line(name = s_name, chunk = s_sym_chunk[i], pos = s_pos)
  }
  # browser()

  ## print
  for (i in seq_along(p_chunk)) {
    if (print_pos) {
      cat(p_sym_chunk[i], "\n")
      cat(p_sym_chunk2[i], "\n")
    }
    cat(p_chunk[i], "\n")
    cat(s_chunk[i], "\n")
    if (print_pos) {
      cat(s_sym_chunk2[i], "\n")
      cat(s_sym_chunk[i], "\n")
    }
    cat(paste(" "), "\n")
  }
  cat(paste(" "), "\n")
  cat(paste(" "), "\n")
  cat(paste(" "), "\n")

  # for (i in 1:length(starts)) {
  #   p_chunk <- substring(pattern, starts[i], starts[i]+linewidth-1)
  #   s_chunk <- substring(subject, starts[i], starts[i]+linewidth-1)
  #
  #   p_res <- p_res + linewidth - Biostrings::countPattern("-",p_chunk)
  #   s_res <- s_res + linewidth - Biostrings::countPattern("-",s_chunk)
  #
  #   if (col_out) {
  #     p_chunk <- .col.letter.fun(p_chunk)
  #     s_chunk <- .col.letter.fun(s_chunk)
  #   }
  #   p_sym_chunk <- substring(p_sym, starts[i], starts[i]+linewidth-1)
  #   s_sym_chunk <- substring(s_sym, starts[i], starts[i]+linewidth-1)
  #
  #   p_sym.line <- paste0(paste(rep(" ", nchar(p_name)), collapse = ""), p_sym_chunk)
  #   s_sym.line <- paste0(paste(rep(" ", nchar(s_name)), collapse = ""), s_sym_chunk)
  #
  #   if (print_pos_end) {
  #     p_line <- paste0(p_name, p_chunk, "  ", p_res + p_start - 1)
  #     s_line <- paste0(s_name, s_chunk, "  ", s_res + s_start - 1)
  #   } else {
  #     p_line <- paste0(p_name, p_chunk, "  ")
  #     s_line <- paste0(s_name, s_chunk, "  ")
  #   }
  #
  #   p_list <- .get_pos_line(name = p_name, chunk = p_sym_chunk, pos = p_pos)
  #   p_pos.line <- p_list[[1]]
  #   p_pos <- p_list[[2]]
  #
  #   s_list <- .get_pos_line(name = s_name, chunk = s_sym_chunk, pos = s_pos)
  #   s_pos.line <- s_list[[1]]
  #   s_pos <- s_list[[2]]
  #
  #   if (print_pos) {
  #     cat(p_pos.line, "\n")
  #     cat(p_sym.line, "\n")
  #   }
  #   cat(p_line, "\n")
  #   cat(s_line, "\n")
  #   if (print_pos) {
  #     cat(s_sym.line, "\n")
  #     cat(s_pos.line, "\n")
  #   }
  #   cat(paste(" "), "\n")
  #
  # }
  # cat(paste(" "), "\n")
  # cat(paste(" "), "\n")
  # cat(paste(" "), "\n")

  if (!is.null(out_file)) {
    closeAllConnections()
  }


}


.get_pos_line <- function(name, chunk, pos) {
  pos_line <- paste(rep(" ", nchar(name)), collapse = "")
  skip = 0
  corr.100.per.line <- F
  for (j in 1:nchar(chunk)) {
    if (stringr::str_sub(chunk, j, j) == "|") {
      pos <- pos + 5
      if (nchar(pos) == 3 & !corr.100.per.line) {
        corr.100.per.line <- T
        pos_line <- stringr::str_sub(pos_line, 1, -2)
      }
      if (pos %% 10 == 0) {
        pos_line <- paste0(pos_line, pos)
        skip <- skip + nchar(pos) - 1
      } else {
        if (skip == 0) {
          pos_line <- paste0(pos_line, " ")
        } else {
          skip <- skip - 1
        }
      }
    } else {
      if (skip == 0) {
        pos_line <- paste0(pos_line, " ")
      } else {
        skip <- skip - 1
      }
    }
  }
  return(list(pos_line, pos))
}

get_symbols <- function(s, start.pos = 1) {
  corr.pos <- 0
  s.sym <- ""
  for (j in start.pos:(length(s) + start.pos - 1)) {
    if (s[j- start.pos + 1] != "-" & (j-corr.pos) %% 5 != 0) {
      s.sym <- c(s.sym, ".")
    } else if (s[j- start.pos + 1] != "-" & (j-corr.pos) %% 5 == 0) {
      s.sym <- c(s.sym, "|")
    } else if (s[j- start.pos + 1] == "-") {
      s.sym <- c(s.sym, " ")
      corr.pos <- corr.pos + 1
    }
  }
  return(strsplit(paste(s.sym, collapse = ""), "")[[1]])
}


.col.letter.fun <- function(x) {


  ## add an option to color by specific aa properties or by groups
  '
  ## from Peptides::aaComp
aa_list <- list(tiny = c("A", "C", "G", "S", "T"),
                small = c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V"),
                aliphatic = c("A", "I", "L", "V"),
                aromatic = c("F", "H", "W", "Y"),
                nonpolar = c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"),
                polar = c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"),
                charged = c("B", "D", "E", "H", "K", "R", "Z"), ## all which are basic or acidic
                basic = c("H", "K", "R"),
                acidic = c("B", "D", "E", "Z"))

aa_df <- stack(aa_list)
aa_df <- aa_df[which(aa_df$values %in% Peptides::aaList()),] #Biostrings::AA_STANDARD
unstack(aa_df[,c(2,1)])

Peptides::aIndex(Peptides::aaList()) # aliphaticness
sort(setNames(Peptides::charge(Peptides::aaList()), Peptides::aaList())) # charge = dist to zero (neg or pos); polar = has charge and nonpolar = has no charge (roughly, not 100 %, e.g. S is polar without charge)
sort(setNames(Peptides::hydrophobicity(Peptides::aaList(), scale = "Eisenberg"), Peptides::aaList())) # polarity: all negative ones are polar, all positives nonpolar
setNames(Peptides::pI(Peptides::aaList()), Peptides::aaList()) # low pI = acidic, high pI = basic
sort(setNames(Peptides::mw(Peptides::aaList()), Peptides::aaList())) # roughly size
  '

  black <- crayon::make_style("black")
  whiter <- crayon::make_style("white")

  x <- strsplit(x, "")[[1]]

  if (all(x == "-")) {
    x <- sapply(x, function(y) {
      crayon::make_style("grey30", bg = T)(whiter(y))
    })
  } else if (igsc:::guess_type(x) == "NT") {

    dark_grey_bg_letters <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "-")
    cols <- RColorBrewer::brewer.pal(6, "Set2")[-c(4,5)]

    x <- sapply(x, function(y) {
      if (y == "A") {
        #crayon::make_style(acp[["A"]], bg = T)(black("A"))
        crayon::make_style(cols[1], bg = T)(black("A"))
      } else if (y == "T" || y == "U") {
        #crayon::make_style(acp[["T"]], bg = T)(black("T"))
        crayon::make_style(cols[2], bg = T)(black(y))
      } else if (y == "C") {
        #crayon::make_style(acp[["C"]], bg = T)(black("C"))
        crayon::make_style(cols[3], bg = T)(black("C"))
      } else if (y == "G") {
        #crayon::make_style(acp[["G"]], bg = T)(black("G"))
        crayon::make_style(cols[4], bg = T)(black("G"))
      } else if (y == "N") {
        #crayon::make_style("grey95", bg = T)(black("N"))
        crayon::make_style("grey40", bg = T)(whiter("N"))
      } else if (y %in% dark_grey_bg_letters) {
        crayon::make_style("grey30", bg = T)(whiter(y))
      } else {
        y
      }
    })

  } else {
    cols <- RColorBrewer::brewer.pal(7, "Set2")[-c(4:5)]

    aa_non_polar <- c("A","V","L","I","M","W","F","Y")
    aa_polar <- c("S","T","N","Q")
    aa_pos_charge <- c("K","R","H")
    aa_neg_charge <- c("D","E")
    aa_special <- c("C", "G", "P")

    x <- sapply(x, function(y) {
      if (y %in% aa_non_polar) {
        crayon::make_style(cols[1], bg = T)(black(y))
      } else if (y %in% aa_polar) {
        crayon::make_style(cols[2], bg = T)(black(y))
      } else if (y %in% aa_pos_charge) {
        crayon::make_style(cols[3], bg = T)(black(y))
      } else if (y %in% aa_neg_charge) {
        crayon::make_style(cols[4], bg = T)(black(y))
      } else if (y %in% aa_special) {
        crayon::make_style(cols[5], bg = T)(black(y))
      } else if (y == "N") {
        crayon::make_style("grey40", bg = T)(whiter("N"))
      } else {
        crayon::make_style("grey30", bg = T)(whiter(y))
      }
    })
  }

  x <- paste(x,collapse = "")

  return(x)
}

col_letters <- function(x, type = NULL, mc.cores = 1) {

  # mc.cores: probably no speedup

  blacker <- crayon::make_style("black")
  whiter <- crayon::make_style("white")

  x <- strsplit(x, "")[[1]]

  if (is.null(type)) {
    type <- igsc:::guess_type(x)
  } else {
    type <- match.arg(type, c("AA", "NT"))
  }

  if (type == "NT") {

    # nucleotides
    x <- style_fun_nt(x = x, blacker = blacker, whiter = whiter, mc.cores = mc.cores)

  } else if (type == "AA") {

    # amino acids
    x <- style_fun_aa(x = x, blacker = blacker, whiter = whiter, mc.cores = mc.cores)

  }

  # if (all(!grepl("[A-Za-z]", x))) {
  #
  #   # only dash or dot, no letters, unspecified
  #   x <- sapply(x, function(y) crayon::make_style("grey30", bg = T)(whiter(y)))
  #
  # }

  x <- paste(x,collapse = "")
  return(x)
}


style_fun_aa <- function(x,
                         blacker = crayon::make_style("black"),
                         whiter = crayon::make_style("white"),
                         # cols only assigned by order yet
                         cols = RColorBrewer::brewer.pal(7, "Set2")[-c(4,5)],
                         mc.cores = getOption("mc.cores", default = 1)) {

  # aa_non_polar <- c("A","V","L","I","M","W","F","Y")
  # aa_polar <- c("S","T","N","Q")
  # aa_pos_charge <- c("K","R","H")
  # aa_neg_charge <- c("D","E")
  # aa_special <- c("C", "G", "P")
  # x <- sapply(x, function(y) {
  #   if (y %in% aa_non_polar) {
  #     crayon::make_style(cols[1], bg = T)(blacker(y))
  #   } else if (y %in% aa_polar) {
  #     crayon::make_style(cols[2], bg = T)(blacker(y))
  #   } else if (y %in% aa_pos_charge) {
  #     crayon::make_style(cols[3], bg = T)(blacker(y))
  #   } else if (y %in% aa_neg_charge) {
  #     crayon::make_style(cols[4], bg = T)(blacker(y))
  #   } else if (y %in% aa_special) {
  #     crayon::make_style(cols[5], bg = T)(blacker(y))
  #   } else if (y == "N") {
  #     crayon::make_style("grey40", bg = T)(whiter(y))
  #   } else {
  #     crayon::make_style("grey30", bg = T)(whiter(y))
  #   }
  # })
  #

  dark_grey_bg_letters <- c("N")  # if you want, can add others later
  style_map <- list(
    # Non-polar
    A = list(cols[1], blacker),
    V = list(cols[1], blacker),
    L = list(cols[1], blacker),
    I = list(cols[1], blacker),
    M = list(cols[1], blacker),
    W = list(cols[1], blacker),
    F = list(cols[1], blacker),
    Y = list(cols[1], blacker),

    # Polar
    S = list(cols[2], blacker),
    T = list(cols[2], blacker),
    N = list(cols[2], blacker),
    Q = list(cols[2], blacker),

    # Positive charge
    K = list(cols[3], blacker),
    R = list(cols[3], blacker),
    H = list(cols[3], blacker),

    # Negative charge
    D = list(cols[4], blacker),
    E = list(cols[4], blacker),

    # Special
    C = list(cols[5], blacker),
    G = list(cols[5], blacker),
    P = list(cols[5], blacker),

    # Special case for ambiguous or unknown or stop codon
    N = list("grey40", whiter),
    "*" = list("grey40", whiter)
  )

  x <- parallel::mclapply(x, function(y) {
    if (y %in% names(style_map)) {
      col <- style_map[[y]][[1]]
      wrapf <- style_map[[y]][[2]]
      y <- crayon::make_style(col, bg = TRUE)(wrapf(y))
    } else if (y %in% dark_grey_bg_letters) {
      y <- crayon::make_style("grey30", bg = TRUE)(whiter(y))
    } else {
     # y <- crayon::make_style("grey30", bg = T)(whiter(y))
    }
    return(y)
  }, mc.cores = mc.cores)

  return(unlist(x))
}

style_fun_nt <- function(x,
                         blacker = crayon::make_style("black"),
                         whiter = crayon::make_style("white"),
                         # cols only assigned by order yet
                         cols = RColorBrewer::brewer.pal(6, "Set2")[-c(4,5)],
                         mc.cores = getOption("mc.cores", default = 1)) {

  # ambuities in DNA / nucleotide code
  # ------------------------------------------------------------
  # Code | Base(s) Represented | Meaning
  # ------------------------------------------------------------
  # M    | A or C              | aMino bases
  # R    | A or G              | puRine bases
  # W    | A or T              | Weak (2 H-bonds)
  # S    | C or G              | Strong (3 H-bonds)
  # Y    | C or T              | pYrimidine bases
  # K    | G or T              | Keto bases
  # V    | A or C or G         | not T (V = "not T")
  # H    | A or C or T         | not G (H = "not G")
  # D    | A or G or T         | not C (D = "not C")
  # B    | C or G or T         | not A (B = "not A")
  # ------------------------------------------------------------
  # A    | A                   | Adenine
  # C    | C                   | Cytosine
  # G    | G                   | Guanine
  # T    | T                   | Thymine
  # U    | U                   | Uracil (in RNA)
  # ------------------------------------------------------------
  dark_grey_bg_letters <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "-")
  style_map <- list(
    A = list(cols[1], blacker),
    T = list(cols[2], blacker),
    U = list(cols[2], blacker),
    C = list(cols[3], blacker),
    G = list(cols[4], blacker),
    N = list("grey40", whiter)
  )

  x <- parallel::mclapply(x, function(y) {
    if (y %in% names(style_map)) {
      col <- style_map[[y]][[1]]
      wrapf <- style_map[[y]][[2]]
      y <- crayon::make_style(col, bg = TRUE)(wrapf(y))
    } else if (y %in% dark_grey_bg_letters) {
      y <- crayon::make_style("grey30", bg = TRUE)(whiter(y))
    } else {
      #nothing
    }
    return(y)
  }, mc.cores = mc.cores)

  return(unlist(x))
}




pad_strings <- function(x, add_default = "  ") {
  lens <- nchar(x)
  max_len <- max(lens)
  x[lens == max_len] <- paste0(x[lens == max_len], add_default)
  new_max <- max(nchar(x))
  x <- sprintf(paste0("%-", new_max, "s"), x)

  return(x)
}

extend_subject_fun <- function(extend_subject, pa, subject, pattern) {
  if (!all(extend_subject == 0)) {
    # work on it
    if (!is.numeric(extend_subject)) {
      stop("extend_subject has to be numeric.")
    }
    if (any(extend_subject < 0)) {
      stop("both extend_subject > 0!")
    }
    if (length(extend_subject) != 2) {
      stop("extend_subject has to be of length 2.")
    }
    # avoid diff length of pattern and subject
    extend_subject[1] <- min(c(pa@subject@range@start-1, extend_subject[1]))
    extend_subject[2] <- min(c(nchar(as.character(pa@subject@unaligned[[1]]))-c(pa@subject@range@start+pa@subject@range@width-1), extend_subject[2]))

    subject <- paste0(substr(as.character(pa@subject@unaligned[[1]]), start = pa@subject@range@start - extend_subject[1], stop = pa@subject@range@start-1), subject)
    subject <- paste0(subject, substr(as.character(pa@subject@unaligned[[1]]), start = pa@subject@range@start+pa@subject@range@width, stop = pa@subject@range@start+pa@subject@range@width-1+extend_subject[2]))

    pattern <- paste0(paste(rep("-", extend_subject[1]), collapse = ""), pattern)
    pattern <- paste0(pattern, paste(rep("-", extend_subject[2]), collapse = ""))
  }

  return(list(subject, pattern))
}

assign_dots_two_seq <- function(match_dots, match_symbol = ".",
                                subject, pattern) {
  if (!is.null(match_dots)) {
    match_dots <- match.arg(match_dots, c("subject", "pattern"), several.ok = T)
    #c(subject, pattern) %<-% strsplit(c(subject, pattern), "")

    same <- pattern == subject
    if ("subject" %in% match_dots) {
      subject[which(same)] <- match_symbol
    }
    if ("pattern" %in% match_dots) {
      pattern[which(same)] <- match_symbol
    }
    # pattern <- paste(pattern, collapse = "")
    # subject <- paste(subject, collapse = "")
  }
  return(list(subject, pattern))
}
