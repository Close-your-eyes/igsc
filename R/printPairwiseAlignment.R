#' Print pairwise alignments of DNA or AA to console or txt file
#'
#' @param alignments one alignment or a list of multiple alignments created with
#' Biostrings::pariwiseAlignment
#' @param linewidth numeric, how many letters of alignment to print per line
#' @param match_dots should matching position be printed as dots for easier
#' recognition of differences; 'subject' to have dots in the subject, 'pattern' to
#' have them in the pattern
#' @param print_pos print the positions of pattern and subject above and below
#' respectively
#' @param print_pos_end print positions on the end of each line
#' @param use_align_starts if TRUE, the first positions are not 1 but refer to
#' position of alignment within the provided sequences
#' @param out_file path to a file where to print results to; if NULL results are
#' printed in console
#' @param col_out color the printed alignment?
#'
#' @return alignment in printed format in console or file
#' @export
#'
#' @examples
printPairwiseAlignment <- function(alignments,
                                   linewidth = 50,
                                   match_dots = NULL,
                                   print_pos = T,
                                   print_pos_end = F,
                                   use_align_starts = T,
                                   col_out = T,
                                   out_file = NULL) {

  if (!is.null(out_file)) {
    if (!grepl("\\.txt", out_file)) {
      print("You may want to save the output to a .txt file. If so, have a file with .txt in the out_file path.")
    }
    sink(out_file)
  }

  lapply(alignments, function(alignment) {

    if (is.null(alignment@subject@unaligned@ranges@NAMES) || is.null(alignment@pattern@unaligned@ranges@NAMES)) {
      print("In order to pass names to a pairwiseAligment object, create XStringSets (not XStrings; X for DNA, RNA or AA) of the sequences.")
    }

    if (is.null(alignment@pattern@unaligned@ranges@NAMES)) {
      p_name <- "pattern"
    } else {
      p_name <- alignment@pattern@unaligned@ranges@NAMES
    }

    if (is.null(alignment@subject@unaligned@ranges@NAMES)) {
      s_name <- "subject"
    } else {
      s_name <- alignment@subject@unaligned@ranges@NAMES
    }

    # make seq names equal in length
    seq.names <- c(p_name, s_name)
    seq.names[which.max(nchar(seq.names))] <- paste0(seq.names[which.max(nchar(seq.names))], "  ")
    seq.names[which.min(nchar(seq.names))] <- paste0(seq.names[which.min(nchar(seq.names))], paste(rep(" ", max(nchar(seq.names)) - min(nchar(seq.names))), collapse = ""))
    p_name <- seq.names[1]
    s_name <- seq.names[2]

    ## extend subject?? - handle gaps!
    pattern <- as.character(Biostrings::pattern(alignment))
    subject <- as.character(Biostrings::subject(alignment))

    if (nchar(subject) != nchar(pattern)) {
      stop("pattern and subject differ in length which cannot be handled.")
    }
    linewidth <- min(c(linewidth, nchar(pattern)))

    if (use_align_starts) {
      p_start <- alignment@pattern@range@start
      s_start <- alignment@subject@range@start
    } else {
      p_start <- 1
      s_start <- 1
    }

    if (!is.null(match_dots)) {
      match_dots <- match.arg(match_dots, c("subject", "pattern"))
      pattern <- strsplit(pattern, "")[[1]]
      subject <- strsplit(subject, "")[[1]]
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
      pattern <- paste(pattern, collapse = "")
      subject <- paste(subject, collapse = "")
    }

    pattern <- strsplit(pattern, "")[[1]]
    subject <- strsplit(subject, "")[[1]]

    p_sym <- paste(.get_symbols(s = pattern, start.pos = p_start), collapse = "")
    s_sym <- paste(.get_symbols(s = subject, start.pos = s_start), collapse = "")

    pattern <- paste(pattern, collapse = "")
    subject <- paste(subject, collapse = "")


    len  <- nchar(pattern)
    starts <- seq(1, len, by = linewidth)
    n <- length(starts)

    p_res <- 0
    s_res <- 0
    p_pos <- p_start - 1
    s_pos <- s_start - 1

    ## this is done to make the previous function work which expects start at 0
    while (p_pos %% 5 != 0) {
      p_pos <- p_pos - 1
    }
    while (s_pos %% 5 != 0) {
      s_pos <- s_pos - 1
    }

    for (i in 1:n) {
      p_chunk <- substring(pattern, starts[i], starts[i]+linewidth-1)
      s_chunk <- substring(subject, starts[i], starts[i]+linewidth-1)
      if (col_out) {
        p_chunk <- .col.letter.fun(p_chunk)
        s_chunk <- .col.letter.fun(s_chunk)
      }
      p_sym.chunk <- substring(p_sym, starts[i], starts[i]+linewidth-1)
      s_sym.chunk <- substring(s_sym, starts[i], starts[i]+linewidth-1)

      p_res <- p_res + linewidth - Biostrings::countPattern("-",p_chunk)
      s_res <- s_res + linewidth - Biostrings::countPattern("-",s_chunk)

      p_sym.line <- paste0(paste(rep(" ", nchar(p_name)), collapse = ""), p_sym.chunk)
      s_sym.line <- paste0(paste(rep(" ", nchar(s_name)), collapse = ""), s_sym.chunk)

      if (print_pos_end) {
        p_line <- paste0(p_name, p_chunk, "  ", p_res + p_start - 1)
        s_line <- paste0(s_name, s_chunk, "  ", s_res + s_start - 1)
      } else {
        p_line <- paste0(p_name, p_chunk, "  ")
        s_line <- paste0(s_name, s_chunk, "  ")
      }

      p_list <- .get_pos_line(name = p_name, chunk = p_sym.chunk, pos = p_pos)
      p_pos.line <- p_list[[1]]
      p_pos <- p_list[[2]]

      s_list <- .get_pos_line(name = s_name, chunk = s_sym.chunk, pos = s_pos)
      s_pos.line <- s_list[[1]]
      s_pos <- s_list[[2]]

      if (print_pos) {
        cat(p_pos.line, "\n")
        cat(p_sym.line, "\n")
      }
      cat(p_line, "\n")
      cat(s_line, "\n")
      if (print_pos) {
        cat(s_sym.line, "\n")
        cat(s_pos.line, "\n")
      }
      cat(paste(" "), "\n")

    }
  })

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

.get_symbols <- function(s, start.pos = 1) {
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

  black <- crayon::make_style("black")
  whiter <- crayon::make_style(rgb(1, 1, 1))

  x <- strsplit(x, "")[[1]]
  if (all(x %in% unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N")))) {

    dark_grey_bg_letters <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B")
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

  } else if (all(x %in% c(Biostrings::AA_ALPHABET, "N"))) {

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
  } else {
    print("AA, DNA or RNA could not be identified.")
  }

  x <- paste(x,collapse = "")

  return(x)
}