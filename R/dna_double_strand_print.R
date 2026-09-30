#' Print a DNA Double Strand
#'
#' Print two DNA strands in aligned blocks to the console or a text file.
#'
#' @param top,bottom A single DNA sequence, supplied as a character string,
#'   DNAString, or DNAStringSet of length one. Supply sequences in display
#'   order: top from 5' to 3', bottom from 3' to 5'. At least one is required.
#'   IUPAC DNA letters, spaces, dashes, and dots are accepted. Spaces and gap
#'   characters retain their positions, allowing overhangs on either side.
#'   A name on a character vector or DNAStringSet is used as the row label.
#' @param linewidth Positive integer; number of alignment columns per block.
#' @param print_ends Logical; print 5' and 3' at the ends of both rows in each
#'   block. These labels indicate direction, including in wrapped blocks.
#' @param print_pipes Logical; print pipes between complementary A/T and C/G
#'   pairs. Mismatches, ambiguous bases, gaps, and overhangs receive spaces.
#' @param col_out Logical; color nucleotides using the same helper as
#'   [pwalign_print()] and [xstringset_print()].
#' @param out_file Path to a text file, or NULL to print to the console.
#' @param align Align strands to the left (default) or right of their shared
#'   display width. Padding is added before wrapping into blocks; existing
#'   spaces and gaps are preserved.
#' @param first_n,last_n Optional positive integer; print only the first or
#'   last n alignment columns, respectively. Supply at most one. Columns are
#'   counted from left to right after alignment, including gaps and overhangs.
#'   Values larger than the display width show the full strands. NULL (the
#'   default) leaves output unrestricted. The returned strands remain complete.
#'
#' @details
#' When one strand is missing, [revcomp_dna()] generates its complement using
#'   the Biostrings implementation, which supports IUPAC ambiguity codes.
#'   Reversal is disabled because the strands are displayed antiparallel.
#'   The supplied strand retains its case; the generated strand is uppercase.
#' When both strands are supplied, their sequences are printed as provided,
#'   without sequence alignment or complementation. Unequal lengths are padded
#'   on the right when `align = "left"`, or on the left when `align = "right"`.
#'   Use spaces or gaps to specify additional overhangs.
#' Row labels and nucleotide colors reuse the helpers shared by
#'   [pwalign_print()] and [xstringset_print()].
#'
#' @return Invisibly, a named character vector with top and bottom strands in
#'   display order, including any generated strand, without added padding.
#' @export
#'
#' @examples
#' dna_double_strand_print(top = "ATGCAA", col_out = FALSE)
#' dna_double_strand_print(bottom = "TACGTT", col_out = FALSE)
#' dna_double_strand_print(top = "  ATGCAA", bottom = "GGTACG",
#'                        col_out = FALSE)
#' dna_double_strand_print(top = "ATGCAA", bottom = "CGTT",
#'                        align = "right", col_out = FALSE)
#' dna_double_strand_print(top = "ATGCAA", first_n = 3, col_out = FALSE)
#' dna_double_strand_print(top = "ATGCAA", last_n = 3, col_out = FALSE)
#' dna_double_strand_print(top = "ATGCAA", print_ends = FALSE,
#'                        print_pipes = FALSE, col_out = FALSE)
dna_double_strand_print <- function(top = NULL,
                                    bottom = NULL,
                                    linewidth = 100,
                                    print_ends = TRUE,
                                    print_pipes = TRUE,
                                    col_out = TRUE,
                                    out_file = NULL,
                                    align = c("left", "right"),
                                    first_n = NULL,
                                    last_n = NULL) {

  align <- match.arg(align)
  if (!is.null(first_n) && !is.null(last_n)) {
    stop("Supply at most one of first_n or last_n.", call. = FALSE)
  }
  limits <- list(first_n = first_n, last_n = last_n)
  for (limit in names(limits)) {
    value <- limits[[limit]]
    if (!is.null(value) &&
        (!is.numeric(value) || length(value) != 1L || is.na(value) ||
         !is.finite(value) || value < 1 || value != floor(value))) {
      stop(limit, " must be NULL or a positive integer.", call. = FALSE)
    }
  }
  if (is.null(top) && is.null(bottom)) {
    stop("Supply at least one of top or bottom.", call. = FALSE)
  }
  if (!is.numeric(linewidth) || length(linewidth) != 1L ||
      is.na(linewidth) || !is.finite(linewidth) || linewidth < 1 ||
      linewidth != floor(linewidth)) {
    stop("linewidth must be a positive integer.", call. = FALSE)
  }
  flags <- list(print_ends = print_ends, print_pipes = print_pipes,
                col_out = col_out)
  for (flag in names(flags)) {
    value <- flags[[flag]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(flag, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.null(out_file) &&
      (!is.character(out_file) || length(out_file) != 1L ||
       is.na(out_file) || !nzchar(out_file))) {
    stop("out_file must be NULL or a single file path.", call. = FALSE)
  }

  row_names <- c("top", "bottom")
  strands <- list(top, bottom)
  for (i in seq_along(strands)) {
    strand <- strands[[i]]
    if (is.null(strand)) next
    if (!is.character(strand) && !methods::is(strand, "DNAString") &&
        !methods::is(strand, "DNAStringSet")) {
      stop(row_names[i], " must be a DNA string.", call. = FALSE)
    }
    strand_name <- names(strand)
    strand <- as.character(strand)
    if (length(strand) != 1L || is.na(strand) || !nzchar(strand) ||
        grepl("[^ACGTRYSWKMBDHVN .-]", toupper(strand)) ||
        !grepl("[ACGTRYSWKMBDHVN]", toupper(strand))) {
      stop(row_names[i], " must contain one nonempty DNA sequence; ",
           "only IUPAC DNA letters, spaces, dashes, and dots are allowed.",
           call. = FALSE)
    }
    if (length(strand_name) == 1L && !is.na(strand_name) &&
        nzchar(strand_name)) {
      row_names[i] <- strand_name
    }
    strands[[i]] <- unname(strand)
  }

  missing_strand <- vapply(strands, is.null, logical(1))
  if (any(missing_strand)) {
    supplied <- strsplit(strands[[which(!missing_strand)]], "")[[1]]
    bases <- grepl("[A-Za-z]", supplied)
    # Complement only the bases, retaining display gaps and overhang spacing.
    supplied[bases] <- strsplit(revcomp_dna(
      paste0(supplied[bases], collapse = ""), fun = "Biostrings", rev = FALSE
    ), "")[[1]]
    strands[[which(missing_strand)]] <- paste0(supplied, collapse = "")
  }
  result <- stats::setNames(unlist(strands, use.names = FALSE), c("top", "bottom"))

  # Pad names with the shared printer helper, and sequences to a common width.
  labels <- igsc:::pad_strings(row_names)
  width <- max(nchar(result))
  padding <- strrep(" ", width - nchar(result))
  padded <- if (align == "left") paste0(result, padding) else paste0(padding, result)
  if (!is.null(first_n)) {
    width <- min(first_n, width)
    padded <- substr(padded, 1, width)
  } else if (!is.null(last_n)) {
    first_column <- max(1, width - last_n + 1)
    padded <- substr(padded, first_column, width)
    width <- width - first_column + 1
  }
  starts <- seq.int(1, width, by = min(linewidth, width))
  chunks <- stringr::str_sub_all(padded, starts, pmin(starts + linewidth - 1, width))

  if (col_out) igsc:::.ensure_packages(c("crayon", "RColorBrewer"))
  destination <- stdout()
  if (!is.null(out_file)) {
    destination <- file(out_file, open = "wt")
    on.exit(close(destination), add = TRUE)
  }

  for (i in seq_along(starts)) {
    block <- vapply(chunks, `[`, character(1), i)
    if (print_pipes) {
      bases <- strsplit(toupper(block), "")
      paired <- bases[[1]] %in% c("A", "C", "G", "T") &
        chartr("ACGT", "TGCA", bases[[1]]) == bases[[2]]
      pipes <- paste0(ifelse(paired, "|", " "), collapse = "")
    }
    if (col_out) block <- vapply(block, igsc:::col_letters, character(1), type = "NT")
    left <- if (print_ends) c("5' ", "3' ") else c("", "")
    right <- if (print_ends) c(" 3'", " 5'") else c("", "")
    lines <- paste0(labels, left, block, right)
    cat(lines[1], "\n", sep = "", file = destination)
    if (print_pipes) {
      cat(strrep(" ", nchar(labels[1]) + nchar(left[1])), pipes, "\n",
          sep = "", file = destination)
    }
    cat(lines[2], "\n", sep = "", file = destination)
    if (i < length(starts)) cat("\n", file = destination)
  }

  invisible(result)
}
