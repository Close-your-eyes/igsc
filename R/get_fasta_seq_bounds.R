#' Locate sequence boundaries in a FASTA file
#'
#' Return the inclusive file-line boundaries of each FASTA record, together
#' with its header, chromosome-style names, and an approximate sequence length.
#' For gzip input, line numbers refer to the decompressed text.
#'
#' @inheritParams get_fasta_names
#'
#' @section Supported formats:
#' Input must be a local plain-text FASTA file (for example, `.fa`, `.fasta`,
#' `.fna`, or `.faa`) or a gzip-compressed FASTA file ending in lowercase `.gz`
#' (for example, `.fa.gz` or `.fasta.gz`). Plain-text files do not need a specific
#' extension. Other compression formats, archives such as `.tar.gz`, and `fst`
#' files are not supported. FASTQ input is not supported. The sequence alphabet
#' is not validated.
#'
#' @details
#' Headers must begin with `">"` in the first column. Sequence identifiers must
#' be unique because [get_fasta_names()] is called with `prep_chr = TRUE`.
#' Records remain in file order; the name columns have chromosome-aware factor
#' levels as described in [get_fasta_names()].
#'
#' The first line is assumed to be a header and the second line a sequence
#' line. The trimmed width of that sequence line is used to estimate every
#' record's length. This estimate assumes a consistent sequence-line width
#' across the file. Short final lines cause overestimates; differing wrap
#' widths, blank lines, and comments can make the estimate inaccurate.
#' `n_seq_lines` counts every physical line after the header through `end_line`,
#' including blank and comment lines. Boundaries do not depend on line width.
#'
#' Line counting uses \code{\link{count_file_lines}}. Line counting and header detection
#' use the external `rg` command when
#' available and otherwise use `grep`. Gzip input uses `rg -z`, which requires
#' a gzip decompressor on `PATH`; the `grep` fallback uses `gunzip -c`.
#' Sequence-line width is read with [vroom::vroom_lines()].
#'
#' @return A data frame with one row per FASTA record, in file order, containing:
#' \describe{
#'   \item{fastaname}{Factor containing the complete header without its leading
#'     `">"`.}
#'   \item{seqname}{Factor containing the first space-delimited header field.}
#'   \item{fctname}{Factor containing `seqname` with `"chr"` removed, used to
#'     define chromosome-aware factor levels.}
#'   \item{start_line}{Numeric one-based line number of the record's header.}
#'   \item{end_line}{Numeric inclusive last line before the next header, or the
#'     final file line for the last record.}
#'   \item{n_seq_lines}{Numeric count of physical lines in the record,
#'     excluding its header.}
#'   \item{bp_approx}{Numeric approximate sequence length, calculated by
#'     multiplying `n_seq_lines` by the first sequence line's trimmed width.
#'     This is a base-count estimate for nucleotide FASTA input, not an exact
#'     sequence length.}
#' }
#' @seealso \code{\link{count_file_lines}}, [get_fasta_names()], [read_fasta()]
#' @export
#'
#' @examples
#' \dontrun{
#' fasta <- tempfile(fileext = ".fa")
#' records <- c(">chr1 chromosome 1", "ACGT", "TGCA", ">chr2", "AAAA")
#' writeLines(records, fasta)
#' get_fasta_seq_bounds(fasta)
#'
#' compressed_fasta <- paste0(fasta, ".gz")
#' con <- gzfile(compressed_fasta, "wt")
#' writeLines(records, con)
#' close(con)
#' get_fasta_seq_bounds(compressed_fasta)
#'
#' unlink(c(fasta, compressed_fasta))
#' }
get_fasta_seq_bounds <- function(file_path) {

  total_lines <- count_file_lines(file_path)

  # get width from second line, assuming this is the first line with sequence
  line_width <- nchar(trimws(vroom::vroom_lines(file = file_path,
                                                skip_empty_rows = T,
                                                skip = 1,
                                                n_max = 1,
                                                progress = F)))

  name_lines <- get_fasta_names(file_path, prep_chr = T)
  name_lines$end_line <- c(name_lines$start_line[-1]-1, total_lines)
  name_lines$n_seq_lines <- name_lines$end_line - name_lines$start_line # exclude name line
  name_lines$bp_approx <- line_width*name_lines$n_seq_lines
  if (all(c("seqname", "fctname") %in% names(name_lines))) {
    name_lines <- name_lines[,c(2:4,1,5:7)]
  }

  return(name_lines)
}

#' Count lines in a text file
#'
#' Count physical lines in a local plain-text or gzip-compressed file without
#' loading the entire file into R. Blank lines are included. A non-empty final
#' line is counted even when it has no terminating newline; an empty file has
#' zero lines. For gzip input, the count refers to decompressed text.
#'
#' @param file_path A single local path to a text-like file. Paths containing
#'   spaces or shell characters are supported, and `~` is expanded to the
#'   user's home directory.
#'
#' @section Supported formats:
#' Any line-oriented plain-text format, including TXT, CSV, TSV, JSON, GTF,
#' FASTA, and FASTQ, can be counted. No specific plain-text extension is
#' required. Gzip-compressed text must end in lowercase `.gz`. Other
#' compression formats, archives such as `.tar.gz`, and binary files such as
#' `fst` are not supported. Contents are not parsed or validated: the result
#' counts physical lines, not records, sequences, or table rows.
#'
#' @details
#' Uses the external `rg` command when available and otherwise uses `grep`.
#' Gzip input uses `rg -z`, which additionally requires `gzip` on `PATH`;
#' the `grep` fallback uses `gunzip -c`. File paths are quoted for shell use.
#'
#' @return A single non-negative numeric line count.
#' @seealso [get_fasta_seq_bounds()]
#' @export
#'
#' @examples
#' \dontrun{
#' path <- tempfile(fileext = ".txt")
#' writeLines(c("first line", "", "third line"), path)
#' count_file_lines(path) # 3, including the blank line
#'
#' compressed_path <- paste0(path, ".gz")
#' con <- gzfile(compressed_path, "wt")
#' writeLines(c("first line", "", "third line"), con)
#' close(con)
#' count_file_lines(compressed_path) # 3 decompressed lines
#'
#' unlink(c(path, compressed_path))
#' }
count_file_lines <- function(file_path) {
  if (!is.character(file_path) || length(file_path) != 1L ||
      is.na(file_path) || !nzchar(file_path)) {
    stop("file_path must be a single non-empty file path.", call. = FALSE)
  }
  file_path <- path.expand(file_path)
  if (!file.exists(file_path) || dir.exists(file_path)) {
    stop("file_path must refer to an existing file.", call. = FALSE)
  }

  count_command <- if (nzchar(Sys.which("rg"))) "rg" else "grep"
  if (!nzchar(Sys.which(count_command))) {
    stop("Line counting requires rg or grep on PATH.", call. = FALSE)
  }
  compressed <- tools::file_ext(file_path) == "gz"
  if (compressed) {
    decompressor <- if (count_command == "rg") "gzip" else "gunzip"
    if (!nzchar(Sys.which(decompressor))) {
      stop("Gzip input requires ", decompressor, " on PATH.", call. = FALSE)
    }
  }

  quoted_path <- shQuote(file_path)
  if (compressed && count_command == "rg") {
    cmd <- paste0("rg -z --count '^' -- ", quoted_path)
  } else if (compressed) {
    cmd <- paste0("gunzip -c -- ", quoted_path, " | grep --count '^'")
  } else {
    cmd <- paste0(count_command, " --count '^' -- ", quoted_path)
  }
  output <- suppressWarnings(system(cmd, intern = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (!status %in% c(0L, 1L)) {
    stop("Line counting failed for: ", file_path, call. = FALSE)
  }
  if (status == 1L && length(output) == 0L) return(0)

  count <- suppressWarnings(as.numeric(trimws(output)))
  if (length(count) != 1L || !is.finite(count) || count < 0) {
    stop("Could not determine the line count for: ", file_path, call. = FALSE)
  }
  return(count)
}
