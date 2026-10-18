#' Write reads to a FASTQ file
#'
#' Export reads supplied as a data frame, individual vectors, or both.
#'
#' @param reads Optional data frame containing `qname`, `seq`, `qual`,
#'   and `strand` columns.
#' @param filepath Output path. A `.gz` suffix enables gzip compression.
#'   Parent directories must exist. Existing files are overwritten.
#' @param qname Optional read names. Overrides `reads$qname`. If unavailable,
#'   names `read_1`, `read_2`, etc. are generated and a message is emitted.
#' @param seq Character vector of DNA sequences. Overrides `reads$seq`.
#' @param qual Character vector of ASCII Phred+33 quality strings.
#'   Overrides `reads$qual`. Each string must match its sequence length.
#' @param strand Optional strand vector containing `"+"`, `"-"`, or `"*"`.
#'   Overrides `reads$strand`. A single value is recycled across reads.
#'   If unavailable, all strands default to `"+"`, with a message.
#'   Individual `NA` values are also replaced with `"+"`, with a message.
#' @param revcomp Logical; reverse-complement minus-strand sequences and
#'   reverse their quality strings. Use `TRUE` for default [get_bam_reads()]
#'   output. Use `FALSE` if sequences already have the desired orientation,
#'   including output extracted with `revcomp_minus_strand = TRUE` and
#'   `revcomp_plus_strand = FALSE`.
#' @param plus_line how to write third line per read? plus only (empty) or
#'   repeat qname?
#'
#' @details
#' The function cannot infer orientation when strand information is absent.
#' Defaulting to `"+"` leaves those sequences and qualities unchanged,
#' even when `revcomp = TRUE`. Explicit `"*"` strands are also unchanged.
#'
#' Each input row or vector element produces one FASTQ record. Records are
#' neither filtered nor deduplicated. Paired reads are written in input order
#' to the same file; they are not separated or arranged into pairs.
#'
#' Missing sequences and qualities are rejected. Alignment information,
#' barcodes, UMIs, and other columns are not exported. Hard-clipped bases
#' cannot be recovered from the supplied sequences.
#'
#' @return The output path.
#' @seealso [get_bam_reads()], [Biostrings::writeXStringSet()]
#'
#' @examples
#' # Example reads in BAM orientation
#' reads <- data.frame(
#'   qname = c("read1", "read2"),
#'   seq = c("AACC", "GGTT"),
#'   qual = c("IJKL", "ABCD"),
#'   strand = c("+", "-")
#' )
#'
#' # Restore sequenced orientation and write compressed FASTQ
#' outfile <- tempfile(fileext = ".fastq.gz")
#' write_reads_fastq(reads, outfile, revcomp = TRUE)
#'
#' exported <- Biostrings::readDNAStringSet(
#'   outfile,
#'   format = "fastq",
#'   with.qualities = TRUE
#' )
#' as.character(exported)
#' # read1: AACC; read2: AACC
#'
#' unlink(outfile)
#'
#' # Write sequences and qualities without changing their orientation
#' outfile <- tempfile(fileext = ".fastq")
#' write_reads_fastq(reads, outfile, revcomp = FALSE)
#' unlink(outfile)
#'
#' \dontrun{
#' # Extract directly in sequenced orientation
#' reads <- get_bam_reads(
#'   bam = "sample.bam",
#'   granges = regions,
#'   scores = FALSE,
#'   revcomp_minus_strand = TRUE,
#'   revcomp_plus_strand = FALSE
#' )
#'
#' # No additional reverse-complementing is needed
#' write_reads_fastq(reads, "sample.fastq.gz", revcomp = FALSE)
#' }
#' @export
write_fastq <- function(reads = NULL,
                        filepath = file.path(getwd(), "reads.fastq.gz"),
                        qname = NULL,
                        seq = NULL,
                        qual = NULL,
                        strand = NULL,
                        revcomp = TRUE,
                        plus_line = c("empty", "qname")) {

  if (!is.null(reads) && !is.data.frame(reads)) {
    stop("`reads` must be a data frame or NULL.", call. = FALSE)
  }

  if (!is.logical(revcomp) || length(revcomp) != 1L || is.na(revcomp)) {
    stop("`revcomp` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(filepath) || length(filepath) != 1L ||
      is.na(filepath) || !nzchar(filepath)) {
    stop("`filepath` must be a non-empty file path.", call. = FALSE)
  }

  # Explicit vectors take precedence over data-frame columns.
  resolve <- function(value, column) {
    if (!is.null(value)) value else reads[[column]]
  }

  plus_line <- rlang::arg_match(plus_line)

  qname <- resolve(qname, "qname")
  seq <- resolve(seq, "seq")
  qual <- resolve(qual, "qual")
  strand <- resolve(strand, "strand")

  if (is.null(seq) || is.null(qual)) {
    stop("Supply both `seq` and `qual`, individually or in `reads`.",
         call. = FALSE)
  }

  if (!is.character(seq) || !is.character(qual)) {
    stop("`seq` and `qual` must be character vectors.", call. = FALSE)
  }

  n <- length(seq)

  if (length(qual) != n) {
    stop("`seq` and `qual` must have equal lengths.", call. = FALSE)
  }

  if (anyNA(seq) || any(!nzchar(seq)) || any(seq == "*")) {
    stop("Read sequences must be present and non-empty.", call. = FALSE)
  }

  if (anyNA(qual) || any(!grepl("^[!-~]+$", qual)) ||
      any(nchar(qual, type = "bytes") != nchar(seq, type = "bytes"))) {
    stop(
      "Missing or invalid qualities, or sequence/quality length mismatch.",
      call. = FALSE
    )
  }

  if (is.null(qname)) {
    message(
      "`qname` is missing; generating names read_1, read_2, ... for ",
      n, " reads."
    )
    qname <- sprintf("read_%d", seq_len(n))
  }

  qname <- as.character(qname)

  if (length(qname) != n || anyNA(qname) || any(!nzchar(qname)) ||
      any(grepl("[[:space:]]", qname))) {
    stop("Supply one non-empty, whitespace-free name per read.",
         call. = FALSE)
  }

  if (is.null(strand)) {
    message("`strand` is missing; assuming '+' for all reads.")
    strand <- rep("+", n)
  } else {
    strand <- as.character(strand)

    if (length(strand) == 1L) {
      strand <- rep(strand, n)
    }

    if (length(strand) != n) {
      stop("`strand` must have length one or match the number of reads.",
           call. = FALSE)
    }

    if (anyNA(strand)) {
      message("Missing `strand` values detected; replacing them with '+'.")
      strand[is.na(strand)] <- "+"
    }

    if (any(!strand %in% c("+", "-", "*"))) {
      stop("`strand` values must be '+', '-', '*', or NA.", call. = FALSE)
    }
  }

  sequences <- Biostrings::DNAStringSet(seq)

  if (revcomp) {
    minus <- which(strand == "-")
    sequences[minus] <- Biostrings::reverseComplement(sequences[minus])
    qual[minus] <- stringi::stri_reverse(qual[minus])
  }

  names(sequences) <- qname

  # Biostrings::writeXStringSet(
  #   sequences,
  #   filepath = filepath,
  #   format = "fastq",
  #   qualities = Biostrings::BStringSet(qual),
  #   compress = grepl("\\.gz$", filepath, ignore.case = TRUE)
  # )

  con <- if (grepl("\\.gz$", filepath, ignore.case = TRUE)) {
    gzfile(filepath, open = "wb", compression = 1)
  } else {
    file(filepath, open = "wb")
  }
  on.exit(close(con), add = TRUE)

  sequences <- as.character(sequences)
  separator <- if (plus_line == "qname") {
    paste0("+", qname)
  } else {
    rep("+", length(qname))
  }
  qname <- paste0("@", qname)
  lines <- as.vector(rbind(
    qname,
    sequences,
    separator,
    qual
  ))

  vroom::vroom_write_lines(lines, file = con)
  return(filepath)
}
