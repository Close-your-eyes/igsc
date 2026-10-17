#' Extract sequence names from a FASTA file
#'
#' Find every FASTA header in a plain-text or gzip-compressed file and return
#' its line number and header text. The leading `">"` is removed from each
#' header.
#'
#' @section Supported formats:
#' Input must be a local plain-text FASTA file (for example, `.fa`, `.fasta`,
#' `.fna`, or `.faa`) or a gzip-compressed FASTA file ending in lowercase `.gz`
#' (for example, `.fa.gz` or `.fasta.gz`). Plain-text files do not need a specific
#' extension. Other compression formats, archives such as `.tar.gz`, and `fst`
#' files are not supported. FASTQ input is not supported. Headers must begin
#' with `">"` in the first column; sequence content is not validated.
#'
#' @details
#' When `prep_chr = TRUE`, the first space-delimited field of each header
#' is treated as a sequence name. The `"chr"` substring is removed to construct
#' an ordering key: numeric names are ordered numerically, followed by
#' non-numeric names in their input order. This is useful for chromosome FASTA
#' files, where lexical ordering would otherwise place `chr10` before `chr2`.
#' This ordering applies to factor levels; rows remain in file order. Sequence
#' identifiers and their ordering keys must be unique when `prep_chr = TRUE`.
#'
#' Header detection uses the external `rg` command when available and falls
#' back to `grep`. Gzip input uses `rg -z`, which requires a gzip decompressor
#' on `PATH`; the `grep` fallback uses `gunzip -c`. Line numbers refer to the
#' decompressed text for gzip input and to the original text for plain input.
#'
#' @param file_path A single local path to a plain-text or gzip-compressed FASTA
#'   file. Gzip input must end in lowercase `.gz`; all other extensions are
#'   treated as plain text. Paths containing spaces or shell characters are
#'   supported, and `~` is expanded to the user's home directory.
#' @param prep_chr Logical; derive chromosome-oriented `seqname` and `fctname`
#'   factor columns and order their factor levels. Rows retain file order.
#'
#' @return A data frame with one row per FASTA header. It always contains
#'   `start_line`, the numeric one-based line number of the header, and
#'   `fastaname`, the complete header without its leading `">"`, as a character
#'   column. When `prep_chr = TRUE`,
#'   `fastaname` is converted to a factor and the result also contains
#'   `seqname`, the first space-delimited field of the header, and `fctname`, the same field with
#'   `"chr"` removed; all three are factors in chromosome-aware order.
#' @seealso [get_fasta_seq_bounds()], [read_fasta()]
#' @export
#'
#' @examples
#' \dontrun{
#' fasta <- tempfile(fileext = ".fa")
#' writeLines(
#'   c(">chr2 chromosome 2", "ACGT", ">chr10 chromosome 10", "TGCA",
#'     ">chrX chromosome X", "AAAA"),
#'   fasta
#' )
#'
#' get_fasta_names(fasta)
#' get_fasta_names(fasta, prep_chr = FALSE)
#'
#' unlink(fasta)
#' }
get_fasta_names <- function(file_path,
                            prep_chr = T) {

  name_command <- if (nzchar(Sys.which("rg"))) "rg" else "grep"
  quoted_path <- shQuote(path.expand(file_path))
  if (tools::file_ext(file_path) == "gz" && name_command == "rg") {
    cmd <- paste0("rg -z -n '^>' -- ", quoted_path)
  } else if (tools::file_ext(file_path) == "gz") {
    cmd <- paste0("gunzip -c -- ", quoted_path,
                  " | grep -n '^>'")
  } else {
    cmd <- paste0(name_command, " -n '^>' -- ", quoted_path)
  }
  name_lines <- system(cmd, intern = T)

  ## factor names
  split <- strsplit(name_lines, "\\:>")
  name_lines <- stats::setNames(sapply(split, "[", 1), sapply(split, "[", 2))
  name_lines <- utils::stack(name_lines)
  names(name_lines)[1:2] <- c("start_line", "fastaname")
  name_lines$start_line <- as.numeric(name_lines$start_line)
  name_lines$fastaname <- as.character(name_lines$fastaname)

  if (prep_chr) {
    name_lines$seqname <- sapply(strsplit(name_lines$fastaname, " "), "[", 1)

    name <- gsub("chr", "", name_lines$seqname)
    name_num <- suppressWarnings(as.numeric(name))
    name_ord1 <- sort(name_num[which(!is.na(name_num))])
    name_ord2 <- name[which(is.na(name_num))]
    name_lines$fctname <- factor(name, c(name_ord1, name_ord2))

    name_lines$seqname <- factor(name_lines$seqname, levels = name_lines$seqname[order(name_lines$fctname)])
    name_lines$fastaname <- factor(name_lines$fastaname, levels = name_lines$fastaname[order(name_lines$fctname)])
  }

  return(name_lines)
}
