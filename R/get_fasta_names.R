#' Extract sequence names from a FASTA file
#'
#' Find every FASTA header in a plain-text or gzip-compressed file and return
#' its line number and header text. The leading `">"` is removed from each
#' header.
#'
#' When `prep_chr = TRUE`, the first whitespace-delimited field of each header
#' is treated as a sequence name. The `"chr"` substring is removed to construct
#' an ordering key: numeric names are ordered numerically, followed by
#' non-numeric names in their input order. This is useful for chromosome FASTA
#' files, where lexical ordering would otherwise place `chr10` before `chr2`.
#'
#' Header detection uses `rg` when available and falls back to `grep`. Reading a
#' file whose extension is `.gz` additionally requires the `gunzip` command.
#'
#' @param file_path Path to a FASTA file. Files with a `.gz` extension are read
#'   as gzip-compressed input; all other extensions are read as plain text.
#' @param prep_chr Logical; derive chromosome-oriented `seqname` and `fctname`
#'   factor columns and order the FASTA names by them.
#'
#' @return A data frame with one row per FASTA header. It always contains
#'   `start_line`, the numeric line number of the header, and `fastaname`, the
#'   complete header without its leading `">"`. When `prep_chr = TRUE`,
#'   `fastaname` is converted to a factor and the result also contains
#'   `seqname`, the first field of the header, and `fctname`, the same field with
#'   `"chr"` removed; all three are factors in chromosome-aware order.
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

  name_lines <- tryCatch(
    {
      if (tools::file_ext(file_path) == "gz") {
        cmd <- paste0("gunzip -c ", file_path, " | rg -n '^>' ")
      } else {
        cmd <- paste0("rg -n '^>' ", file_path)
      }
      system(cmd, intern = T)
    },
    error = function(err) {
      if (tools::file_ext(file_path) == "gz") {
        cmd <- paste0("gunzip -c ", file_path, " | grep -n '^>' ")
      } else {
        cmd <- paste0("grep -n '^>' ", file_path)
      }
      system(cmd, intern = T)
    }
  )

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
