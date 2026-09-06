#' Write a GTF data frame to a file
#'
#' Write the rows of a data frame as tab-delimited GTF records, optionally
#' processing the attribute column for identifier uniqueness and compressing the
#' result. Header lines are written before the records, and missing parent
#' directories are created recursively.
#'
#' Columns are written in their existing order; the function does not rearrange
#' or validate them against the nine-column GTF specification. Consequently,
#' `gtf_df` must include `start` and `end` columns and should already have the
#' intended GTF column order.
#'
#' When `gzip = TRUE`, the external `gzip` command compresses `file`, removes
#' the uncompressed file, and creates a file whose path is `paste0(file,
#' ".gz")`.
#'
#' @param gtf_df A data frame containing GTF records. It must contain `start`
#'   and `end` columns; all columns are written in their current order.
#' @param file Output path for the uncompressed GTF file. Its parent directory
#'   is created when necessary.
#' @param header Character vector of header lines to prepend to the file. Supply
#'   `character(0)` to write no header.
#' @param check_unique Logical; process the GTF attribute column with
#'   [process_gtf_attribute_col()] before writing, ensuring that `gene_name`,
#'   `gene_id`, and `transcript_id` values are unique where required.
#' @param gzip Logical; compress the completed file with the external `gzip`
#'   command.
#' @param verbose Logical; emit progress messages while processing attributes
#'   and report the requested output path after writing.
#'
#' @return `NULL`, invisibly. This function is called for its file-writing side
#'   effect.
#' @export
#'
#' @examples
#' gtf <- data.frame(
#'   seqname = c("chr1", "chr1"),
#'   source = "example",
#'   feature = c("gene", "transcript"),
#'   start = c(1, 1),
#'   end = c(100, 100),
#'   score = ".",
#'   strand = "+",
#'   frame = ".",
#'   attribute = c(
#'     'gene_id "g1"; gene_name "G1";',
#'     'gene_id "g1"; transcript_id "t1"; gene_name "G1";'
#'   )
#' )
#'
#' gtf_file <- tempfile(fileext = ".gtf")
#' write_gtf(
#'   gtf_df = gtf,
#'   file = gtf_file,
#'   check_unique = FALSE,
#'   verbose = FALSE
#' )
#' readLines(gtf_file)
#' unlink(gtf_file)
write_gtf <- function(gtf_df,
                      file,
                      header = "##gtf file",
                      check_unique = T,
                      gzip = F,
                      verbose = T) {

  if (missing(file)) {
    stop("file missing.")
  }
  if (missing(gtf_df)) {
    stop("gtf_df missing.")
  }

  dir.create(dirname(file), recursive = T, showWarnings = F)

  gtf_df <- gtf_df |>
    dplyr::mutate(start = as.character(start), end = as.character(end)) |>
    tibble::as_tibble()

  if (check_unique) {
    # fix duplicate gene_name, gene_id, transcript_id
    gtf_df <- process_gtf_attribute_col(
      gtf_df,
      attr_as = "kv",
      verbose = verbose)[["gtf"]]
  }

  out <- vroom::vroom_write_lines(c(header,
                                    apply(gtf_df, 1, paste, collapse = "\t")),
                                  file = file)
  if (gzip) {
    system(paste0("gzip ", file))
  }
  if (verbose) {
    message(file)
  }
}
