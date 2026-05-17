#' Write a gtf data frames to file
#'
#' @param gtf_df gtf data frame
#' @param header vector of header lines
#' @param file file path where to save
#' @param check_unique check for uniqueness of gene_name, gene_id,
#' transcript_id
#' @param verbose
#' @param gzip
#'
#' @returns
#' @export
#'
#' @examples
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
