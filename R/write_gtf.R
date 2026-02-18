#' Write a gtf data frames to file
#'
#' @param gtf_df gtf data frame
#' @param header header lines
#' @param file file path where to save
#'
#' @returns
#' @export
#'
#' @examples
write_gtf <- function(gtf_df,
                      file,
                      header = "##gtf file") {

  if (missing(file)) {
    stop("file missing.")
  }

  dir.create(dirname(file), recursive = T, showWarnings = F)

  gtf_df <- gtf_df |>
    dplyr::mutate(start = as.character(start), end = as.character(end)) |>
    tibble::as_tibble()

  out <- vroom::vroom_write_lines(c(header,
                                    apply(gtf_df, 1, paste, collapse = "\t")),
                                  file = file)
}
