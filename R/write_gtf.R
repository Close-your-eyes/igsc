#' Title
#'
#' @param gtf_df_list
#' @param header
#' @param file
#'
#' @returns
#' @export
#'
#' @examples
write_gtf <- function(gtf_df_list,
                      header = "##gtf file",
                      file) {
  if (missing(file)) {
    stop("file missing.")
  }
  if (is.list(gtf_df_list)) {
    gtf_df_list <- purrr::reduce(sapply(gtf_df_list, "[", 2), dplyr::bind_rows)
  }
  gtf_df_list <- gtf_df_list |>
    dplyr::mutate(start = as.character(start), end = as.character(end)) |>
    tibble::as_tibble()

  vroom::vroom_write_lines(c(header,
                             apply(gtf_df_list, 1, paste, collapse = "\t")),
                           file = file)
}
