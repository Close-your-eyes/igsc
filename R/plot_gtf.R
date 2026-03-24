#' Plot ranges from a GTF file to see genome organization
#'
#' @param gtf_df data frame from igsc::read_gtf
#' @param y y axis column
#' @param color color column; overlapping is computed internally
#'
#' @returns list of ggplot and data frame
#' @export
#'
#' @examples
#'\dontrun{
#' gtf <- igsc::read_gtf(gtf_path)$gtf |>
#'   dplyr::filter(feature == "CDS")
#' plot_gtf(gtf, color = "overlapping")
#' }
plot_gtf <- function(gtf_df,
                     y = "transcript_id",
                     color = c("strand", "overlapping")) {

  if ("seqname" %in% names(gtf_df) && length(unique(gtf_df[["seqname"]]))>1) {
    stop("more than one seqname found. makes not sense.")
  }

  if (!"start" %in% names(gtf_df) || !"end" %in% names(gtf_df)) {
    stop("start and end column needed in gtf_df.")
  }

  color <- rlang::arg_match(color)

  # start and end are fixed
  gtf_df[[y]] <- forcats::fct_reorder(gtf_df[[y]], gtf_df[["start"]])

  if (color == "strand") {
    if (!color %in% names(gtf_df)) {
      stop("color column needed in gtf_df.")
    }
    plot <- ggplot2::ggplot(gtf_df) +
      ggplot2::geom_segment(ggplot2::aes(
        x = start,
        xend = end,
        y = !!rlang::sym(y),
        yend = !!rlang::sym(y),
        color = !!rlang::sym(color)
      )) +
      colrr::theme_material() +
      ggplot2::labs(x = "position")

  } else if (color == "overlapping") {

    ## separate grouping column for streaks of overlapping T/F needed
    gtf_df[[y]] <- as.character(gtf_df[[y]])
    gtf_df$row <- make.unique(gtf_df[[y]])

    gtf_df2 <- utils::stack(brathering::seq2(stats::setNames(gtf_df$start, gtf_df[["row"]]), gtf_df$end)) |>
      dplyr::add_count(values) |>
      dplyr::mutate(overlapping = n > 1) |>
      dplyr::arrange(ind, values) |>
      dplyr::mutate(ind = as.character(ind)) |>
      dplyr::left_join(gtf_df[,c("row", y)], by = c("ind" = "row"))

    r <- rle(gtf_df2$overlapping)
    gtf_df2$ovlp_run <- rep(seq_along(r$lengths), r$lengths)
    gtf_df2$group <- paste0(gtf_df2$ind, "__", gtf_df2$ovlp_run)
    temp <- dplyr::slice_min(gtf_df2, order_by = values, n = 1, by = y) |>
      dplyr::arrange(values)
    gtf_df2[[y]] <- factor(gtf_df2[[y]], level = temp[[y]])

    plot <- ggplot2::ggplot(gtf_df2, ggplot2::aes(x = values, y = !!rlang::sym(y), color = !!rlang::sym(color))) +
      ggplot2::geom_line(ggplot2::aes(group = group)) +
      colrr::theme_material() +
      ggplot2::labs(x = "position", y = y)
    gtf_df <- gtf_df2
  }

  return(list(plot = plot, data = gtf_df))
}
