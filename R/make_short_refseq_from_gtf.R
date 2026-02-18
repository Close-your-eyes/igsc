#' Make shorter refseqs from gtf subset and original full refseqs
#'
#' @param gtf_df subsetted gtf data frame
#' @param refseqs named ref sequence character vector; names must be in
#' gtf_df$seqname.
#' @param overhang_whole how much overhang in refseq
#'
#' @returns
#' @export
#'
#' @examples
#' \dontrun{
#' out <- make_short_refseq_from_gtf(hugtf, genome)
#' }
make_short_refseq_from_gtf <- function(gtf_df,
                                       refseqs,
                                       overhang_whole = c(0,0)) {

  if (any(!c("seqname", "start", "end", "gene_name") %in% names(gtf_df))) {
    stop("gtf_df has to have the columns at least: seqname, start, end, gene_name.")
  }
  if (is.null(names(refseqs))) {
    stop("names(refseqs) cannot be NULL.")
  }
  if (any(!unique(gtf_df$seqname) %in% names(refseqs))) {
    stop("not all seqname from gtf_df in names(refseqs).")
  }

  gtf_df <- dplyr::mutate(gtf_df, seqname2 = paste0(seqname, "_", gene_name))

  refseqs <- purrr::map_chr(split(gtf_df, gtf_df$seqname2), ~get_seq_from_refseq(
    refseq = refseqs[unique(.x$seqname)],
    gtf_df = .x,
    #names = paste0(unique(.x$seqname), "_", unique(.x$gene_name)),
    what = "whole",
    overhang_whole = overhang_whole
  ))

  ## correct coordinates
  gtf_df <- purrr::map_dfr(split(gtf_df, gtf_df$seqname2), function(x) {
    min_coord <- min(x$start, x$end) - 1 - overhang_whole[1]
    x$start <- x$start-min_coord
    x$end <- x$end-min_coord
    return(x)
  })

  gtf_df <- gtf_df |>
    dplyr::mutate(seqname = seqname2) |>
    dplyr::select(-seqname2)

  message("use write_gtf and write_fasta to save gtf and refseqs.")

  return(list(gtf_df = gtf_df, refseqs = refseqs))
}
