#' Title
#'
#' @param genome
#' @param gtf
#'
#' @returns
#' @export
#'
#' @examples
rotate_genome_and_gtf <- function(genome, gtf, verbose = T) {

  genome_length = nchar(genome)

  # cut <- pick_best_cut(gtf, genome_length = genome_length) # igsc:::
  cut <- detect_wrap_and_cut(df = gtf, genome_length = genome_length, verbose = verbose)
  if (!is.null(cut)) {
    gtf <- igsc:::rotate_coords(gtf, cut = cut$cut_position, genome_length = genome_length)
    # gtf <- fix_duplicate_rows(gtf) # done in rotate_coords
    genome <- igsc:::rotate_genome_string(genome = genome, cut = cut$cut_position)
  }

  gtf <- igsc:::make_kv_attr_col(gtf, verbose = verbose)

  return(list(genome = genome, gtf = gtf))
}
