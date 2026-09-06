#' Rotate a circular genome and its GTF annotations
#'
#' Detect whether annotations wrap across the origin of a circular genome and,
#' when necessary, rotate the genome sequence and annotation coordinates to a
#' shared cut position. This keeps the sequence and its GTF records in the same
#' coordinate system.
#'
#' The genome length is calculated with [base::nchar()] and supplied to
#' [detect_wrap_and_cut()]. If that function returns a cut, the GTF coordinates
#' and genome string are rotated together. If no cut is detected, both remain in
#' their original coordinate system. In either case, the GTF attribute column
#' is converted to key-value form before the result is returned.
#'
#' @param genome A single character string containing the complete circular
#'   genome sequence.
#' @param gtf A data frame containing GTF annotations for `genome`, in the format
#'   accepted by [detect_wrap_and_cut()]. Its coordinates must refer to the
#'   supplied, unrotated genome sequence.
#' @param verbose Logical; emit informational messages while detecting a cut and
#'   preparing the GTF attribute column.
#'
#' @return A named list with two elements: `genome`, the original or rotated
#'   genome string, and `gtf`, the corresponding GTF data frame with original or
#'   rotated coordinates and a key-value attribute column.
#' @export
rotate_genome_and_gtf <- function(genome,
                                  gtf,
                                  verbose = T) {

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
