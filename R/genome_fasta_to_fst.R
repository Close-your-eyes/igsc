#' Title
#'
#' @param genome_file
#' @param out_folder
#' @param overwrite
#' @param fst_compression
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
genome_fasta_to_fst <- function(genome_file,
                                out_folder = file.path(dirname(genome_file), "genome_fa_fst"),
                                overwrite = F,
                                fst_compression = 100,
                                verbose = T,
                                simplify_names = T) {
  seq_bounds <- get_fasta_seq_bounds(genome_file)
  dir.create(out_folder, showWarnings = F)

  for (x in seq_along(seq_bounds$name)) {
    if (simplify_names) {
      name <- strsplit(seq_bounds$name[x], " ")[[1]][1]
    } else {
      name <- seq_bounds$name[x]
    }

    out_path <- file.path(out_folder, paste0(name, ".fst"))
    if (verbose) {
      message(name, " ", x, "/", length(seq_bounds$name), " ≈", format(seq_bounds$bp_approx[x], big.mark = ","), " bp")
    }
    skip <- file.exists(out_path)
    if (!skip && overwrite) {
      skip <- F
    }
    if (!skip) {
      #ind <- which(seq_bounds$name == name)
      refseq <- read_fasta(genome_file,
                           start_line = seq_bounds[x, "start_line"],
                           end_line = seq_bounds[x, "end_line"])

      refseq <- data.frame(x = strsplit(x = refseq, split = "", fixed = T)[[1]])
      names(refseq) <- name
      fst::write_fst(refseq, path = out_path, compress = fst_compression)
    }
  }
}
