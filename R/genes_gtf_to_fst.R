#' Title
#'
#' @param gtf_file
#' @param subfolder
#' @param compression
#'
#' @return
#' @export
#'
#' @examples
genes_gtf_to_fst <- function(gtf_file,
                             subfolder = "genes_fst",
                             compression = 50) {
  if (missing(gtf_file) || length(gtf_file) == 0) {
    stop("path to gtf_file missing.")
  }
  if (!file.exists(gtf_file)) {
    stop("gtf_file not found.")
  }
  path <- file.path(dirname(gtf_file), subfolder)
  dir.create(path, showWarnings = F)
  seqnames <-
    read_gtf(file_path = gtf_file, process_attr_col = F)[["gtf"]] %>%
    dplyr::distinct(seqname) %>%
    dplyr::pull(seqname)
  # separate fst file for each seqname from gtf file
  # one file with info which gene is on which chr
  genes_chr <- purrr::map_dfr(seqnames, function(x) {
    print(x)
    y <- read_gtf(file_path = gtf_file, seqnames = x)[["gtf"]]
    fst::write_fst(y, path = file.path(path, paste0(x, "_gtf.fst")), compress = compression)
    return(dplyr::distinct(y, seqname, gene_id, gene_name))
  })
  fst::write_fst(genes_chr, path = file.path(dirname(gtf_file), "genes_on_chr.fst"), compress = compression)
  if (tools::file_ext(gtf_file) != "gz") {
    system(paste0("gzip -1 ", gtf_file))
  }
}
