#' Title
#'
#' @param gtf_file
#'
#' @return
#' @export
#'
#' @examples
genes_gtf_to_fst <- function(gtf_file,
                             subfolder = "genes_gtf_fst") {
  path <- file.path(dirname(gtf_file), subfolder)
  seqnames <-
    read_gtf(file_path = gtf_file, process_attr_col = F)[["gtf"]] %>%
    dplyr::distinct(seqname) %>%
    dplyr::pull(seqname)
  # separate fst file for each seqname from gtf file
  purrr::map(seqname, function(x) {
    y <- read_gtf(file_path = gtf_file, seqnames = x)[["gtf"]]
    fst::write_fst(y, path = file.path(path, paste0("genes_gft_", x, ".fst")), compress = 100)
    return(NULL)
  }, .progress = T)
  # one file with info which gene is on which chr
  fst_files <- list.files(path, full.names = T)
  genes_chr <- purrr::map_dfr(fst_files, function(x) dplyr::distinct(fst::read_fst(path = x, columns = c("seqname", "gene_name"))))
  fst::write_fst(genes_chr, path = file.path(dirname(path), "genes_chr.fst"), compress = 100)
}
