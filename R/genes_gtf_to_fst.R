#' Split GTF annotations into per-sequence fst files
#'
#' Partition a GTF file by `seqname` and write each partition as a separate
#' [fst::write_fst()] file. A second lookup file records the distinct
#' combinations of sequence name, gene ID, and gene name.
#'
#' Per-sequence files are written to `subfolder` beneath the directory containing
#' `gtf_file` and are named `<seqname>_gtf.fst`. The lookup table is written as
#' `genes_on_chr.fst` beside the input GTF, not inside `subfolder`. Existing
#' output files may be replaced by [fst::write_fst()].
#'
#' After all outputs are written, an input file whose extension is not `.gz` is
#' compressed with the external command `gzip -1`. On success, this creates
#' `paste0(gtf_file, ".gz")` and removes the original uncompressed GTF file.
#'
#' @param gtf_file Path to the source GTF file. It must be readable by
#'   [read_gtf()] and provide `seqname`, `gene_id`, and `gene_name` columns after
#'   its attribute column is processed.
#' @param subfolder Name of the output directory to create beneath
#'   `dirname(gtf_file)` for the per-sequence files.
#' @param compression Numeric fst compression level passed as `compress` to
#'   [fst::write_fst()], conventionally from `0` to `100`.
#'
#' @return Invisibly, the exit status returned by `gzip` when an uncompressed
#'   input is compressed; otherwise `NULL`. This function is primarily called
#'   for its file-writing side effects.
#' @export
#'
#' @examples
#' \dontrun{
#' genes_gtf_to_fst(
#'   gtf_file = "reference/genes.gtf",
#'   subfolder = "genes_fst",
#'   compression = 50
#' )
#'
#' list.files("reference/genes_fst", pattern = "[.]fst$")
#' }
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
    read_gtf(file_path = gtf_file, process_attr_col = F)[["gtf"]] |>
    dplyr::distinct(seqname) |>
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
