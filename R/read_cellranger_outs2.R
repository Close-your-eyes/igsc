#' Read  and preprocess the contig_annotation file from 10X Genomics' cellranger output and return a concatenated data frame
#'
#' This function has been reduced from read_cellranger_outs. Most important change is that CDR3 sequences
#' from equal clonotype_ids are not shared accross the clonotype. In principle this function will only read the
#' filtered_contig_annotations.csv and filtered_contig.fasta and joins them.
#'
#' @param vdj_outs_path named vector of paths to the outs-folder (or any folders containing the necessary file, specified above);
#' names will be added as "sample"-column to the output data frame
#'
#' @return data frame (cl_long)
#' @export
#'
#' @examples
read_cellranger_outs2 <- function(vdj_outs_path) {


  if (is.null(names(vdj_outs_path))) {
    stop("Please provide names with vdj_outs_path, indicating the sample name.")
  }

  ff <- c("filtered_contig_annotations.csv", "filtered_contig.fasta")

  for (i in vdj_outs_path) {
    if (any(!ff %in% list.files(i))) {
      stop(paste(ff[which(!ff %in% list.files(i))], collapse = ", "), " not found in ", i.)
    }
  }

  filt_cont_ann <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig_annotations\\.csv$", recursive = T, full.names = T)
  filt_cont_fast <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig\\.fasta$", recursive = T, full.names = T)

  contig_annotations <-
    purrr::map_dfr(filt_cont_ann, vroom::vroom, delim = ",", show_col_types = F, col_types = vroom::cols(productive = vroom::col_character()), .id = "sample") %>%
    dplyr::rename("contig_id" = contig_id, "clonotype_id" = raw_clonotype_id, "consensus_id" = raw_consensus_id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(consensus_id = gsub(paste0(clonotype_id, "_"), "", consensus_id)) %>%
    dplyr::mutate(consensus_id = gsub("_", "", consensus_id)) %>%
    dplyr::mutate(high_confidence = tolower(high_confidence), full_length = tolower(full_length), productive = tolower(productive), is_cell = tolower(is_cell))

  contig_fasta <-
    purrr::map_dfr(filt_cont_fast, function(x) utils::stack(igsc:::read_fasta(x)), .id = "sample") %>%
    dplyr::rename("contig_seq" = values, "contig_id" = ind) %>%
    tidyr::separate(contig_id, into = c("barcode", "temp", "temp2"), sep = "_", remove = F) %>%
    dplyr::select(-c(temp, temp2))

  cl_long <-
    contig_annotations %>%
    dplyr::full_join(contig_fasta, by = c("contig_id", "barcode", "sample")) %>%
    dplyr::mutate(barcode = stringr::str_replace(barcode, "-1$", ""))


  return(cl_long)
}

rcsv <- function(x) {dplyr::mutate(utils::read.csv(x, sep = ","), sample = basename(dirname(x)))}
rfasta <- function(x, vn) {dplyr::mutate(utils::stack(igsc:::read_fasta(x)), sample = basename(dirname(x))) %>% dplyr::rename(!!vn := values)}


