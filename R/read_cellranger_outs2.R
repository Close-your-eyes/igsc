#' Read necessary file from 10X Genomics' cellranger output and return a concatenated data frame
#'
#' Required files in each subfolder of vdj_outs_path are consensus_annotations.csv, filtered_contig_annotations.csv, consensus.fasta, concat_ref.fasta, filtered_contig.fasta.
#' If one or more are missing this subfolder is excluded. clonotype_id and consensus_id are suffixed with _cr for cellranger; entries from consensus.fasta becomes consensus_seq_cr.
#' ...
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

  ff <- c("concat_ref.fasta", "consensus_annotations.csv", "consensus.fasta", "filtered_contig_annotations.csv", "filtered_contig.fasta")
  pattern <- "consensus_annotations\\.csv$|filtered_contig_annotations\\.csv$|consensus\\.fasta$|concat_ref\\.fasta$|filtered_contig\\.fasta$"

  for (i in vdj_outs_path) {
    if (any(!ff %in% list.files(i))) {
      stop(paste(ff[which(!ff %in% list.files(i))], collapse = ", "), " not found in ", i.)
    }
  }


  cons_ann <- sapply(vdj_outs_path, list.files, pattern = "consensus_annotations\\.csv$", recursive = T, full.names = T)
  filt_cont_ann <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig_annotations\\.csv$", recursive = T, full.names = T)
  cons_fast <- sapply(vdj_outs_path, list.files, pattern = "consensus\\.fasta$", recursive = T, full.names = T)
  cons_ref_fast <- sapply(vdj_outs_path, list.files, pattern = "concat_ref\\.fasta$", recursive = T, full.names = T)
  filt_cont_fast <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig\\.fasta$", recursive = T, full.names = T)


  consensus_annotations <-
    purrr::map_dfr(cons_ann, utils::read.csv, sep = ",", .id = "sample") %>%
    dplyr::rename("clonotype_id_cr" = clonotype_id, "consensus_id_cr" = consensus_id)  %>%
    dplyr::mutate(consensus_id_cr = stringr::str_extract(consensus_id_cr, "[:digit:]{1,}$"))

  consensus_fasta <-
    purrr::map_dfr(cons_fast, function(x) utils::stack(igsc:::read_fasta(x)), .id = "sample") %>%
    dplyr::rename("consensus_seq_cr" = values) %>%
    tidyr::separate(ind, into = c("clonotype_id_cr", "temp", "consensus_id_cr"), sep = "_") %>%
    dplyr::select(-temp)

  consensus_ref_fasta <-
    purrr::map_dfr(cons_ref_fast, function(x) utils::stack(igsc:::read_fasta(x)), .id = "sample") %>%
    dplyr::rename("ref_seq_cr" = values) %>%
    tidyr::separate(ind, into = c("clonotype_id_cr", "temp1", "temp2", "refseq_id_cr"), sep = "_") %>%
    dplyr::select(-c(temp1, temp2))

  contig_annotations <-
    purrr::map_dfr(filt_cont_ann, utils::read.csv, sep = ",", .id = "sample") %>%
    dplyr::rename("contiq_id_cr" = contig_id, "clonotype_id_cr" = raw_clonotype_id, "consensus_id_cr" = raw_consensus_id) %>%
    dplyr::mutate(contiq_id_cr = stringr::str_extract(contiq_id_cr, "[:digit:]{1,}$"))

  contig_fasta <-
    purrr::map_dfr(filt_cont_fast, function(x) utils::stack(igsc:::read_fasta(x)), .id = "sample") %>%
    dplyr::rename("contiq_seq_cr" = values) %>%
    tidyr::separate(ind, into = c("barcode", "temp", "contiq_id_cr"), sep = "_") %>%
    dplyr::select(-c(temp))


  consensus_data <-
    consensus_annotations %>%
    dplyr::full_join(consensus_fasta, by = c("clonotype_id_cr", "consensus_id_cr", "sample")) %>%
    dplyr::full_join(consensus_ref_fasta, by = c("clonotype_id_cr" = "clonotype_id_cr", "consensus_id_cr" = "refseq_id_cr", "sample" = "sample")) %>%
    dplyr::mutate(clonotype_id_cr = stringr::str_extract(clonotype_id_cr, "[:digit:]{1,}$"))

  contig_data_barcodes <-
    contig_annotations %>%
    dplyr::full_join(contig_fasta, by = c("contiq_id_cr", "barcode", "sample")) %>%
    dplyr::distinct(sample, clonotype_id_cr, barcode) %>%
    dplyr::mutate(clonotype_id_cr = stringr::str_extract(clonotype_id_cr, "[:digit:]{1,}$")) %>%
    dplyr::mutate(barcode = stringr::str_replace(barcode, "-1$", ""))

  cl_long <-
    dplyr::left_join(consensus_data, contig_data_barcodes, by = c("clonotype_id_cr", "sample")) %>%
    dplyr::rename("V_cr" = v_gene, "D_cr" = d_gene, "J_cr" = j_gene, "C_cr" = c_gene, "CDR3_aa_cr" = cdr3, "CDR3_nt_cr" = cdr3_nt)

  return(cl_long)
}

rcsv <- function(x) {dplyr::mutate(utils::read.csv(x, sep = ","), sample = basename(dirname(x)))}
rfasta <- function(x, vn) {dplyr::mutate(utils::stack(igsc:::read_fasta(x)), sample = basename(dirname(x))) %>% dplyr::rename(!!vn := values)}


