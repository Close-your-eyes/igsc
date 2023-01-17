#' Read necessary file from 10X Genomics' cellranger output and return a concatenated data frame
#'
#' Caution: This function will add CDR3 sequences etc. uniformly to equal clonotype_ids even though
#' if e.g. TRA was detected in every cell that has been assigned with clonotype_id i.
#' Use read_cellranger_outs2 for a different behaviour, namely that this kind of
#' imputation is not performed.
#'
#' Required files in each subfolder of vdj_path are consensus_annotations.csv, filtered_contig_annotations.csv, consensus.fasta, concat_ref.fasta, filtered_contig.fasta.
#' If one or more are missing this subfolder is excluded. clonotype_id and consensus_id are suffixed with _cr for cellranger; entries from consensus.fasta becomes consensus_seq_cr.
#' ...
#'
#' @param vdj_path path to the parent folder containing one or multiple subfolders with cellrangers VDJ outs
#'
#' @return data frame (cl_long)
#' @export
#'
#' @examples
read_cellranger_outs <- function(vdj_path) {

  vdj_dirs <- list.dirs(vdj_path, recursive = F)
  if (length(vdj_dirs) > 1) {vdj_dirs <- vdj_dirs[-1]}

  ff <- c("concat_ref.fasta", "consensus_annotations.csv", "consensus.fasta", "filtered_contig_annotations.csv", "filtered_contig.fasta")
  pattern <- "consensus_annotations\\.csv$|filtered_contig_annotations\\.csv$|consensus\\.fasta$|concat_ref\\.fasta$|filtered_contig\\.fasta$"

  vdj_dirs <- sapply(vdj_dirs, USE.NAMES = F, function(x) {
    if (length(list.files(x, pattern = pattern)) == 5) {
      return(x)
    } else if (length(list.files(x, pattern = pattern) < 5)) {
      f <- list.files(x, pattern = pattern)
      print(paste0(paste(ff[which(!ff %in% f)], collapse = ", "), " not found in ", basename(x)))
      return(NULL)
    } else {
      f <- list.files(x, pattern = pattern)
      print(paste0(paste(ff, collapse = ", "), " not uniquely identified in ", basename(x), ":"))
      print(f)
      return(NULL)
    }
  })

  cons_ann <- list.files(vdj_path, "consensus_annotations\\.csv$", recursive = T, full.names = T)
  filt_cont_ann <- list.files(vdj_path, "filtered_contig_annotations\\.csv$", recursive = T, full.names = T)
  cons_fast <- list.files(vdj_path, "consensus\\.fasta$", recursive = T, full.names = T)
  cons_ref_fast <- list.files(vdj_path, "concat_ref\\.fasta$", recursive = T, full.names = T)
  filt_cont_fast <- list.files(vdj_path, "filtered_contig\\.fasta$", recursive = T, full.names = T)

  consensus_annotations <-
    dplyr::bind_rows(lapply(cons_ann, rcsv)) %>%
    dplyr::rename("clonotype_id_cr" = clonotype_id, "consensus_id_cr" = consensus_id)  %>%
    dplyr::mutate(consensus_id_cr = stringr::str_extract(consensus_id_cr, "[:digit:]{1,}$"))

  consensus_fasta <-
    dplyr::bind_rows(lapply(cons_fast, rfasta, vn = "consensus_seq_cr")) %>%
    tidyr::separate(ind, into = c("clonotype_id_cr", "temp", "consensus_id_cr"), sep = "_") %>%
    dplyr::select(-temp)

  consensus_ref_fasta <-
    dplyr::bind_rows(lapply(cons_ref_fast, rfasta, vn = "ref_seq_cr")) %>%
    tidyr::separate(ind, into = c("clonotype_id_cr", "temp1", "temp2", "refseq_id_cr"), sep = "_") %>%
    dplyr::select(-c(temp1, temp2))

  contig_annotations <-
    dplyr::bind_rows(lapply(filt_cont_ann, rcsv)) %>%
    dplyr::rename("contiq_id_cr" = contig_id, "clonotype_id_cr" = raw_clonotype_id, "consensus_id_cr" = raw_consensus_id) %>%
    dplyr::mutate(contiq_id_cr = stringr::str_extract(contiq_id_cr, "[:digit:]{1,}$"))

  contig_fasta <-
    dplyr::bind_rows(lapply(filt_cont_fast, rfasta, vn = "contiq_seq_cr")) %>%
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
rfasta <- function(x, vn) {dplyr::mutate(utils::stack(read_fasta(x)), sample = basename(dirname(x))) %>% dplyr::rename(!!vn := values)}


