read_cellranger_outs <- function(vdj_path) {

  vdj_path <- "/Volumes/AG_Hiepe/Christopher.Skopnik/2019_scRNAseq/R_scRNAseq/2019_SLE_LN/data/Sequencing_data/arranged_for_Seurat/filtered_feature_bc_matrix/VDJ"

  vdj_dirs <- list.dirs(vdj_path, recursive = F)
  if (length(vdj_dirs) > 1) {
    vdj_dirs <- vdj_dirs[-1]
  }
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

  rcsv <- function(x) {dplyr::mutate(read.csv(x, sep = ","), sample = basename(dirname(x)))}
  rfasta <- function(x, vn) {dplyr::mutate(stack(read_fasta(x)), sample = basename(dirname(x))) %>% dplyr::rename(!!vn := values)}

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
    dplyr::rename("Barcode" = barcode, "contiq_id_cr" = contig_id, "clonotype_id_cr" = raw_clonotype_id, "consensus_id_cr" = raw_consensus_id) %>%
    dplyr::mutate(contiq_id_cr = stringr::str_extract(contiq_id_cr, "[:digit:]{1,}$"))

  contig_fasta <-
    dplyr::bind_rows(lapply(filt_cont_fast, rfasta, vn = "contiq_seq_cr")) %>%
    tidyr::separate(ind, into = c("Barcode", "temp", "contiq_id_cr"), sep = "_") %>%
    dplyr::select(-c(temp))


  consensus_data <-
    consensus_annotations %>%
    dplyr::full_join(consensus_fasta, by = c("clonotype_id_cr" = "clonotype_id_cr", "consensus_id_cr" = "consensus_id_cr", "sample" = "sample")) %>%
    dplyr::full_join(consensus_ref_fasta, by = c("clonotype_id_cr" = "clonotype_id_cr", "consensus_id_cr" = "refseq_id_cr", "sample" = "sample")) %>%
    dplyr::mutate(clonotype_id_cr = stringr::str_extract(clonotype_id_cr, "[:digit:]{1,}$"))

  contig_data_barcodes <-
    contig_annotations %>%
    dplyr::full_join(contig_fasta, by = c("contiq_id_cr" = "contiq_id_cr", "Barcode" = "Barcode", "sample" = "sample")) %>%
    dplyr::distinct(sample, clonotype_id_cr, Barcode) %>%
    dplyr::mutate(clonotype_id_cr = stringr::str_extract(clonotype_id_cr, "[:digit:]{1,}$")) %>%
    dplyr::mutate(Barcode = stringr::str_replace(Barcode, "-1$", ""))

  cl_long <-
    dplyr::left_join(consensus_data, contig_data_barcodes, by = c("clonotype_id_cr", "sample")) %>%
    dplyr::mutate(clonotype_id_cr = paste0(sample, "_", clonotype_id_cr)) %>%
    dplyr::rename("V_cr" = v_gene, "D_cr" = d_gene, "J_cr" = j_gene, "C_cr" = c_gene, "CDR3_aa_cr" = cdr3, "CDR3_nt_cr" = cdr3_nt)

  return(cl_long)
}
