#' Title
#'
#' @param tcrs data.frame of TCRs from sequencing, first column indicating the chain TRA or TRB, second the aa sequence of CDR3
#' @param vdjdb data.frame of TCRs from vdjdb, first column indicating the chain TRA or TRB, second the aa sequence of CDR3
#'
#' @return
#' @export
#'
#' @examples
vdj_hits <- function(vdjdb,
                     tcrs,
                     vdj_tr_col = "Gene",
                     tcr_tr_col = "chain",
                     vdj_cdr3_col = "CDR3",
                     tcr_cdr3_col = "CDR3_aa_cr") {

  # ... stringdist nthread
  tcrs <-
    cl_long[,which(names(cl_long) %in% c("chain", "CDR3_aa_cr"))] %>%
    dplyr::distinct()

  vdjdb <-
    read.csv("/Users/vonskopnik/Documents/scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv", sep = "\t", header = T) %>%
    dplyr::filter(Species == "HomoSapiens") %>%
    dplyr::distinct(CDR3, Gene)

  # Biostrings pDict approach and similar are not applicable

  vdjdb <- vdjdb %>% dplyr::group_by(Gene) %>% dplyr::slice(1:5)
  tcrs <- tcrs %>% dplyr::group_by(chain) %>% dplyr::slice(1:5)
  matches <- mapply(stringdist::stringdistmatrix, split(vdjdb$CDR3, vdjdb$Gene), split(tcrs$CDR3_aa_cr, tcrs$chain), method = "lv", useNames = "strings", SIMPLIFY = F)
  matches <- sapply(matches, reshape2::melt, simplify = F, c("CDR3", "CDR3_aa_cr"), value.name = "lv")
  matches <- do.call(rbind, matches)
  matches$TR <- sapply(strsplit(rownames(matches), "\\."), "[", 1)
  rownames(matches) <- NULL


}
