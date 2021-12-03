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
                     tcr_cdr3_col = "CDR3_aa_cr",
                     separate_TRA_TRB_only = T,
                     ...) { # ... stringdist nthread

  if (vdj_tr_col == tcr_tr_col) {
    stop("Please avoid same names for tcr_tr_col and vdj_tr_col.")
  }
  if (vdj_cdr3_col == tcr_cdr3_col) {
    stop("Please avoid same names for tcr_cdr3_col and vdj_cdr3_col")
  }


  tcrs <- cl_long
  tcrs <- dplyr::distinct(tcrs, !!sym(tcr_tr_col), !!sym(tcr_cdr3_col))

#/Volumes/AG_Hiepe/Christopher.Skopnik/2019_scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv
  #"/Users/vonskopnik/Documents/scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv"
  vdjdb <-
    read.csv("/Volumes/AG_Hiepe/Christopher.Skopnik/2019_scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv", sep = "\t", header = T) %>%
    dplyr::filter(Species == "HomoSapiens") %>%
    dplyr::distinct(!!sym(vdj_cdr3_col), !!sym(vdj_tr_col))

  vdjdb <- vdjdb %>% dplyr::group_by(Gene) %>% dplyr::slice(1:5)
  tcrs <- tcrs %>% dplyr::group_by(chain) %>% dplyr::slice(1:5)
  if (separate_TRA_TRB_only) {
    matches <- mapply(stringdist::stringdistmatrix,
                      split(vdjdb$CDR3, vdjdb$Gene)[c("TRA", "TRB")],
                      split(tcrs$CDR3_aa_cr, tcrs$chain)[c("TRA", "TRB")],
                      method = "lv",
                      useNames = "strings",
                      SIMPLIFY = F,
                      USE.NAMES = T)
    matches <- sapply(matches, reshape2::melt, simplify = F, c("CDR3", "CDR3_aa_cr"), value.name = "lv")
    matches <- do.call(rbind, matches)
    matches[,vdj_tr_col] <- sapply(strsplit(rownames(matches), "\\."), "[", 1)
    matches[,tcr_tr_col] <- sapply(strsplit(rownames(matches), "\\."), "[", 1)
  } else {

  }

  split(vdjdb$CDR3, vdjdb$Gene)


  rownames(matches) <- NULL


}
