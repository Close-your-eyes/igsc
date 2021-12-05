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
                     TRAxTRB = T,
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
    read.csv("/Users/vonskopnik/Documents/scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv", sep = "\t", header = T) %>%
    dplyr::filter(Species == "HomoSapiens") %>%
    dplyr::distinct(!!sym(vdj_cdr3_col), !!sym(vdj_tr_col))

  vdjdb <- vdjdb %>% dplyr::group_by(Gene) %>% dplyr::slice(1:10)
  tcrs <- tcrs %>% dplyr::group_by(chain) %>% dplyr::slice(1:10)

  if (!identical(sort(unique(vdjdb[,vdj_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(vdj_tr_col, " column of vdjdb should only contain TRA and/or TRB."))
  }
  if (!identical(sort(unique(tcrs[,tcr_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(tcr_tr_col, " column of tcrs should only contain TRA and/or TRB."))
  }

  matches <- mapply(stringdist::stringdistmatrix,
                    split(vdjdb[,vdj_cdr3_col,drop=T], vdjdb[,vdj_tr_col,drop=T])[sort(unique(vdjdb[,vdj_tr_col,drop=T]))],
                    split(tcrs[,tcr_cdr3_col,drop=T], tcrs[,tcr_tr_col,drop=T])[sort(unique(tcrs[,tcr_tr_col,drop=T]))],
                    method = "lv",
                    useNames = "strings",
                    SIMPLIFY = F)
  names(matches) <- paste(sort(unique(vdjdb[,vdj_tr_col,drop=T])), sort(unique(tcrs[,tcr_tr_col,drop=T])), sep = "_")
  matches <- sapply(matches, reshape2::melt, simplify = F, c(vdj_cdr3_col, tcr_cdr3_col), value.name = "lv")
  matches <- do.call(rbind, matches)
  matches[,vdj_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 1)
  matches[,tcr_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 2)
  rownames(matches) <- NULL

  if (TRAxTRB) {
    matches2 <- mapply(stringdist::stringdistmatrix,
                      split(vdjdb[,vdj_cdr3_col,drop=T], vdjdb[,vdj_tr_col,drop=T])[rev(sort(unique(vdjdb[,vdj_tr_col,drop=T])))],
                      split(tcrs[,tcr_cdr3_col,drop=T], tcrs[,tcr_tr_col,drop=T])[sort(unique(tcrs[,tcr_tr_col,drop=T]))],
                      method = "lv",
                      useNames = "strings",
                      SIMPLIFY = F)
    names(matches2) <- paste(rev(sort(unique(vdjdb[,vdj_tr_col,drop=T]))), sort(unique(tcrs[,tcr_tr_col,drop=T])), sep = "_")
    matches2 <- sapply(matches2, reshape2::melt, simplify = F, c(vdj_cdr3_col, tcr_cdr3_col), value.name = "lv")
    matches2 <- do.call(rbind, matches2)
    matches2[,vdj_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches2), "\\."), "[", 1), "_"), "[", 1)
    matches2[,tcr_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches2), "\\."), "[", 1), "_"), "[", 2)
    rownames(matches2) <- NULL
    matches <- rbind(matches, matches2)
  }




}
