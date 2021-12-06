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
                     TRAxTRB = F,
                     max_lvdist = 3,
                     PAM30_similarity = F,
                     mapply_fun = mapply,
                     nthread = getOption("sd_num_thread"),
                     ...) {

  if (vdj_tr_col == tcr_tr_col) {
    stop("Please avoid same names for tcr_tr_col and vdj_tr_col.")
  }
  if (vdj_cdr3_col == tcr_cdr3_col) {
    stop("Please avoid same names for tcr_cdr3_col and vdj_cdr3_col")
  }
  if (!identical(sort(unique(vdjdb[,vdj_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(vdj_tr_col, " column of vdjdb should only contain TRA and/or TRB."))
  }
  if (!identical(sort(unique(tcrs[,tcr_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(tcr_tr_col, " column of tcrs should only contain TRA and/or TRB."))
  }

  mapply_fun <- match.fun(mapply_fun)

  tcrs <- cl_long
  tcrs <- dplyr::distinct(tcrs, !!sym(tcr_tr_col), !!sym(tcr_cdr3_col))

  #/Volumes/AG_Hiepe/Christopher.Skopnik/2019_scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv
  #"/Users/vonskopnik/Documents/scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv"
  vdjdb <-
    read.csv("/Users/vonskopnik/Documents/scRNAseq/R_scRNAseq/2019_SLE_LN/data/20211020_vdjdb.tsv", sep = "\t", header = T) %>%
    dplyr::filter(Species == "HomoSapiens") %>%
    dplyr::distinct(!!sym(vdj_cdr3_col), !!sym(vdj_tr_col))

  vdjdb <- vdjdb %>% dplyr::group_by(Gene) %>% dplyr::slice(1:200)
  tcrs <- tcrs %>% dplyr::group_by(chain) %>% dplyr::slice(1:200)



  matches <- do.call(rbind, lapply(unique(c(F, TRAxTRB)),
                                   vdjdb_tcrs_match_fun,
                                   vdjdb = vdjdb,
                                   tcrs = tcrs,
                                   vdj_tr_col = vdj_tr_col,
                                   tcr_tr_col = tcr_tr_col,
                                   vdj_cdr3_col = vdj_cdr3_col,
                                   tcr_cdr3_col = tcr_cdr3_col,
                                   nthread = nthread))
  matches <- matches[which(matches$lv <= max_lvdist),]

  if (PAM30_similarity) {
    matches$PAM30 <- mapply_fun(Biostrings::pairwiseAlignment,
                                matches[,vdj_cdr3_col],
                                matches[,tcr_cdr3_col],
                                substitutionMatrix = "PAM30",
                                scoreOnly = T,
                                ...)
    ref1 <- mapply_fun(Biostrings::pairwiseAlignment,
                       matches[,vdj_cdr3_col],
                       matches[,vdj_cdr3_col],
                       substitutionMatrix = "PAM30",
                       scoreOnly = T,
                       ...)
    ref2 <- mapply_fun(Biostrings::pairwiseAlignment,
                       matches[,tcr_cdr3_col],
                       matches[,tcr_cdr3_col],
                       substitutionMatrix = "PAM30",
                       scoreOnly = T,
                       ...)
    matches$PAM30_corr <- matches$PAM30/(mapply(function(x, y) mean(c(x, y)), ref1, ref2))
  }


  matches2 <- dplyr::left_join(matches, tcrs) %>% dplyr::left_join(vdjdb) %>% dplyr::distinct(CDR3, CDR3_aa_cr, lv, Gene, chain, clonotype_id_cr, sample,
                                                                                              Species, Epitope, Epitope.gene, Epitope.species, Score)

  ggplot(matches2, aes(x = chain, y = lv, fill = Epitope.species, color = Gene)) +
    geom_jitter(width = 0.15, height = 0.15, shape = 21) +
    scale_y_continuous(breaks = fcexpr::int_breaks(n = 3)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = c("#000000", "#999999")) +
    facet_grid(rows = vars(Score), cols = vars(clonotype_id_cr))



}

vdjdb_tcrs_match_fun <- function(sort_desc,
                                 vdjdb,
                                 tcrs,
                                 vdj_tr_col,
                                 tcr_tr_col,
                                 vdj_cdr3_col,
                                 tcr_cdr3_col,
                                 nthread) {
  matches <- mapply(stringdist::stringdistmatrix,
                    split(vdjdb[,vdj_cdr3_col,drop=T], vdjdb[,vdj_tr_col,drop=T])[sort(unique(vdjdb[,vdj_tr_col,drop=T]), decreasing = sort_desc)],
                    split(tcrs[,tcr_cdr3_col,drop=T], tcrs[,tcr_tr_col,drop=T])[sort(unique(tcrs[,tcr_tr_col,drop=T]))],
                    method = "lv",
                    useNames = "strings",
                    SIMPLIFY = F,
                    nthread = nthread)
  names(matches) <- paste(sort(unique(vdjdb[,vdj_tr_col,drop=T]), decreasing = sort_desc), sort(unique(tcrs[,tcr_tr_col,drop=T])), sep = "_")
  matches <- sapply(matches, reshape2::melt, simplify = F, c(vdj_cdr3_col, tcr_cdr3_col), value.name = "lv")
  matches <- do.call(rbind, matches)
  matches[,vdj_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 1)
  matches[,tcr_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 2)
  rownames(matches) <- NULL
  return(matches)
}
