#' Compare CDR3s from sequencing to known CDR3s in the vdjdb
#'
#' Matching CDR3s from sequencing with CDR3s of known antigen specificity may help to narrow potential antigen specificities
#' of one's sequenced TCRs. Export a tsv of all or selected entries here https://vdjdb.cdr3.net/search. Provide that tsv as
#' data frame (vdjdb). From vdjdb and tcrs distinct rows of tr_cols and cdr3_col are compared by levensthein distance
#' and if selected (PAM30_similarity) by amino acid similarity. Returned data frame may be used to join original columns
#' from vdjdb and tcrs followed by plotting with ggplot2 (see example). Original column names (vdj_tr_col, tcr_tr_col and
#' vdj_cdr3_col, tcr_cdr3_col) are kept in the returned data frame to allow easy subsequent joining of vdjdb and tcrs. It requires
#' though that these columns have different names in vdjdb and tcrs.
#'
#' @param tcrs data frame of TCR data from sequencing, first column indicating the chain TRA or TRB, second the aa sequence of CDR3
#' @param vdjdb data frame of TCR data from vdjdb, first column indicating the chain TRA or TRB, second the aa sequence of CDR3;
#' if not provided base::system.file("extdata", "vdjdb.tsv.tar.gz", package = "igsc") is used
#' @param vdj_tr_col column (name) which codes the TCR chain (TRA or TRB) in the vdjdb reference data frame
#' @param tcr_tr_col column (name) which codes the TCR chain (TRA or TRB) in the tcrs data frame
#' @param vdj_cdr3_col column (name) which codes the CDR3s in the vdjdb reference data frame
#' @param tcr_cdr3_col column (name) which codes the CDR3s in the tcrs data frame
#' @param TRAxTRB logical, apart from comparing TRA vs TRA and TRB vs TRB between tcrs and vdjdb also check for TRA vs TRB?
#' @param max_lvdist maximum levensthein distance between CDR3s for return data frame (also filtering before PAM30 similarity calculation)
#' @param PAM30_similarity logical whether to calculate the similarity between CDR3s based on PAM30 substitution matrix
#' @param mapply_fun mapply function for PAM30_similarity, suggested are mapply, pbapply::pbmapply, parallel::mcmapply
#' @param nthread number of threads (cores) to use for levensthein distance calculation
#' @param ... arguments passed to mapply_fun, most relevant: mc.cores (e.g parallel::detectCores())
#'
#' @return data frame of matched CDR3s from tcrs and vdjdb, column lv indicates the levensthein
#' between CDR3 from vdjdb and tcrs
#' @export
#'
#' @examples
#' \dontrun{
#' vdjdb <- dplyr::filter(read.csv("vdjdb.tsv", sep = "\t", header = T) , Species == "HomoSapiens")
#' matches <- vdjdb_hits(vdjdb, tcrs, mapply_fun = parallel::mcmapply, mc.cores = parallel::detectCores(), nthread = parallel::detectCores())
#' matches <- matches %>% dplyr::left_join(tcrs) %>% dplyr::left_join(vdjdb) %>% dplyr::distinct(complex.id, CDR3, CDR3_aa_cr, lv, Gene, chain, clonotype_id_cr, sample, Epitope, Epitope.gene, Epitope.species, Score)
#'
#' ggplot(matches, aes(x = chain, y = lv, fill = Epitope.species, color = Gene)) +
#' geom_jitter(width = 0.15, height = 0.15, shape = 21) +
#' scale_y_continuous(breaks = fcexpr::int_breaks(n = 3)) +
#' theme_bw() +
#' theme(panel.grid.major.x = element_blank()) +
#' scale_color_manual(values = c("#000000", "#999999")) +
#' guides(fill = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) +
#' facet_grid(rows = vars(Score), cols = vars(clonotype_id_cr))
#' }
vdjdb_hits <- function(vdjdb,
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

  if (missing(vdjdb)) {
    utils::untar(base::system.file("extdata", "vdjdb.tsv.tar.gz", package = "igsc"), exdir = tempdir())
    vdjdb <- read.csv(file.path(tempdir(), "vdjdb.tsv"), sep = "\t", header = T)
  }

  if (missing(tcrs)) {
    stop("Please provide a tcrs data frame.")
  }

  # check if columns exist

  if (!identical(sort(unique(vdjdb[,vdj_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(vdj_tr_col, " column of vdjdb should only contain TRA and/or TRB."))
  }
  if (!identical(sort(unique(tcrs[,tcr_tr_col,drop=T])), c("TRA", "TRB"))) {
    stop(paste0(tcr_tr_col, " column of tcrs should only contain TRA and/or TRB."))
  }

  mapply_fun <- match.fun(mapply_fun)
  vdjdb <- dplyr::distinct(vdjdb, !!sym(vdj_cdr3_col), !!sym(vdj_tr_col))
  tcrs <- dplyr::distinct(tcrs, !!sym(tcr_tr_col), !!sym(tcr_cdr3_col))

  # what if only TRA or TRB is provided in one of the two? - duplicate results due to mapply?
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

  return(matches)

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

