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
#' @param tcrs data frame of TCR data from sequencing, tcr_cdr3_col and tcr_tr_col must be present; e.g. cl_long from clonotype prep
#' @param vdjdb data frame of TCR data from vdjdb; vdj_cdr3_col and vdj_tr_col must be present;
#' if not provided system.file("extdata", "vdjdb.tsv.rds", package = "igsc") is used (table downloaded 06/12/2021)
#' @param vdj_tr_col column which codes the TCR chain in the vdjdb reference data frame (should contain TRA and/or TRB only)
#' @param tcr_tr_col column which codes the TCR chain in the tcrs data frame (should contain TRA and/or TRB only)
#' @param vdj_cdr3_col column which codes the CDR3s in the vdjdb reference data frame
#' @param tcr_cdr3_col column which codes the CDR3s in the tcrs data frame
#' @param TRAxTRB logical, apart from comparing TRA vs TRA and TRB vs TRB between tcrs and vdjdb also check for TRA vs TRB?
#' @param max_lvdist maximum levensthein distance between CDR3s for return data frame (also filtering before PAM30 similarity calculation)
#' @param PAM30_similarity logical whether to calculate the similarity (pairwise alignment) between CDR3s based on PAM30 substitution matrix
#' @param mapply_fun mapply function for PAM30_similarity, suggested are mapply, pbapply::pbmapply, parallel::mcmapply
#' @param ... arguments passed to mapply_fun and lapply_fun, most relevant: mc.cores (e.g parallel::detectCores())
#' @param lapply_fun lapply function for levensthein distance calculation, suggested are lapply, pblapply::pblapply, parallel::mclapply
#'
#' @return data frame of matched CDR3s from tcrs and vdjdb, column lv indicates the levensthein
#' between CDR3 from vdjdb and tcrs
#' @export
#'
#' @examples
#' \dontrun{
#' # tcrs data frame with sequencing data from TCRs has to be created (cl_long is appropriate)
#' matches <- vdjdb_hits(vdjdb, tcrs, mapply_fun = parallel::mcmapply, mc.cores = parallel::detectCores(), nthread = parallel::detectCores())
#' matches_append <- matches %>% dplyr::left_join(tcrs) %>% dplyr::left_join(vdjdb)
#'
#' df <-
#' df %>%
#' dplyr::distinct(CDR3, CDR3_aa_cr, lv, Gene, chain, cl_name, patient, Epitope, Epitope.gene, Epitope.species, Score) %>%
#' tidyr::drop_na() %>%
#' dplyr::filter(cl_name %in% c("Madalyn", "Mikaylla", "Morgan")) %>%
#' dplyr::filter(lv < 3) %>%
#' dplyr::filter(Score > 0)
#'
#' # use jitter for many hits (when beeswarms become overlapping)
#' ggplot(df, aes(x = Score, y = lv, color = Epitope.species)) +
#'  ggbeeswarm::geom_beeswarm(groupOnX = F, cex = 4, size = 2.5) +
#'  #geom_jitter(width = 0.15, height = 0.15, size = 3) +
#'  theme_bw() +
#'  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), strip.background = element_rect(fill = "white"), text = element_text(family = "Courier")) +
#'  scale_x_continuous(breaks = fcexpr::int_breaks(n = 3), limits = c(0.8,3.2)) +
#'  scale_y_reverse(breaks = fcexpr::int_breaks(n = 3)) +
#'  guides(fill = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) +
#'  facet_grid(rows = vars(chain), cols = vars(patient, cl_name))
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
                       lapply_fun = lapply,
                       ...) {


  if (vdj_tr_col == tcr_tr_col) {
    stop("Please avoid same names for tcr_tr_col and vdj_tr_col.")
  }
  if (vdj_cdr3_col == tcr_cdr3_col) {
    stop("Please avoid same names for tcr_cdr3_col and vdj_cdr3_col")
  }

  if (missing(vdjdb)) {
    vdjdb <- readRDS(system.file("extdata", "vdjdb.tsv.rds", package = "igsc"))
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
  lapply_fun <- match.fun(lapply_fun)
  vdjdb <- dplyr::distinct(vdjdb, !!rlang::sym(vdj_cdr3_col), !!rlang::sym(vdj_tr_col))
  tcrs <- dplyr::distinct(tcrs, !!rlang::sym(tcr_tr_col), !!rlang::sym(tcr_cdr3_col))

  # what if only TRA or TRB is provided in one of the two? - duplicate results due to mapply?
  matches <- do.call(rbind, lapply(unique(c(F, TRAxTRB)),
                                   .vdjdb_tcrs_match_fun,
                                   vdjdb = vdjdb,
                                   tcrs = tcrs,
                                   vdj_tr_col = vdj_tr_col,
                                   tcr_tr_col = tcr_tr_col,
                                   vdj_cdr3_col = vdj_cdr3_col,
                                   tcr_cdr3_col = tcr_cdr3_col,
                                   max_lvdist = max_lvdist,
                                   lapply_fun = lapply_fun,
                                   ...))

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

.vdjdb_tcrs_match_fun <- function(sort_desc,
                                 vdjdb,
                                 tcrs,
                                 vdj_tr_col,
                                 tcr_tr_col,
                                 vdj_cdr3_col,
                                 tcr_cdr3_col,
                                 max_lvdist,
                                 lapply_fun,
                                 ...) {

  tcrs_split <- split(tcrs, ceiling(1:nrow(tcrs)/100))

  matches <- lapply_fun(tcrs_split, function(x) {
    #print(which(mapply(identical, tcrs_split, list(x))))
    matches <- mapply(stringdist::stringdistmatrix,
                      split(vdjdb[,vdj_cdr3_col,drop=T], vdjdb[,vdj_tr_col,drop=T])[sort(unique(vdjdb[,vdj_tr_col,drop=T]), decreasing = sort_desc)],
                      split(x[,tcr_cdr3_col,drop=T], x[,tcr_tr_col,drop=T])[sort(unique(x[,tcr_tr_col,drop=T]))],
                      method = "lv",
                      useNames = "strings",
                      SIMPLIFY = F,
                      nthread = 1)
    matches <- sapply(matches, function(y) y[which(apply(y, 1, min) <= max_lvdist), which(apply(y, 2, min) <= max_lvdist), drop = F], simplify = F)
    matches <- sapply(matches, reshape2::melt, simplify = F, c(vdj_cdr3_col, tcr_cdr3_col), value.name = "lv")
    names(matches) <- paste(sort(unique(vdjdb[,vdj_tr_col,drop=T]), decreasing = sort_desc), sort(unique(x[,tcr_tr_col,drop=T])), sep = "_")
    matches <- do.call(rbind, matches)
    matches <- matches[which(matches$lv <= max_lvdist),]
    if (nrow(matches) == 0) {
      return(NULL)
    }
    matches[,vdj_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 1)
    matches[,tcr_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 2)
    rownames(matches) <- NULL
    return(matches)
  }, ...)

  matches <- do.call(rbind, matches)
  rownames(matches) <- NULL

  return(matches)
}

