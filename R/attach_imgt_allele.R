#' Match best matching IMGT reference alleles (TCR V and J) and attache respective columns to cl_long
#'
#' @param cl_long clonotype data in long format, preferentially from igsc::read_cellranger_outs
#' @param imgt_ref imgt reference data frame
#' @param pick.by match disambiguate alleles from IMGT by best match (alignment) or just randomly (random);
#' random is offered to speed up the process in case of many alleles
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
attach_imgt_alleles <- function(cl_long,
                                imgt_ref,
                                pick.by = "alignment",
                                lapply_fun = lapply,
                                ...) {

  if (missing(imgt_ref)) {
    stop("Provide an imgt_ref dataframe, e.g. readRDS(system.file('extdata', 'IMGT_ref/human/hs.rds', package = 'igsc'))")
  }
  if (missing(cl_long)) {
    stop("Provide a cl_long dataframe, preferentially prepared from igsc::read_cellranger_out")
  }

  if (length(intersect(c("consensus_seq_cr", "clonotype_id_cr", "sample", "V_cr", "J_cr"), names(cl_long))) < 5) {
    stop("cl_long must have the following columns: consensus_seq_cr, clonotype_id_cr, sample, V_cr, J_cr.")
  }

  pick.by <- match.arg(pick.by, c("alignment", "random"))
  lapply_fun <- match.fun(lapply_fun)

  uniques <- lapply(c("V", "J"), function(i) {
    matches <- match_imgt(unique(cl_long[,paste0(i, "_cr"),drop=T]))
    cl_long[[paste0(i, "_imgt")]] <- matches[cl_long[,paste0(i, "_cr"),drop=T]]

    unique <- dplyr::distinct(cl_long, consensus_seq_cr, clonotype_id_cr, sample, !!sym(paste0(i, "_imgt")))

    unique[[paste0(i, "_imgt")]] <- unlist(lapply_fun(1:nrow(unique), function(x) {
      als <- unique(unlist(unique[x,paste0(i, "_imgt")]))
      if (length(als) == 1) {
        return(unname(als))
      } else {
        if (pick.by == "alignment") {
          seq <- unique[x,"consensus_seq_cr"]
          refs <- imgt_ref[which(imgt_ref$Allele %in% als),]
          max.ind <- which.max(Biostrings::pairwiseAlignment(subject = Biostrings::DNAString(seq), pattern = Biostrings::DNAStringSet(refs[,"seq.nt"]), type = "local", scoreOnly = T))
          max.ind <- sample(1:length(max.ind), 1) # if multiple matches with equal score, pick one randomly
          return(refs[max.ind,"Allele"])
        } else if (pick.by == "random") {
          refs <- imgt_ref[which(imgt_ref$Allele %in% als),]
          return(refs[sample(seq(1,nrow(refs)), 1),"Allele"])
        }
      }
    }, ...))
    return(unique)
  })

  cl_long <-
    cl_long %>%
    dplyr::left_join(uniques[[1]], by = c("consensus_seq_cr", "clonotype_id_cr", "sample")) %>%
    dplyr::left_join(uniques[[2]], by = c("consensus_seq_cr", "clonotype_id_cr", "sample"))

  return(cl_long)
}


match_imgt <- function(alleles) {
  sapply(alleles, function(x) {
    ref_al <- imgt_ref$Allele[which(grepl(stringr::str_extract(x, "^[:alpha:]{1,}[:digit:]{1,}"), imgt_ref$Allele))]
    dists <- utils::adist(x, ref_al, ignore.case = T)
    min_val <- min(dists)
    return(unique(ref_al[which(dists == min_val)]))
  })
}
