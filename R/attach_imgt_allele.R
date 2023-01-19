#' Find best matching IMGT reference alleles for TCR V and J gene segments
#'
#' Entries from seq are made unique and are matched against reference alleles in imgt_ref.
#' If names are provided to seq (or respective V and J columns in case of a data frame) matching will be
#' quicker as possibles alleles are initially narrowed down by string matching.
#'
#' @param seq data frame with clonotype data in long format, e.g. cl_long preferentially from igsc::read_cellranger_outs or
#' a named vector of consensus sequences, named by V and J segments, separated by "___"; e.g. stats::setNames(object = cl_long$consensus_seq_cr, nm = paste0(cl_long$V_cr, '___', cl_long$J_cr))
#' @param imgt_ref IMGT reference data frame created with igsc::imgt_tcr_segment_prep or a named vector of sequences and respective allele names;
#' e.g. stats::setNames(imgt_ref$seq.nt, imgt_ref$Allele)
#' @param pick.by match disambiguate alleles from IMGT by best match (alignment) or just randomly (random);
#' random is intended to speed up the process for testing or when exact IMGT alleles are not necessary
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return a data frame of unique seq entries in column
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' imgt_ref <- readRDS(system.file("extdata", "IMGT_ref/human/hs.rds", package = "igsc"))
#' ata <- attach_imgt_alleles(seq = cl_long, imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)
#' cl_long <-
#' dplyr::left_join(cl_long, ata, by = c("consensus_seq_cr" = "seq")) %>%
#' tidyr::separate(VJ_IMGT, into = c("V_imgt", "J_imgt"), sep = "___")
#'
#' ## or pass vector, named or unnamed
#' ata <- attach_imgt_alleles(seq = cl_long$consensus_seq_cr, imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)
#'
#' ata <- attach_imgt_alleles(seq = stats::setNames(object = cl_long$consensus_seq_cr, nm = paste0(cl_long$V_cr, '___', cl_long$J_cr)),
#' imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)
#' }
attach_imgt_alleles <- function(seq,
                                imgt_ref,
                                pick.by = "alignment",
                                lapply_fun = lapply,
                                seq_col = "consensus_seq_cr",
                                V_col = "V_cr",
                                J_col = "J_cr",
                                ...) {

  if (missing(seq)) {
    stop("Provide a data frame for seq (e.g. cl_long), preferentially prepared with igsc::read_cellranger_out or a named vector, e.g.
         stats::setNames(object = seq$consensus_seq_cr, nm = paste0(seq$V_cr, '___', seq$J_cr))")
  }

  if (missing(imgt_ref)) {
    stop("Provide an imgt_ref dataframe, e.g. readRDS(system.file('extdata', 'IMGT_ref/human/hs.rds', package = 'igsc') for human or
         readRDS(system.file('extdata', 'IMGT_ref/mouse/mm.rds', package = 'igsc')) for mouse; or a named vector of reference sequences.")
  }

  if (is.data.frame(seq)) {
    if (any(!c(seq_col, V_col, J_col) %in% names(seq))) {
      stop("seq data frame has to have ", seq_col , " ", V_col, " ", J_col, " columns.")
    }
    un <- which(!duplicated(seq[,seq_col,drop=T]))
    seq <- stats::setNames(object = seq[un,seq_col,drop=T], nm = paste0(seq[un,V_col,drop=T], "___", seq[un,J_col,drop=T]))
  } else if (is.character(seq)) {
    if (!is.null(names(seq))) {
      if(!all(lengths(strsplit(names(seq)[which(!is.na(names(seq)))], "___")) == 2)) {
        stop("Names of seq have to look like this: 'TRBV11-3___TRBJ2-5' or 'TRAV___TRAJ' or NA; V- and J- segments have to be split by '___'.")
      }
    }
    seq <- seq[which(!duplicated(seq))]
  } else {
    stop("seq has to be data frame or named character vector.")
  }


  if (is.data.frame(imgt_ref)) {
    if (any(!c("seq.nt", "Allele") %in% names(imgt_ref))) {
      stop("imgt_ref data frame has to have a 'seq.nt' and 'Allele' column.")
    }
    ref_seq <- stats::setNames(imgt_ref$seq.nt, imgt_ref$Allele)
  } else if (is.character(imgt_ref)) {
    if (is.null(names(imgt_ref))) {
      stop("vector of sequences has to have names.")
    }
    ref_seq <- imgt_ref
  } else {
    stop("imgt_ref has to be data frame or named character vector.")
  }

  pick.by <- match.arg(pick.by, c("alignment", "random"))
  lapply_fun <- match.fun(lapply_fun)

  if (!is.null(names(seq))) {
    v_m <- .match_imgt_str(seq_names = unique(sapply(strsplit(names(seq), "___"), "[", 1)), ref_seq_names = names(ref_seq))
    j_m <- .match_imgt_str(seq_names = unique(sapply(strsplit(names(seq), "___"), "[", 2)), ref_seq_names = names(ref_seq))
  } else {
    v_m <- NULL
    j_m <- NULL
  }

  browser()
  names(seq) <- unlist(lapply_fun(seq_along(seq), function(x) {
    print(x)
    if (!is.null(names(seq[x])) && !is.na(names(seq[x]))) {
      TRV_name <- strsplit(names(seq[x]), "___")[[1]][1]
      if (grepl("TR[[:alpha:]]{1}V", TRV_name)) {
        v_sn <- v_m[[TRV_name]]
      } else {
        v_sn <- ""
      }
      TRJ_name <- strsplit(names(seq[x]), "___")[[1]][2]
      if (grepl("TR[[:alpha:]]{1}J", TRJ_name)) {
        j_sn <- j_m[[TRJ_name]]
      } else {
        j_sn <- ""
      }
    } else {
      v_sn <- grep("V", names(ref_seq), value = T, ignore.case = T)
      j_sn <- grep("J", names(ref_seq), value = T, ignore.case = T)
    }


    vj <- lapply(list(v_sn, j_sn), function (y) {
      if (length(y) == 1) {
        return(y)
      } else {
        if (pick.by == "alignment") {
          max.ind <- which.max(Biostrings::pairwiseAlignment(subject = Biostrings::DNAString(seq[x]), pattern = Biostrings::DNAStringSet(ref_seq[y]), type = "local", scoreOnly = T))
          max.ind <- max.ind[sample(seq_along(max.ind), 1)] # if multiple matches with equal score, pick one randomly
          return(y[max.ind])
        } else if (pick.by == "random") {
          return(y[sample(seq_along(y),1)])
        }
      }
    })
    vj <- unlist(vj)
    if (length(vj) != 2) {
      stop("vj should have length of 2.")
    }
    vj <- paste0(vj[1], "___", vj[2])
    return(vj)
  }, ...))
  seq <- utils::stack(seq)
  names(seq) <- c("seq", "VJ_IMGT")
  return(seq)
}


.match_imgt_str <- function(seq_names, ref_seq_names) {
  sapply(seq_names, function(x) {
    ref_al <- grep(stringr::str_extract(x, "^[:alpha:]{1,}[:digit:]{1,}"), ref_seq_names, ignore.case = T, value = T)
    dists <- utils::adist(x, ref_al, ignore.case = T)
    min_val <- min(dists)
    return(unique(ref_al[which(dists == min_val)]))
  })
}
