#' Align IMGTs reference sequences to your TCR sequence
#'
#' With Biostrings::pairwiseAlignment() the borders of the V(D)J-segments in your TCR sequence can be determined.
#' The sequence within these border may then be used for further steps. In case of mouse, IMGT gene segments from Mus musculus are preferred.
#'
#' @param chain character which TCR chain to align, currently only TRA and TRB
#' @param TCR a named vector; which of your TCRs to align, the name indicates which column name of cl_long to use, the value indicates what to look for in that column
#' @param imgt_ref the prepared data.frame of IMGT references
#' @param sequence_col name of the column to pull sequences from
#' @param cl_long clonotype data frame long format
#' @param cl_wide clonotype data frame wide format
#' @param C_allele optional name of the constant allele to use for the alignment, must be the entry in the "Allele"-column of the imgt_ref data frame
#' @param type type of alignment
#' @param ...  arguments to pwalign_multi
#'
#' @return a list
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
align_imgt_ref_to_TCR_seq <- function(chain,
                                      TCR,
                                      cl_long,
                                      #cl_wide,
                                      imgt_ref,
                                      sequence_col = "consensus_seq",
                                      C_allele = NULL,
                                      type = "local",
                                      ...) {

  if (!requireNamespace("BiocManager", quietly = T)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("Biostrings", quietly = T)){
    BiocManager::install("Biostrings")
  }
  if (!requireNamespace("DECIPHER", quietly = T)){
    BiocManager::install("DECIPHER")
  }

  chain <- match.arg(chain, c("TRA", "TRB"))
  cl_long <- as.data.frame(cl_long)
  #cl_wide <- as.data.frame(cl_wide)

  V_imgt.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "V_imgt"])
  if (length(V_imgt.name) > 1) {
    message("More than one V allele: ", paste(V_imgt.name, collapse = ", "), ".")
    message("Splitting output by those. Double check results, please.")
  }

  J_imgt.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "J_imgt"])
  if (length(J_imgt.name) > 1) {
    message("More than one J allele: ", paste(J_imgt.name, collapse = ", "), ".")
    message("Splitting output by those. Double check results, please.")
  }

  cl_long_sub <-
    cl_long %>%
    dplyr::filter(chain == !!chain) %>% # https://stackoverflow.com/questions/34219912/how-to-use-a-variable-in-dplyrfilter, https://stackoverflow.com/questions/40169949/filter-dataframe-using-global-variable-with-the-same-name-as-column-name
    dplyr::filter(!!rlang::sym(names(TCR)) == TCR) %>%
    dplyr::select(V_imgt, J_imgt) %>%
    dplyr::distinct(V_imgt, J_imgt)

  pairs <- character(nrow(cl_long_sub))
  for (i in 1:nrow(cl_long_sub)) {
    pairs[i] <- paste(cl_long_sub$V_imgt[i], cl_long_sub$J_imgt[i], sep = ",")
  }

  #pairs <- paste(strsplit(unique(cl_wide[which(cl_wide[,names(TCR)] == TCR),paste0("V_imgt_", chain)]), ",")[[1]], strsplit(unique(cl_wide[which(cl_wide[,names(TCR)] == TCR),paste0("J_imgt_", chain)]), ",")[[1]], sep = ",")

  out <- lapply(pairs, function(x) {
    v <- strsplit(x, ",")[[1]][1]
    j <- strsplit(x, ",")[[1]][2]
    raw.cs <- cl_long[Reduce(intersect, list(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR), which(cl_long[,"V_imgt"] == v), which(cl_long[,"J_imgt"] == j))), sequence_col]
    names(raw.cs) <- paste0("seq_", seq_along(raw.cs))
    if (length(raw.cs) > 1) {
      raw.cs <- collapse_duplicate_sequences(seq_set = raw.cs)

      if (length(raw.cs) > 1) {
        consensus_seq <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(raw.cs), verbose = F)
        consensus_seq <- DECIPHER::ConsensusSequence(consensus_seq, includeTerminalGaps = F, minInformation = 0.8, threshold = 0.1)
        consensus_seq <- stats::setNames(as.character(consensus_seq), "consensus")
      } else {
        consensus_seq <- stats::setNames(raw.cs, "consensus")
      }


      raw.cs <- c(raw.cs, consensus_seq)
      al_df <- check_ref_seq_for_matches(seq_set = DECIPHER::AlignSeqs(Biostrings::DNAStringSet(raw.cs), verbose = F),
                                         ref_seq_name = "consensus",
                                         pos_col = "position",
                                         seq_col = "seq",
                                         name_col = "seq.name")

      p1 <- algnmt_plot(al_df) + ggplot2::ggtitle(paste(paste0(unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), names(TCR)]), "_", v, "_", j), collapse = ", "))
      cs <- consecutive_distinct_seq(consensus_seq, seq_type = "NT")
      p1 <- p1 + ggplot2::geom_vline(xintercept = cs[["limits"]], linetype = "dashed")
      cs <- cs[["seq"]]
      names(cs) <- paste0(TCR, "_consensus")

    } else if (length(raw.cs) == 1)  {
      cs <- raw.cs
      names(cs) <- TCR
      p1 <- NULL
    } else {
      stop("No sequence found. Does the TCR exist in the respective column? Or is it missing the alpha or beta chain?")
    }
    return(list(p1, cs))
  })

  p1 <- sapply(out, "[", 1)
  names(p1) <- pairs
  cs <- unlist(sapply(out, "[", 2))
  names(cs) <- pairs

  # imgt ref seqs
  imgt_v_allele_seq <- unlist(lapply(V_imgt.name, function(x) {
    imgt_v_allele_seq <- imgt_ref[which(imgt_ref$Allele %in% x),] %>% dplyr::distinct(seq.nt, meta, Allele)
    if (nrow(imgt_v_allele_seq) > 1) {
      if (sum(grepl("Mus musculus", imgt_v_allele_seq$meta)) == 1) {
        imgt_v_allele_seq <- imgt_v_allele_seq[which(grepl("Mus musculus", imgt_v_allele_seq$meta)),]
      } else {
        imgt_v_allele_seq <- imgt_v_allele_seq[sample(1:nrow(imgt_v_allele_seq), 1),]
      }
    }
    return(imgt_v_allele_seq[,"seq.nt",drop=T])
  }))
  names(imgt_v_allele_seq) <- V_imgt.name

  imgt_j_allele_seq <- unlist(lapply(J_imgt.name, function(x) {
    imgt_j_allele_seq <- imgt_ref[which(imgt_ref$Allele %in% x),] %>% dplyr::distinct(seq.nt, meta, Allele)
    if (nrow(imgt_j_allele_seq) > 1) {
      if (sum(grepl("Mus musculus", imgt_j_allele_seq$meta)) == 1) {
        imgt_j_allele_seq <- imgt_j_allele_seq[which(grepl("Mus musculus", imgt_j_allele_seq$meta)),]
      } else {
        imgt_j_allele_seq <- imgt_j_allele_seq[sample(1:nrow(imgt_j_allele_seq), 1),]
      }
    }
    return(imgt_j_allele_seq[,"seq.nt",drop=T])
  }))
  names(imgt_j_allele_seq) <- J_imgt.name

  if (is.null(C_allele)) {
    C.allele.seq <- imgt_ref[intersect(which(grepl(chain, imgt_ref$Allele)), which(grepl("C", imgt_ref$Allele)))[1], "seq.nt"]
    names(C.allele.seq) <- imgt_ref[intersect(which(grepl(chain, imgt_ref$Allele)), which(grepl("C", imgt_ref$Allele)))[1], "Allele"]
    if (length(C.allele.seq) > 1) {
      C.allele.seq <- C.allele.seq[sample(1:length(C.allele.seq), 1)]
    }
  } else {
    C.allele.seq <- imgt_ref[which(imgt_ref$Allele == C_allele), "seq.nt"]
    names(C.allele.seq) <- C_allele
  }

  # what about one V allele but 2 J allele or vice versa?
  p2 <- lapply(pairs, function (x) {
    v <- strsplit(x, ",")[[1]][1]
    j <- strsplit(x, ",")[[1]][2]

    pwalign_multi(subject = Biostrings::DNAStringSet(cs[x]),
                                           patterns = Biostrings::DNAStringSet(c(imgt_v_allele_seq[v],imgt_j_allele_seq[j], C.allele.seq)),
                                           type = type,
                                           algnmt_plot_args = list(add_length_suffix = T, subject.lim.lines = T, subject_name = names(cs[x])))
  })
  names(p2) <- pairs

  p3 <- lapply(pairs, function (x) {
    v <- strsplit(x, ",")[[1]][1]
    j <- strsplit(x, ",")[[1]][2]

    pwalign_multi(subject = Biostrings::DNAStringSet(cs[x]),
                                           patterns = Biostrings::DNAStringSet(c(imgt_v_allele_seq[v],imgt_j_allele_seq[j])),
                                           type = type,
                                           algnmt_plot_args = list(add_length_suffix = T, subject.lim.lines = T, subject_name = names(cs[x])))
  })
  names(p3) <- pairs

  TCR.seq <- unlist(lapply(pairs, function (x) stringr::str_sub(cs[x], start = p3[[x]][["min.max.subject.position"]][1], end = p3[[x]][["min.max.subject.position"]][2])))
  names(TCR.seq) <- paste0(TCR, "_", chain, "_", pairs)

  return(list("consensus.alingments" = p1, "consensus.seqs.uncut" = cs, "VJC.alignments" = p2, "VJ.alignments" = p3, "V.to.J.TCR.seqs" = TCR.seq))
}
