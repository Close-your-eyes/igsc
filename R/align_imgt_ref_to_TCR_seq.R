#' Align IMGTs reference sequences to your TCR sequence
#'
#' With Biostrings::pairwiseAlignment() the borders of the V(D)J-segments in your TCR sequence can be determined.
#' The sequence within these border may then be used for further steps.
#'
#' @param chain character which TCR chain to align, currently only TRA and TRB
#' @param TCR a named vector; which of your TCRs to align, the name indicates which column name of cl_long to use, the value indicates what to look for in that column
#' @param cl_long the prepared clonotype data frame in long format
#' @param imgt_ref the prepared data.frame of IMGT references
#'
#' @return a list
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
align_imgt_ref_to_TCR_seq <- function(chain, TCR, cl_long, imgt_ref) {
  chain <- match.arg(chain, c("TRA", "TRB"))

  raw.cs <- cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "cell.ranger.consensus.seq"]

  if (length(raw.cs) > 1) {
    p1 <- MultipleSequenceAlignmentDECIPHER(input.set = raw.cs) + ggplot2::ggtitle(paste(unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), names(TCR)]), collapse = ", "))
    cs <- DECIPHER::ConsensusSequence(DECIPHER::AlignSeqs(Biostrings::DNAStringSet(raw.cs), verbose = F))
    cs <- consecutive.disambiguate.consensus.seq.from.decipher.consensus(cs)
    names(cs) <- paste0(TCR, "_consensus")
  } else {
    cs <- raw.cs
    names(cs) <- TCR
    p1 <- NULL
  }

  ## pull out reference sequences
  V.allele.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "V.allele"])
  if (length(V.allele.name) > 1) {
    print(paste0("More than one V allele: ", paste(V.allele.name, collapse = ", "), "."))
    V.allele.name <- V.allele.name[sample(1:length(V.allele.name), 1)]
    print(paste0("Picking: ", V.allele.name, "."))
  }
  V.allele.seq <- imgt_ref[which(imgt_ref$Allele == V.allele.name), "seq.nt"]
  names(V.allele.seq) <- V.allele.name

  J.allele.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "J.allele"])
  if (length(J.allele.name) > 1) {
    print(paste0("More than one J allele: ", paste(J.allele.name, collapse = ", "), "."))
    J.allele.name <- J.allele.name[sample(1:length(J.allele.name), 1)]
    print(paste0("Picking: ", J.allele.name, "."))
  }
  J.allele.seq <- imgt_ref[which(imgt_ref$Allele == J.allele.name), "seq.nt"]
  names(J.allele.seq) <- J.allele.name

  if (chain == "TRA") {
    C.allele.seq <- imgt_ref[which(grepl("TRAC", imgt_ref$Allele)), "seq.nt"]
    names(C.allele.seq) <- "TRAC*01"
  } else if (chain == "TRB") {
    C.allele.seq <- imgt_ref[which(grepl("TRBC1\\*01", imgt_ref$Allele)), "seq.nt"]
    names(C.allele.seq) <- "TRBC1*01"
  }

  p2 <- MultiplePairwiseAlignmentsToOneSubject(subject = Biostrings::DNAStringSet(cs), patterns = Biostrings::DNAStringSet(c(V.allele.seq, J.allele.seq, C.allele.seq)), type = "local", pattern.positions.size = 4, print.subject.min.max = T)
  p3 <- MultiplePairwiseAlignmentsToOneSubject(subject = Biostrings::DNAStringSet(cs), patterns = Biostrings::DNAStringSet(c(V.allele.seq, J.allele.seq)), type = "local", pattern.positions.size = 4, print.subject.min.max = T)

  TCR.seq <- stringr::str_sub(cs, start = p3[["min.max.subject.position"]][1], end = p3[["min.max.subject.position"]][2])
  names(TCR.seq) <- paste0(TCR, "_", chain)
  return(list("consensus.alingment" = p1, "consensus.seq.uncut" = cs, "VJC.alignment" = p2, "VJ.alignment" = p3, "V.to.J.TCR.seq" = TCR.seq))
}
