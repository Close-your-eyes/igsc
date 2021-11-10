#' Align IMGTs reference sequences to your TCR sequence
#'
#' With Biostrings::pairwiseAlignment() the borders of the V(D)J-segments in your TCR sequence can be determined.
#' The sequence within these border may then be used for further steps.
#'
#' @param chain character which TCR chain to align, currently only TRA and TRB
#' @param TCR a named vector; which of your TCRs to align, the name indicates which column name of cl_long to use, the value indicates what to look for in that column
#' @param cl_long the prepared clonotype data frame in long format
#' @param imgt_ref the prepared data.frame of IMGT references
#' @param sequence_col name of the column to pull sequences from, defaults to "cell.ranger.consensus.seq"
#' @param C_allele optional name of the constant allele to use for the alignment, must be the entry in the "Allele"-column of the imgt_ref data frame
#'
#' @return a list
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
align_imgt_ref_to_TCR_seq <- function(chain, TCR, cl_long, imgt_ref, sequence_col = "cell.ranger.consensus.seq", C_allele) {
  chain <- match.arg(chain, c("TRA", "TRB"))

  V.allele.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "V.allele"])
  if (length(V.allele.name) > 1) {
    print(paste0("More than one V allele: ", paste(V.allele.name, collapse = ", "), "."))
    print("Splitting output by those. Double check results, please.")
  }

  J.allele.name <- unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), "J.allele"])
  if (length(J.allele.name) > 1) {
    print(paste0("More than one J allele: ", paste(J.allele.name, collapse = ", "), "."))
    print("Splitting output by those. Double check results, please.")
  }

  out <- lapply(V.allele.name, function(x) {
    raw.cs <- cl_long[Reduce(intersect, list(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR), which(cl_long[,"V.allele"] == x))), sequence_col]
    if (length(raw.cs) > 1) {
      J_al <- unique(cl_long[Reduce(intersect, list(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR), which(cl_long[,"V.allele"] == x))), "J.allele"])
      p1 <- MultipleSequenceAlignmentDECIPHER(input.set = raw.cs) + ggplot2::ggtitle(paste(paste0(unique(cl_long[intersect(which(cl_long$chain == chain), which(cl_long[,names(TCR)] == TCR)), names(TCR)]), "_", x, "_", J_al), collapse = ", "))
      cs <- DECIPHER::ConsensusSequence(DECIPHER::AlignSeqs(Biostrings::DNAStringSet(raw.cs), verbose = F))
      cs <- consecutive.disambiguate.consensus.seq.from.decipher.consensus(cs)
      names(cs) <- paste0(TCR, "_consensus")
    } else if (length(raw.cs) == 1)  {
      cs <- raw.cs
      names(cs) <- TCR
      p1 <- NULL
    } else {
      stop("No sequence found. Does the TCR exist in the respective column?")
    }
    return(list(p1, cs))
  })

  p1 <- sapply(out, "[", 1)
  names(p1) <- paste0(V.allele.name,"_",J.allele.name)
  cs <- unlist(sapply(out, "[", 2))
  names(cs) <- paste0(V.allele.name,"_",J.allele.name)

  imgt_ref[which(imgt_ref$Allele == V.allele.name), ] %>% dplyr::distinct(seq.nt)
  imgt_ref[which(imgt_ref$Allele == J.allele.name), ]

  # imgt ref seqs
  V.allele.seq <- imgt_ref[which(imgt_ref$Allele == V.allele.name),] %>% dplyr::distinct(seq.nt, meta, Allele)
  if (nrow(V.allele.seq) > 1) {
    if (sum(grepl("Mus musculus", V.allele.seq$meta)) == 1) {
      V.allele.seq <- V.allele.seq[which(grepl("Mus musculus", V.allele.seq$meta)),]
    } else {
      V.allele.seq <- V.allele.seq[sample(1:nrow(V.allele.seq), 1),]
    }
  }
  V.allele.seq <- V.allele.seq[,"seq.nt",drop=T]
  names(V.allele.seq) <- V.allele.name

  J.allele.seq <- imgt_ref[which(imgt_ref$Allele == J.allele.name),] %>% dplyr::distinct(seq.nt, meta, Allele)
  if (nrow(J.allele.seq) > 1) {
    if (sum(grepl("Mus musculus", J.allele.seq$meta)) == 1) {
      J.allele.seq <- J.allele.seq[which(grepl("Mus musculus", J.allele.seq$meta)),]
    } else {
      J.allele.seq <- J.allele.seq[sample(1:nrow(J.allele.seq), 1),]
    }
  }
  J.allele.seq <- J.allele.seq[,"seq.nt",drop=T]
  names(J.allele.seq) <- J.allele.name

  if (missing(C_allele)) {
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
  p2 <- lapply(seq_along(V.allele.seq), function (x) {
    MultiplePairwiseAlignmentsToOneSubject(subject = Biostrings::DNAStringSet(cs[x]), patterns = Biostrings::DNAStringSet(c(V.allele.seq[x], J.allele.seq[x], C.allele.seq)), type = "local", pattern.positions.size = 4, print.subject.min.max = T)
  })
  names(p2) <- paste0(V.allele.name,"_",J.allele.name)

  p3 <- lapply(seq_along(V.allele.seq), function (x) {
    MultiplePairwiseAlignmentsToOneSubject(subject = Biostrings::DNAStringSet(cs[x]), patterns = Biostrings::DNAStringSet(c(V.allele.seq[x], J.allele.seq[x])), type = "local", pattern.positions.size = 4, print.subject.min.max = T)
  })
  names(p3) <- paste0(V.allele.name,"_",J.allele.name)

  TCR.seq <- unlist(lapply(seq_along(V.allele.seq), function (x) {
    stringr::str_sub(cs[x], start = p3[[x]][["min.max.subject.position"]][1], end = p3[[x]][["min.max.subject.position"]][2])
  }))
  names(TCR.seq) <- paste0(TCR, "_", chain, "(", V.allele.name, "_", J.allele.name, ")")

  return(list("consensus.alingments" = p1, "consensus.seqs.uncut" = cs, "VJC.alignments" = p2, "VJ.alignments" = p3, "V.to.J.TCR.seqs" = TCR.seq))
}
