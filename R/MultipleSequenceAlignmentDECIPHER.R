#' Have a Multiple Sequence Alignment calculated with DECIPHER, plotted as ggplot2 object
#'
#' @param input.set Vector of DNA sequences as character or DNAStringSet()
#' @param add.consensus.seq logicle of a consensus sequence should be added to the alignment
#' @param print.disamb.cons.seq.border logicle if to print the borders of longest stretch of disambiguate consensus sequence
#' @param collapse.duplicate.seqs collapse identical sequences in the input.set to one sequence before alignemnt
#' @param open.browser open the browser via DECIPHERs command DECIPHER::BrowseSeqs in addition to having the plot returned
#'
#' @return a ggplot object of the multiple sequence alignment
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#'
#' }
MultipleSequenceAlignmentDECIPHER <- function(input.set,
                                              add.consensus.seq = T,
                                              print.disamb.cons.seq.border = T,
                                              collapse.duplicate.seqs = T,
                                              open.browser = F) {
  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"Biostrings" %in% rownames(utils::installed.packages())) {BiocManager::install("Biostrings")}
  if (!"DECIPHER" %in% rownames(utils::installed.packages())) {BiocManager::install("DECIPHER")}


  alignment.colour.palette <- c("A" = "#ffafaf",
                                "T" = "#fed7af",
                                "C" = "#afffaf",
                                "G" = "#afffff",
                                "-" = "#ffffff",
                                "match" = "#999999",
                                "mismatch" = "#a51515",
                                "gap" = "#9248d4",
                                "insertion" = "#000000",
                                "ambiguous" = "#E69F00")

  if (add.consensus.seq) {
    if ("consensus" %in% names(input.set)) {
      stop("consensus cannot be the name of either of the input sequences.")
    }
    if (is.character(input.set)) {
      cons.seq <- DECIPHER::ConsensusSequence(DECIPHER::AlignSeqs(Biostrings::DNAStringSet(input.set), verbose = F))
    } else {
      cons.seq <- DECIPHER::ConsensusSequence(DECIPHER::AlignSeqs(input.set, verbose = F))
    }
  }

  if (is.character(input.set)) {
    if (is.null(names(input.set))) {
      names(input.set) <- seq(1,length(input.set),1)
    }
    if (collapse.duplicate.seqs) {
      if (length(unique((input.set))) != length(input.set))  {
        unique.input.set <- unique(input.set)
        n.dup <- sapply(unique.input.set, function(x) {length(which(x == input.set))})
        names(unique.input.set) <- paste("unique seq", seq(1,length(unique.input.set)), " (", n.dup, ") ", sep = "")
        input.set <- unique.input.set
      }
    }
    if (length(input.set) > 1) {
      input.set <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(input.set), verbose = F)
    } else {
      input.set <- Biostrings::DNAStringSet(input.set)
    }
  }

  if (add.consensus.seq) {
    input.set[length(input.set)+1] <- cons.seq
    names(input.set)[length(input.set)] <- "consensus"
  }

  if (open.browser) {
    DECIPHER::BrowseSeqs(input.set)
  }

  out <- utils::stack(strsplit(as.character(input.set), ""))
  names(out) <- c("seq", "seq.name")
  out$position <- rep(1:nchar(as.character(input.set[1])), length(input.set))
  ## best way here is to use group_by from dplyr, even though a loop would also work
  out <-
    out %>%
    dplyr::group_by(position) %>%
    dplyr::mutate(case = dplyr::case_when(nlevels(as.factor(seq)) == 1 ~ "match",
                                          (nlevels(as.factor(seq)) > 1 & !"-" %in% seq) ~ "mismatch",
                                          (nlevels(as.factor(seq)) > 1 & "-" %in% seq) ~ "gap")) %>%
    dplyr::ungroup()

  out$case <- ifelse(out$seq.name == "consensus", out$seq, out$case)
  out$case <- ifelse(out$case == "gap" & out$seq.name != "consensus", out$seq, out$case)
  out$case <- ifelse(out$case %in% c(names(Biostrings::IUPAC_CODE_MAP[-c(1:4)]), "+"), "ambiguous", out$case)
  out$case <- factor(out$case, levels = names(alignment.colour.palette))
  out$seq.name <- factor(out$seq.name, levels = c("consensus", names(input.set)[-length(input.set)]))

  g2 <- ggplot2::ggplot(out, ggplot2::aes(x = position, y = seq.name, fill = case)) +
    ggplot2::geom_raster() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = alignment.colour.palette)

  if (print.disamb.cons.seq.border & add.consensus.seq) {
    cons.seq.disam <- consecutive.disambiguate.consensus.seq.from.decipher.consensus(cons.seq)
    cons.seq.disam.pos <- Biostrings::pairwiseAlignment(subject = cons.seq, pattern = cons.seq.disam, type = "local")
    cons.seq.disam.pos <- c(cons.seq.disam.pos@subject@range@start, cons.seq.disam.pos@subject@range@start + cons.seq.disam.pos@subject@range@width - 1)
    cons.seq.disam.pos <- data.frame(pos = cons.seq.disam.pos, name = c("dismab.cons.seq.low.lim", "dismab.cons.seq.up.lim"))
    g2 <- g2 + ggplot2::geom_vline(data = cons.seq.disam.pos,  ggplot2::aes(xintercept = pos), linetype = "dashed") # + scale_color_manual(values = custom.scale[c(5,8)]
  }

  return(g2)
}


consecutive.disambiguate.consensus.seq.from.decipher.consensus <- function (consensus) {
  # returns the longest stretch of disambiguate nucleotides from a multiple sequence alignment made with DECIPHER::AlignSeqs
  if (class(consensus) != "DNAStringSet") {
    stop("Please provide a consesus sequence as DNAStringSet, generated with DECIPHER::ConsensusSequence")
  }
  return(names(which.max(sapply(strsplit(as.character(consensus), paste(c(names(Biostrings::IUPAC_CODE_MAP[-c(1:4)])), collapse = "|"))[[1]], nchar))))
}
