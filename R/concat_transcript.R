#' Pull and concatenate exon ranges into one sequence
#'
#' Based on gtf_df which should only contain lines of exon definitions for
#' one gene, the respective sequences are pulled from refseq and concatenated
#' into the coding sequence (CDS).
#'
#' https://www.biostars.org/p/3423/
#'
#' @param gtf_df subsetted gtf data frame
#' @param refseq reference sequence to pull from
#'
#' @return
#' @export
#'
#' @examples
concat_transcript <- function(gtf_df,
                         refseq,
                         refseq_strand = c("+", "-")) {

  ## argument
  refseq_strand <- match.arg(refseq_strand, c("+", "-")) # test "for" and "rev"

  # check for needed columns start and end
  if (any(!c("exon_number", "seqname", "start", "end", "transcript_id", "feature", "strand") %in% names(gtf_df))) {
    stop("gtf_df has to have the columns at least: seqname, exon_number, start, end, transcript_id, feature, strand.")
  }
  # check for exon lines only
  # if (any(gtf_df[["feature"]] != "exon")) {
  #   stop("gtf_df should only have lines with feature == 'exon'.")
  # }
  # check for uniqueness of transcript_id
  if (length(unique(gtf_df[["transcript_id"]])) > 1) {
    stop("More than one transcript_id found in gtf_df.")
  }
  # check equality of seqname in gtf_df and name of refseq
  if (!is.null(names(refseq))) {
    if (names(refseq) != unique(gtf_df[["seqname"]])) {
      message("name of refseq and seqname in gtf_df are unequal: ", names(refseq), " vs. ", unique(gtf_df[["seqname"]]), ".")
    }
  }
  # print message if strand is minus; then alignment against sequences from ncbi requires igsc:::revcompDNA
  # if ("strand" %in% names(gtf_df) && unique(gtf_df[["strand"]]) == "-") {
  #   message("gene is on minus strand. for alignment against sequences from e.g. NCBI, the reverse complement is require, e.g. with igsc:::revcompDNA.")
  # }

  # check diffs on - and + strand, name utrs by 3' and 5', label seqlist or attribute with 3' and 5' end
  # what are antisense genes

  # for minus strand: order everything in reverse (for intron range inference below!)
  # for plus strand: check
  strand <- unique(gtf_df$strand)
  gtf_df <- gtf_df[order(gtf_df$start, decreasing = all(gtf_df$start > gtf_df$end)),]
  seqpos_intron <- stats::setNames(igsc:::seq2(gtf_df_exon$end[-length(gtf_df_exon$end)] + 1 , gtf_df_exon$start[-1] - 1),
                                   paste0("intron", sprintf("%02d", as.numeric(gtf_df_exon$exon_number[-1]))))
  intron_ranges <- lapply(seqpos_intron, range)
  if (strand == "-") {
    start_fun <- match.fun("max")
    end_fun <- match.fun("min")
  } else if (strand == "+") {
    start_fun <- match.fun("min")
    end_fun <- match.fun("max")
  }
  gtf_df_intron <- data.frame(feature = "intron",
                              feature2 = names(intron_ranges),
                              start = unname(sapply(intron_ranges, start_fun)),
                              end = unname(sapply(intron_ranges, end_fun)))
  gtf_df <- dplyr::bind_rows(gtf_df, gtf_df_intron)

  seqnames <- paste0(gtf_df$feature, gsub("NA", "", sprintf("%02d", as.numeric(gtf_df$exon_number))))
  seqnames[which(grepl("codon", seqnames))] <- stringr::str_sub(seqnames[which(grepl("codon", seqnames))],1,-3)
  # change utr naming later
  # add introns to gtf_df
  seqlist <- mapply(substr, start = stats::setNames(gtf_df$start, seqnames), stop = gtf_df$end, x = refseq)
  if (refseq_strand != unique(gtf_df$strand)) {
    seqlist <- lapply(seqlist, igsc:::revcompDNA)
  }
  seqpos <- stats::setNames(igsc:::seq2(gtf_df$start, gtf_df$end), seqnames)
  gtf_df$feature2 <- seqnames
  # introns
  gtf_df_exon <- gtf_df[which(gtf_df$feature == "exon"),]
  gtf_df_exon$exon_number[-length(gtf_df_exon$exon_number)]
  seqpos <- c(seqpos, seqpos_intron)

  ## add seq for plus and minus? position for plus and minus?
  subject_df <- data.frame(seq = strsplit(seqlist[["transcript"]], "")[[1]],
                          position_genome = seqpos[["transcript"]],
                          position = seqpos[["transcript"]] - min(seqpos[["transcript"]]) + 1)
  intron_exon_df <-
    do.call(rbind, seqpos_stack[which(grepl("exon|intron", names(seqpos_stack)))]) %>%
    dplyr::rename("intron_exon" = pattern)
  intron_exon_df$intron_exon <- stringr::str_sub(intron_exon_df$intron_exon, 1, -3)
  subject_df <- dplyr::left_join(subject_df, intron_exon_df, by = "position_genome")

  seqpos_stack <- list()
  j <- 0
  for (i in which(names(seqpos) != "transcript")) {
    j <- j + 1
    seqpos_stack[[j]] <- utils::stack(seqpos[i])
    names(seqpos_stack[[j]]) <- c("position_genome", "pattern")
    seqpos_stack[[j]]$pattern <- as.character(seqpos_stack[[j]]$pattern)
  }
  names(seqpos_stack) <- names(seqpos)[which(names(seqpos) != "transcript")]
  seqpos_df <-
    dplyr::bind_rows(seqpos_stack) %>%
    tidyr::nest(.key = "pattern", .by = position_genome)
  subject_df <- dplyr::left_join(subject_df, seqpos_df, by = "position_genome")
  subject_df_unnest <-
    subject_df %>%
    tidyr::unnest(cols = pattern)



  ### transcript
  seq_transcript <- seqlist[["transcript"]]



  #seqpos_stack2 <- dplyr::bind_rows(seqpos_stack[which(!grepl("exon|intron", names(seqpos_stack)))])
  #names(seqpos_stack2)[2] <- "pattern2"
  #algnmt_df <- dplyr::left_join(algnmt_df, seqpos_stack2, by = "position_genome")

  #sort(table(seqpos_stack2$position_genome), decreasing = T)
  #mpa <- MultiplePairwiseAlignmentsToOneSubject(subject = seqlist["transcript"], patterns = seqlist[which(names(seqlist) != "transcript")])






  if ("CDS" %in% gtf_df$feature) {
    seqlist_CDS <- seqlist[which(grepl("CDS", names(seqlist)))]
    seq_CDS <- paste(seqlist_CDS, collapse = "")
    nt_cum_CDS <- cumsum(nchar(seqlist_CDS))
    starts_CDS <- c(0, nt_cum_CDS) + 1
    ranges_CDS <- paste0(starts_CDS[-length(starts_CDS)], "..", nt_cum_CDS)

    names(ranges_CDS) <- paste0("CDS", sprintf("%02d", seq_along(ranges_CDS)))
    nt_CDS <- nchar(seqlist_CDS)
    codon_cum_CDS <- cumsum(nt_CDS/3)
    intron_phase <- ifelse(dplyr::near(codon_cum_CDS[-length(codon_cum_CDS)] %% 1, 0), "0",
                           ifelse(dplyr::near(codon_cum_CDS[-length(codon_cum_CDS)] %% 1, 1/3), "1", "2"))
    names(intron_phase) <- gsub("CDS", "intron", names(intron_phase))

    gtf_df_CDS <- gtf_df[which(gtf_df$feature == "CDS"),]
    nt_intron <- lengths(igsc:::seq2(gtf_df_CDS$start[-length(gtf_df_CDS$start)], gtf_df_CDS$end[-1])) - 2
    nt_sum_intron <- sum(nt_intron)

    #seqlist[which(grepl("codon", names(seqlist)))]
    if (any(grepl("codon", gtf_df$feature))) {
      gtf_codons <- gtf_df[which(grepl("codon", gtf_df$feature)),]
      min_CDS_pos <- min(unlist(gtf_df[which(grepl("CDS", gtf_df$feature)),c("start", "end")]))
      min_CDS_pos <- min_CDS_pos - 3 # by definition here, start codon is not part of CDS

      ## think about this. +/- stand?? max_CDS_pos
      gtf_codons[which(gtf_codons$feature == "start_codon"),"start"] <- gtf_codons[which(gtf_codons$feature == "start_codon"),"start"] - min_CDS_pos + 1 - nt_sum_intron
      gtf_codons[which(gtf_codons$feature == "start_codon"),"end"] <- gtf_codons[which(gtf_codons$feature == "start_codon"),"end"] - min_CDS_pos + 1 - nt_sum_intron
      gtf_codons[which(gtf_codons$feature == "end_codon"),"start"] <- gtf_codons[which(gtf_codons$feature == "end_codon"),"start"] - min_CDS_pos + 1
      gtf_codons[which(gtf_codons$feature == "end_codon"),"end"] <- gtf_codons[which(gtf_codons$feature == "end_codon"),"end"] - min_CDS_pos + 1
    }
    attributes(seq_CDS) <- list(range = ranges_CDS, nt = nt_CDS, nt_cumsum = nt_cum_CDS, codon_cumsum = codon_cum_CDS, intron_phase = intron_phase)
  }

  if ("exon" %in% gtf_df$feature) {
    seqlist_exon <- seqlist[which(grepl("exon", names(seqlist)))]
    seq_exon <- paste(seqlist_exon, collapse = "")
    nt_cum_exon <- cumsum(nchar(seqlist_exon))
    starts_exon <- c(0, nt_cum_exon) + 1
    ranges_exon <- paste0(starts_exon[-length(starts_exon)], "..", nt_cum_exon)

    names(ranges_exon) <- paste0("exon", sprintf("%02d", seq_along(ranges_exon)))
    nt_exon <- nchar(seqlist_exon)
    attributes(seq_exon) <- list(range = ranges_exon, nt = nt_exon, nt_cum = nt_cum_exon)
  }


  if ("transcript" %in% names(seqlist)) {
    seq <- seqlist[["transcript"]]
  }

  igsc:::revcompDNA(seqlist[["CDS01"]])

  seq <- paste(seqlist, collapse = "")
  exon_cum <- cumsum(nchar(seqlist))
  starts <- c(0, exon_cum) + 1
  exon_ranges <- paste0(starts[-length(starts)], "..", exon_cum)
  exon_names <- paste0("exon", sprintf("%02d", as.numeric(gtf_df$exon_number)))
  names(exon_ranges) <- exon_names
  exon_lengths <- nchar(seqlist)
  names(exon_lengths) <- exon_names
  n_codons <- cumsum(exon_lengths/3)
  names(n_codons) <- exon_names
  attributes(seq) <- list(ranges = exon_ranges, lengths = exon_lengths, codons = n_codons)
  names(seqlist) <- exon_names

  return(list(seq = seq, exons = seqlist))
}


