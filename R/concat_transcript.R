#' Pull and concatenate exon ranges into one sequence
#'
#' Based on gtf_df which should only contain lines of exon definitions for
#' one gene, the respective sequences are pulled from refseq and concatenated
#' into the coding sequence (CDS).
#' Stop codon added to CDS - add option
#'
#' https://www.biostars.org/p/3423/, https://pmc.ncbi.nlm.nih.gov/articles/PMC6099125/, https://en.wikipedia.org/wiki/Sense_(molecular_biology)
#'
#' @param gtf_df subsetted gtf data frame
#' @param refseq reference sequence in 5' to 3' direction, genomic sequence: plus or minus strand
#' @param refseq_strand which strand does refseq represent?
#' @param stop_codon_to_CDS include stop codon in CDS?
#' @param run_test compare returned sequences to references?
#' @param ... arguments to get_genome_seq
#'
#' @return
#' @export
#'
#' @examples
concat_transcript <- function(gtf_df,
                              refseq,
                              refseq_strand = c("+", "-"),
                              stop_codon_to_CDS = F,
                              run_test = F,
                              ...) {


  # # CD8A - strand
  # gtf_df <- readRDS("/Users/vonskopnik/Documents/2024_igsc_testing/cd8a_gtf.rds")
  # refseq <- readRDS("/Users/vonskopnik/Documents/2024_igsc_testing/chr2_seq.rds")
  # refseq_strand <- "+"
  # stop_codon_to_CDS <- F
  #
  # # CD4 + strand
  # gtf_df <- readRDS("/Users/vonskopnik/Documents/2024_igsc_testing/cd4_gtf.rds")
  # refseq <- readRDS("/Users/vonskopnik/Documents/2024_igsc_testing/chr12_seq.rds")
  # refseq_strand <- "+"
  # stop_codon_to_CDS <- F

  # check for needed columns start and end
  if (any(!c("exon_number", "seqname", "start", "end", "transcript_id", "feature", "strand") %in% names(gtf_df))) {
    stop("gtf_df has to have the columns at least: seqname, exon_number, start, end, transcript_id, feature, strand.")
  }
  # check for uniqueness of transcript_id
  if (length(unique(gtf_df[["transcript_id"]])) > 1) {
    stop("More than one transcript_id found in gtf_df.")
  }
  # check for uniqueness of transcript_id
  if (length(unique(gtf_df[["strand"]])) > 1) {
    stop("Strand column has more than one level in gtf_df.")
  }

  feat1 <- sapply(c("transcript", "exon", "CDS", "start_codon", "stop_codon", "UTR"), grepl, x = unique(gtf_df[["feature"]]), simplify = F)
  if (any(!sapply(feat1, any))) {
    stop("feature columns needs to have at least one of each of these: transcript exon CDS start_codon stop_codon UTR")
  }

  # check equality of seqname in gtf_df and name of refseq (from fasta file)
  if (!is.null(names(refseq))) {
    if (names(refseq) != unique(gtf_df[["seqname"]])) {
      message("name of refseq and seqname in gtf_df are unequal: ", names(refseq), " vs. ", unique(gtf_df[["seqname"]]), ".")
    }
  }

  # Start/end: In both cases, gene on plus or minus strand, end is larger
  # then start. This makes immediate sense for the plus strand but for the minus
  # strand, genes actually go from a larger to a smaller position (5'->3').

  refseq_strand <- match.arg(refseq_strand, c("+", "-"))
  # always pull sequences from the plus strand, as start and end refer to that strand
  # if one applies revcomp to the plus strand to have the minus strand, positions do no longer match
  if (refseq_strand == "-") {
    refseq <- revcompDNA(refseq, fun = "rcpp")
  }

  gtf_df$exon_number <- as.numeric(gtf_df$exon_number)
  strand <- unique(gtf_df$strand)
  gtf_df <- add_introns(gtf_df = gtf_df, strand = strand, n_exon = max(gtf_df$exon_number, na.rm = T))
  # add alternative feature names
  minpos <- min(gtf_df$start, gtf_df$end)
  gtf_df <-
    gtf_df %>%
    dplyr::mutate(feature2 = paste0(feature, gsub("NA", "", sprintf("%02d", exon_number))), .after = feature) %>%
    dplyr::mutate(feature2 = ifelse(grepl("codon", feature2), stringr::str_sub(feature2, 1, -3), feature2)) %>%
    dplyr::mutate(start_transcript = start - minpos + 1, end_transcript = end - minpos + 1, .after = end)
  if (strand == "+") {
    gtf_df <- dplyr::arrange(gtf_df, start)
  } else if (strand == "-") {
    gtf_df <- dplyr::bind_rows(gtf_df[which(gtf_df$feature == "transcript"),],
                               dplyr::arrange(gtf_df[-which(gtf_df$feature == "transcript"),], -start))
  }

  # change utr naming
  start_codon_start <- gtf_df[which(gtf_df$feature == "start_codon"), "start",drop=T]
  if (strand == "+") {
    UTR_order <- c("5", "3")
  } else if (strand == "-") {
    UTR_order <- c("3", "5")
  }
  gtf_df[intersect(which(gtf_df$start < start_codon_start), which(gtf_df$feature == "UTR")),"feature2"] <- paste0(UTR_order[1], gtf_df[intersect(which(gtf_df$start < start_codon_start), which(gtf_df$feature == "UTR")),"feature2",drop=T])
  gtf_df[intersect(which(gtf_df$start > start_codon_start), which(gtf_df$feature == "UTR")),"feature2"] <- paste0(UTR_order[2], gtf_df[intersect(which(gtf_df$start > start_codon_start), which(gtf_df$feature == "UTR")),"feature2",drop=T])


  # get sequences from refseq
  #seqlist <- mapply(substr, start = stats::setNames(gtf_df$start, gtf_df$feature2), stop = gtf_df$end, x = refseq)
  seqlist <- purrr::map2_chr(.x = stats::setNames(gtf_df$start, gtf_df$feature2), .y = gtf_df$end, substr, x = refseq)
  pos_gen <- stats::setNames(igsc:::seq2(gtf_df$start, gtf_df$end), gtf_df$feature2)
  #pos_gen_rel <- stats::setNames(igsc:::seq2(gtf_df$start_transcript, gtf_df$end_transcript), gtf_df$feature2)

  # sequences are always derived from plus strand (see above)
  # so if the gene is on minus strand, revcomp or rev is required for sequences and positions, respectively
  if (strand == "-") {
    seqlist <- revcompDNA(seqlist, fun = "rcpp")
    pos_gen <- lapply(pos_gen, rev)
    #pos_gen_rel <- lapply(pos_gen_rel, rev)
  }

  direction_5to3 <- "5' --> 3'"
  attributes(seqlist) <- list(names = names(seqlist), strand = strand, direction = direction_5to3)

  ## add seq for plus and minus strand
  transcript_df <- data.frame(seq1 = strsplit(seqlist[["transcript"]], "")[[1]],
                              #seq2 = strsplit(revcompDNA(seqlist[["transcript"]], rev = F), "")[[1]],
                              position_genome = pos_gen[["transcript"]],
                              position = sort(pos_gen[["transcript"]] - min(pos_gen[["transcript"]]) + 1)) # , decreasing = decr_pos_gen
  names(transcript_df)[1] <-
    if (strand == "+") {
      c("seq_plus_5to3")
      #c("seq_plus_5to3", "seq_minus_3to5")
    } else {
      c("seq_minus_5to3")
      #c("seq_minus_5to3", "seq_plus_3to5")
    }
  strand_coding <- names(transcript_df)[1]
  gtf_df[[names(transcript_df)[1]]] <- seqlist
  #gtf_df[[names(transcript_df)[2]]] <- revcompDNA(seqlist)

  pos_gen_stack <- make_pos_gen_stack(pos_gen)

  # gene on plus strand: position and position_genome both increase
  # gene on misnus strand: to have left-to-right reading, position increase while position_genome decreases
  transcript_df <- dplyr::left_join(transcript_df,
                                    dplyr::bind_rows(pos_gen_stack[which(grepl("exon|intron", names(pos_gen_stack)))]) %>%
                                      dplyr::rename("intron_exon" = pattern) %>%
                                      dplyr::mutate(intron_exon = stringr::str_sub(intron_exon, 1, -3)),
                                    by = "position_genome")
  transcript_df_unnest <-
    dplyr::left_join(transcript_df,
                     dplyr::bind_rows(pos_gen_stack[which(names(pos_gen_stack) != "transcript")]) %>%
                       tidyr::nest(.key = "pattern", .by = position_genome),
                     by = "position_genome") %>%
    tidyr::unnest(cols = pattern)
  #transcript_df <- dplyr::distinct(transcript_df_unnest, seq_plus, seq_minus, position_genome, position, intron_exon)
  algnmt_df <-
    dplyr::left_join(transcript_df,
                     dplyr::bind_rows(pos_gen_stack) %>%
                       tidyr::nest(.key = "pattern", .by = position_genome), by = "position_genome") %>%
    tidyr::unnest(cols = pattern)

  ### transcript
  groups <- c("UTR", "codon", "exon", "CDS", "intron")
  groups2 <- c("exon", "CDS")
  ranges_transcript <- lapply(stats::setNames(groups, paste0("range_", groups)), function(x) {
    df <- transcript_df_unnest[which(grepl(x, transcript_df_unnest$pattern)),]
    df <- df %>% dplyr::group_by(pattern) %>% dplyr::summarise(range = paste0(min(position), "..", max(position)))
    return(stats::setNames(df$range, df$pattern))
  })
  nts <- lapply(stats::setNames(groups, paste0("nt_", groups)), function(x) {
    nchar(seqlist[which(grepl(x, names(seqlist)))])
  })
  nts_cumsum <- lapply(stats::setNames(groups2, paste0("nt_cumsum_", groups2)), function(x) {
    cumsum(nchar(seqlist[which(grepl(x, names(seqlist)))]))
  })
  codon_cumsum <- nts_cumsum$nt_cumsum_CDS/3
  intron_phase <- ifelse(dplyr::near(codon_cumsum[-length(codon_cumsum)] %% 1, 0), "0",
                         ifelse(dplyr::near(codon_cumsum[-length(codon_cumsum)] %% 1, 1/3), "1", "2"))
  names(intron_phase) <- gsub("CDS", "intron", names(intron_phase))
  seq_transcript <- seqlist[["transcript"]]
  attributes(seq_transcript) <- c(list(strand = strand, direction = direction_5to3),
                                  ranges_transcript, nts, nts_cumsum,
                                  list(codon_cumsum = codon_cumsum, intron_phase = intron_phase))

  df0_transcript <-
    purrr::map_dfr(ranges_transcript, utils::stack) %>%
    dplyr::relocate(ind, values) %>%
    dplyr::mutate(ind = as.character(ind)) %>%
    dplyr::rename(value = values, feature = ind) %>%
    tidyr::separate(value, into = c("start", "end"), sep = "\\.\\.") %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))
  df1 <-
    purrr::map_dfr(nts, utils::stack) %>%
    dplyr::mutate(ind = as.character(ind)) %>%
    dplyr::rename(nt = values, feature = ind)
  df2 <-
    purrr::map_dfr(nts_cumsum, utils::stack) %>%
    dplyr::mutate(ind = as.character(ind)) %>%
    dplyr::rename(nt_cumsum = values, feature = ind) %>%
    dplyr::mutate(temp = stringr::str_sub(feature, 1, -3))
  df2 <- split(df2, df2$temp)
  df2 <- purrr::map(df2, function(x) {
    names(x)[1] <- paste0(names(x)[1], "_", x[["temp"]][1])
    x <- x[,-3]
    return(x)
  })
  df3 <-
    utils::stack(codon_cumsum) %>%
    dplyr::mutate(ind = as.character(ind), values = as.numeric(values)) %>%
    dplyr::rename(codon_cumsum = values, feature = ind)
  df4 <-
    utils::stack(intron_phase) %>%
    dplyr::mutate(ind = as.character(ind), values = as.integer(values)) %>%
    dplyr::rename(intron_phase = values, feature = ind)
  df_transcript_meta <-
    df0_transcript %>%
    dplyr::left_join(df1, by = "feature") %>%
    dplyr::left_join(df2[[1]], by = "feature") %>%
    dplyr::left_join(df2[[2]], by = "feature") %>%
    dplyr::left_join(df3, by = "feature") %>%
    dplyr::left_join(df4, by = "feature")
  # df_transcript_meta_long <-
  #   tidyr::pivot_longer(df_transcript_meta, cols = -feature) %>%
  #   tidyr::drop_na()
  # tidyr::pivot_wider(df_transcript_meta_long) # reverse

  # in df_transcript_meta: leave UTRs with trailing exon number as UTR may spread over more than one exon
  # in df_exon_meta the trailing number can be removed, as without introns the UTR from multiple exons becomes adjacent
  # in CDS, UTR are irrelevant

  ### exon
  exon_df_nest <-
    transcript_df_unnest %>%
    dplyr::filter(intron_exon == "exon") %>%
    # remove trailing UTR number here to collapse UTR across multiple exon (e.g. in CD4 gene)
    dplyr::mutate(pattern = ifelse(grepl("UTR", pattern), stringr::str_sub(pattern, 1, -3), pattern)) %>%
    tidyr::nest(.key = "pattern", .by = names(transcript_df)) %>%
    dplyr::mutate(position = dplyr::row_number())
  seq_exon <- paste(exon_df_nest[[strand_coding]], collapse = "")
  exon_df_nest <- tidyr::unnest(exon_df_nest, cols = pattern)
  #max(diff(exon_df_nest$position))
  groups <- c("UTR", "codon", "exon", "CDS")
  ranges_exon <- lapply(stats::setNames(groups, paste0("range_", groups)), function(x) {
    df <- exon_df_nest[which(grepl(x, exon_df_nest$pattern)),]
    df <- df %>% dplyr::group_by(pattern) %>% dplyr::summarise(range = paste0(min(position), "..", max(position)))
    return(stats::setNames(df$range, df$pattern))
  })
  attributes(seq_exon) <- c(list(strand = strand, direction = direction_5to3),
                            ranges_exon, nts[which(names(nts) != "nt_intron")], nts_cumsum,
                            list(codon_cumsum = codon_cumsum, intron_phase = intron_phase))

  df0_exon <-
    purrr::map_dfr(ranges_exon, utils::stack) %>%
    dplyr::relocate(ind, values) %>%
    dplyr::mutate(ind = as.character(ind)) %>%
    dplyr::rename(value = values, feature = ind) %>%
    tidyr::separate(value, into = c("start", "end"), sep = "\\.\\.") %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))
  df_exon_meta <-
    df0_exon %>%
    dplyr::left_join(df1, by = "feature") %>%
    dplyr::left_join(df2[[1]], by = "feature") %>%
    dplyr::left_join(df2[[2]], by = "feature") %>%
    dplyr::left_join(df3, by = "feature") %>%
    dplyr::mutate(nt = ifelse(is.na(nt), end - start + 1, nt))
  # df_exon_meta_long <-
  #   tidyr::pivot_longer(df_exon_meta, cols = -feature) %>%
  #   tidyr::drop_na()

  ### CDS
  if (stop_codon_to_CDS) {
    grep_str <- "CDS|stop_codon"
  } else {
    grep_str <- "CDS"
  }
  CDS_pos <- transcript_df_unnest[which(grepl(grep_str, transcript_df_unnest$pattern)),"position",drop=T]
  CDS_df_nest <-
    transcript_df_unnest %>%
    dplyr::filter(position %in% CDS_pos) %>%
    dplyr::filter(!grepl("exon", pattern)) %>%
    tidyr::nest(.key = "pattern", .by = names(transcript_df)) %>%
    dplyr::mutate(position = dplyr::row_number())
  seq_CDS <- paste(CDS_df_nest[[strand_coding]], collapse = "")
  CDS_df_nest <- tidyr::unnest(CDS_df_nest, cols = pattern)

  groups <- c("codon", "CDS")
  ranges_CDS <- lapply(stats::setNames(groups, paste0("range_", groups)), function(x) {
    df <- CDS_df_nest[which(grepl(x, CDS_df_nest$pattern)),]
    df <- df %>% dplyr::group_by(pattern) %>% dplyr::summarise(range = paste0(min(position), "..", max(position)))
    return(stats::setNames(df$range, df$pattern))
  })

  attributes(seq_CDS) <- c(list(strand = strand, direction = direction_5to3),
                           ranges_CDS, nts[which(names(nts) %in% c("nt_CDS", "nt_codon"))],
                           nts_cumsum[which(names(nts_cumsum) %in% c("nt_cumsum_CDS"))],
                           list(codon_cumsum = codon_cumsum, intron_phase = intron_phase))

  df0_CDS <-
    purrr::map_dfr(ranges_CDS, utils::stack) %>%
    dplyr::relocate(ind, values) %>%
    dplyr::mutate(ind = as.character(ind)) %>%
    dplyr::rename(value = values, feature = ind) %>%
    tidyr::separate(value, into = c("start", "end"), sep = "\\.\\.") %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))
  df_CDS_meta <-
    df0_CDS %>%
    dplyr::left_join(df1, by = "feature") %>%
    dplyr::left_join(df2[[1]], by = "feature") %>%
    dplyr::left_join(df3, by = "feature")
  # df_CDS_meta_long <-
  #   tidyr::pivot_longer(df_CDS_meta, cols = -feature) %>%
  #   tidyr::drop_na()

  #mpa <- MultiplePairwiseAlignmentsToOneSubject(subject = seqlist["transcript"], patterns = seqlist[which(names(seqlist) != "transcript")])

  # reduce redundancy in gtf_df before return
  cols_reduce <- c("seqname", "source", "gene_id", "gene_name", "transcript_id", "transcript_name")
  redundant_info <- unlist(gtf_df[1,cols_reduce])
  names(redundant_info)[1] <- "seqname_gtf"
  redundant_info <- c(redundant_info, seqname_refseq = names(refseq))

  # how to elegantly add seq_CDS to CDS_df_nest as putative ref; removal is easy outside the function
  CDS_df_nest <-
    dplyr::bind_rows(CDS_df_nest,
                     dplyr::distinct(CDS_df_nest, dplyr::across(-pattern)) %>%
                       dplyr::mutate(pattern = "CDS")) %>%
    dplyr::arrange(position, pattern)
  exon_df_nest <-
    dplyr::bind_rows(exon_df_nest,
                     dplyr::distinct(exon_df_nest, dplyr::across(-pattern)) %>%
                       dplyr::mutate(pattern = "exon")) %>%
    dplyr::arrange(position, pattern)

  return_list <- c(
    list(pre_mRNA = list(seq = seq_transcript, features = df_transcript_meta, seq_df = algnmt_df),
         mRNA = list(seq = seq_exon, features = df_exon_meta, seq_df = exon_df_nest),
         CDS = list(seq = seq_CDS, features = df_CDS_meta, seq_df = CDS_df_nest),
         seqlist = seqlist,
         gtf2 = gtf_df[,-which(names(gtf_df) %in% c(cols_reduce, "score"))]),
    as.list(redundant_info))

  # check transcript seq with get_genome_seq

  if (run_test) {
    trans_range <- range(return_list[["pre_mRNA"]][["seq_df"]] %>%
                           dplyr::filter(pattern == "transcript") %>%
                           #dplyr::arrange(position_genome) %>%
                           dplyr::pull(position_genome))
    df_seq <- paste(return_list[["pre_mRNA"]][["seq_df"]] %>%
                      dplyr::filter(pattern == "transcript") %>%
                      #dplyr::arrange(position_genome) %>%
                      dplyr::pull(!!rlang::sym(strand_coding)), collapse = "")
    test_seq <- get_genome_seq(chromosome = gtf_df$seqname[1],
                               start = trans_range[1],
                               end = trans_range[2],
                               ...)
    if (strand == "-") {
      test_seq <- revcompDNA(test_seq)
    }

    out <-
      return_list[["pre_mRNA"]][["seq_df"]] %>%
      dplyr::filter(grepl("codon", pattern))

    if (paste(out[which(out$pattern == "start_codon"), 1, drop = T], collapse = "") != "ATG") {
      message("start codon on coding strand is not ATG.")
    }
    if (!paste(out[which(out$pattern == "stop_codon"), 1, drop = T], collapse = "") %in% c("TGA", "TAG", "TAA")) {
      message("stop codon on coding strand is neither of TGA, TAG, TAA.")
    }
    if (!identical(seqlist[["transcript"]], test_seq)) {
      message("seqlist and test unequal.")
    }
    if (!identical(df_seq, test_seq)) {
      message("df_seq and test unequal.")
    }
  }



  return(return_list)

}



add_introns <- function(gtf_df, strand, n_exon) {
  gtf_df <- gtf_df[order(gtf_df[["start"]], decreasing = ifelse(strand == "+", F, T)),]
  gtf_df_exon <- gtf_df[which(gtf_df$feature == "exon"),]
  # for minus strand: exons in reverse to infer introns
  gtf_df_exon <- gtf_df_exon[order(gtf_df_exon[["start"]]),]
  # add introns to gtf_df
  pos_gen_intron <- stats::setNames(igsc:::seq2(gtf_df_exon$end[-length(gtf_df_exon$end)] + 1 , gtf_df_exon$start[-1] - 1),
                                    paste0("intron", sprintf("%02d", gtf_df_exon$exon_number[-ifelse(strand == "+", n_exon, 1)])))
  intron_ranges <- lapply(pos_gen_intron, range)
  gtf_df_intron <- data.frame(feature = "intron",
                              # min,max: same for minus and plus strand
                              start = unname(sapply(intron_ranges, min)),
                              end = unname(sapply(intron_ranges, max)),
                              strand = strand,
                              exon_number = as.numeric(substr(names(intron_ranges), 7, 9)))
  gtf_df <-
    dplyr::bind_rows(gtf_df, gtf_df_intron) %>%
    dplyr::select(-index) %>%
    dplyr::mutate(dplyr::across(-exon_number, ~ tidyr::replace_na(., get_mode(.))))
  #tidyr::fill(seqname, source, score, strand, gene_id, gene_name)
  return(gtf_df)
}

make_pos_gen_stack <- function(pos_gen) {
  pos_gen_stack <- list()
  j <- 0
  for (i in names(pos_gen)) {
    j <- j + 1
    pos_gen_stack[[j]] <- utils::stack(pos_gen[i])
    names(pos_gen_stack[[j]]) <- c("position_genome", "pattern")
    pos_gen_stack[[j]]$pattern <- as.character(pos_gen_stack[[j]]$pattern)
  }
  names(pos_gen_stack) <- names(pos_gen)
  return(pos_gen_stack)
}

get_mode <- function(x) {
  x <- na.omit(x) # Remove NA values
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}



