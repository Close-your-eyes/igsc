#' Title
#'
#' @param gtf_file
#' @param genome_file
#' @param bam_file
#' @param seqname
#' @param gene
#' @param n_reads_sample
#'
#' @return
#' @export
#'
#' @examples
plot_sc_reads <- function(gtf_file,
                          genome_file,
                          bam_file,
                          seqname = "chr2",
                          gene = "CD8A",
                          n_reads_sample = 40) {
  gene_info <-
    read_gtf(file_path = gtf_file, seqname_filter =  seqname)[["gtf"]] %>%
    dplyr::filter(gene_name == gene) %>%
    dplyr::filter(feature == "gene")

  seq_bounds <- get_fasta_seq_bounds(genome_file)
  refseq <- read_fasta(genome_file,
                       start_line = seq_bounds[which(seq_bounds$name == seqname), "start_line"],
                       end_line = seq_bounds[which(seq_bounds$name == seqname), "end_line"])

  refseq_range <- GenomicRanges::GRanges(seqnames = seqname, strand = gene_info$strand[1],
                                         ranges = IRanges::IRanges(start = gene_info$start[1], end = gene_info$end[1]))
  reads <- scexpr::reads_from_bam(bam_file, genomic_ranges = refseq_range, revcomp_minus_strand = F)
  reads_sub <- dplyr::slice_sample(reads, n = n_reads_sample)

  pattern_df <- purrr::pmap_dfr(list(reads_sub$cigar,
                                     reads_sub$pos,
                                     reads_sub$seq,
                                     reads_sub$qname),
                                function(x,y,z,a) cigar_to_position(cigar = x, start = y, seq = z, name = a))
  gene_refseq <- stringi::stri_sub(refseq,
                                   from = gene_info$start[1],
                                   to = gene_info$end[1])
  algnmt_df <- data.frame(seq = strsplit(gene_refseq, "")[[1]],
                          position = gene_info$start[1]:gene_info$end[1],
                          seq.name = gene)
  algnmt_df <- dplyr::bind_rows(algnmt_df, pattern_df)
  #algnmt_df <- dplyr::filter(algnmt_df, dplyr::between(position, min(reads_sub$pos)-20, min(reads_sub$pos)+200))
  plot <- algnmt_plot(algnmt = algnmt_df, algnmt_type = "NT", ref = gene)

  return(list(plot = plot, algnmt_df = algnmt_df, reads = reads, gene_info = gene_info, gene_refseq = gene_refseq))
}
