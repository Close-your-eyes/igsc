#' Plot reads from a bam file from scRNAseq against the genomic reference
#'
#' This function should help to memorize the workflow of plotting reads from a
#' selected gene against the genomic reference sequence. E.g., to inspect
#' their location or to see mismatches. It is assumed that Cellranger pipeline
#' was used withe gtf_file and genome_file as input and bam_file as output.
#' gtf and genome files are references downloaded from ensemble, for instance.
#'
#' @param gtf_file path to the genes.gtf file
#' @param genome_file path to the genome.fa file
#' @param bam_file path to the possorted_genome_bam.bam file
#' @param seqname chromosome the gene is on, used to read this chromosome only
#' from genome_file; is passed to igsc::get_fasta_seq_bounds and then
#' to igsc::read_fasta; used to filter the gtf_file
#' @param n_reads_sample how many sampled reads to plot
#' @param gene_name used to filter from the gtf_file
#' @param feature used to filter from the gtf_file
#' @param ... arguments to igsc::algnmt_plot
#'
#' @return a list with a ggplot and intermediate results
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' result <- plot_sc_reads(gtf_file = "path1/genes.gtf",
#' genome_file = "path2/genome.fa",
#' bam_file = "path3/possorted_genome_bam.bam",
#' seqname = "chr2",
#' gene_name = "CD8A",
#' feature = "gene")
#' }
plot_sc_reads <- function(gtf_file,
                          genome_file,
                          bam_file,
                          seqname,
                          gene_name,
                          feature,
                          n_reads_sample = 40,
                          ...) {
  gene_info <-
    read_gtf(file_path = gtf_file,
             seqnames = seqname,
             features = feature)[["gtf"]] %>%
    dplyr::filter(gene_name == !!gene_name)
    #dplyr::filter(feature == !!feature)

  if (nrow(gene_info) == 0) {
    stop("no rows left in gtf data frame, check seqname, feature and gene_name.")
  }
  if (nrow(gene_info) != 1) {
    print(head(gene_info))
    stop("gtf data frame should have one row only after filtering. check.")
  }
  # gene info should have one line only

  seq_bounds <- get_fasta_seq_bounds(genome_file)
  if (!seqname %in% seq_bounds$name) {
    stop("seqname not found in genome file.")
  }
  refseq <- read_fasta(genome_file,
                       start_line = seq_bounds[which(seq_bounds$name == seqname), "start_line"],
                       end_line = seq_bounds[which(seq_bounds$name == seqname), "end_line"])

  refseq_range <- GenomicRanges::GRanges(seqnames = seqname, strand = gene_info$strand[1],
                                         ranges = IRanges::IRanges(start = gene_info$start[1], end = gene_info$end[1]))
  reads <- igsc::get_bam_reads(file_path = bam_file, genomic_ranges = refseq_range, revcomp_minus_strand = F)
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
                          seq.name = gene_name)
  algnmt_df <- dplyr::bind_rows(algnmt_df, pattern_df)
  #algnmt_df <- dplyr::filter(algnmt_df, dplyr::between(position, min(reads_sub$pos)-20, min(reads_sub$pos)+200))
  plot <- algnmt_plot(algnmt = algnmt_df, algnmt_type = "NT", ref = gene_name, subject_name = gene_name, ...)

  return(list(plot = plot, algnmt_df = algnmt_df, reads = reads, gene_info = gene_info, gene_refseq = gene_refseq))
}
