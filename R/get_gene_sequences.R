#' Extract transcript sequences for a gene
#'
#' Read the annotations for one gene from a GTF file, split them by transcript,
#' and reconstruct sequence data for each transcript from a reference genome
#' FASTA file.
#'
#' The function locates the FASTA record named by the GTF `seqname`, reads that
#' record once, and passes each transcript's annotation rows and the reference
#' sequence to [concat_transcript()]. All selected transcript annotations are
#' therefore expected to refer to the same sequence record. Each individual
#' transcript must contain no more than one unique `seqname`.
#'
#' @param gtf_path Path to a GTF annotation file accepted by [read_gtf()]. The
#'   selected records must contain `transcript_id` and `seqname` columns.
#' @param genome_path Path to the reference genome FASTA file corresponding to
#'   `gtf_path`. FASTA record names must match the GTF `seqname` values and be
#'   discoverable by [get_fasta_seq_bounds()].
#' @param gene_name A single gene name used to filter the GTF annotations.
#'
#' @return Invisibly, a named list with one element per `transcript_id`. Each
#'   element is the sequence-data object returned by [concat_transcript()] for
#'   that transcript.
#' @export
#'
#' @examples
#' \dontrun{
#' transcripts <- get_gene_sequences(
#'   gtf_path = "reference/genes.gtf",
#'   genome_path = "reference/genome.fa",
#'   gene_name = "GZMB"
#' )
#'
#' names(transcripts)
#' transcripts[[1]]
#' }
get_gene_sequences <- function(gtf_path,
                               genome_path,
                               gene_name) {
  exondf <- read_gtf(file_path = gtf_path,
                     gene_names = gene_name)[["gtf"]]

  exonlist <- split(exondf, exondf$transcript_id)
  seq_bounds <- get_fasta_seq_bounds(genome_path)

  # all trancript versions for one gene will be on same seqname, so read it once outside of purrr::map
  refseq <- read_fasta(genome_path,
                       start_line = seq_bounds[which(seq_bounds$seqname == exonlist[[1]]$seqname[1]), "start_line"],
                       end_line = seq_bounds[which(seq_bounds$seqname == exonlist[[1]]$seqname[1]), "end_line"])

  out <- purrr::map(exonlist, function(x) {
    if (length(unique(x$seqname)) > 1) {
      stop("more than one seqname not allowed.")
    }
    # refseq <- read_fasta(genome_path,
    #                      start_line = seq_bounds[which(seq_bounds$seqname == x$seqname[1]), "start_line"],
    #                      end_line = seq_bounds[which(seq_bounds$seqname == x$seqname[1]), "end_line"])
    data <- concat_transcript(gtf_df = x,
                              refseq = refseq,
                              refseq_strand = "+",
                              run_test = T)
    return(data)
  })

}
