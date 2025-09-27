#' Title
#'
#' @param gtf_path
#' @param genome_path
#' @param gene_name
#'
#' @returns
#' @export
#'
#' @examples
#'\dontrun{
#' out <- get_gene_sequences(gtf_path = "/Users/chris/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genes.gtf",
#'                           genome_path = "/Users/chris/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genome.fa",
#'                           gene_name = "GZMB")
#' }
get_gene_sequences <- function(gtf_path,
                               genome_path,
                               gene_name) {
  exondf <- read_gtf(file_path = gtf_path,
                     gene_names = gene_name,
                     use_fun = "r")[["gtf"]]

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
