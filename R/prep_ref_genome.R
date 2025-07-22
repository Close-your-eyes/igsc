#' Title
#'
#' @param path
#' @param genome_fasta
#' @param genes_gtf
#'
#' @return
#' @export
#'
#' @examples
prep_ref_genome <- function(path,
                            genome_fasta = "genome.fa",
                            genes_gtf = "genes.gtf") {

  genes_gtf_to_fst(gtf_file = file.path(path, genes_gtf))
  genome_fasta_to_fst(genome_file = file.path(path, genome_fasta))

}

