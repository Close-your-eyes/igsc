#' Get gene, transcript, amino acid sequences from biomaRt
#'
#' Retrieve canonical sequences. coding = cds (exons), cDNA = exons + UTR,
#' transcript_exon_intron = introns+exons, gene_exon_intron = full genomic span,
#' gene_exon = exons in genomic order, not revcomp
#'
#' @param hgnc_symbol gene symbol
#' @param dataset passed to biomaRt::useEnsembl
#'
#' @returns data frame
#' @export
#'
#' @examples
#'granzymes <- c("GZMA","GZMB","GZMH","GZMK","GZMM")
#'out <- get_sequences_from_biomart(granzymes)
get_sequences_from_biomart <- function(hgnc_symbol, dataset = "hsapiens_gene_ensembl") {

  mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # attr <- biomaRt::listAttributes(mart)

  df <- biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_peptide_id",
      "transcript_is_canonical",
      "strand"
    ),
    filters = "hgnc_symbol",
    values = hgnc_symbol,
    mart = mart
  ) |>
    dplyr::filter(transcript_is_canonical == 1)

  seq_types <- c(
    cds       = "coding", # CDS (coding DNA sequence)
    cdna      = "cdna", # cDNA (exons + UTR, no introns)
    # Transcript sequence with exons+introns (exons=UPPERCASE, introns=lowercase)
    # isoform-specific transcript sequence, exons/introns mapped exactly to that transcript.
    transcript_exon_intron = "transcript_exon_intron",
    gene_exon_intron = "gene_exon_intron",   # full genomic span of the gene, with all exons marked, regardless of transcript isoforms.
    peptide   = "peptide", # amino acid seq
    utr5      = "5utr",
    utr3      = "3utr",
    # genomic order
    # For a + strand gene, they’re listed left to right (smallest → largest coordinate).
    # For a – strand gene, they’re still listed in genomic order (smallest → largest coordinate), which means they’ll look “backwards” relative to transcription direction.
    gene_exon = "gene_exon"
  )

  sequences <- purrr::map(seq_types, ~ biomaRt::getSequence(
    id     = df$ensembl_transcript_id,
    type   = "ensembl_transcript_id",
    seqType= .x,
    mart   = mart
  ))
  sequences[["gene_exon"]] <- tidyr::nest(sequences[["gene_exon"]], gene_exon = gene_exon)

  df <- purrr::reduce(c(list(df = df), sequences), dplyr::left_join, by = "ensembl_transcript_id")

  return(df)
}

