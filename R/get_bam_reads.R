#' Read BAM alignments overlapping genomic ranges
#'
#' Extract alignments with [Rsamtools::scanBam()] and return one data-frame row
#' per alignment found in each requested range. The function does not filter by
#' mapping quality, alignment flag, or cell barcode.
#'
#' @details
#' A coordinate-sorted, indexed BAM is needed for range queries. The sequence
#' names in `granges` must match those in the BAM header (for example, `chr6`
#' and `6` are different names). If `granges` is omitted, only the first 10,000
#' bases of the BAM reference sequence with the most mapped reads are queried;
#' this is not a whole-BAM extraction. Range coordinates are one-based and
#' inclusive. The strand stored in `granges` does not restrict alignments.
#'
#' A read overlapping two input ranges can appear twice. Paired ends and
#' secondary or supplementary alignments can also share a `qname`. This function
#' does not deduplicate or create unique read names. For HLA-specific extraction
#' with additional filtering, see [get_hla_reads()]. Tight gene or exon ranges
#' may miss reads aligned nearby or to a related locus.
#'
#' The default tags include `CR`/`CB` (raw/corrected cell barcode), `UR`/`UB`
#' (raw/corrected UMI), `AS` (alignment score), and `NH` (number of reported
#' alignments), among others. Tag columns are present only when the BAM provides
#' them.
#'
#' BAM stores sequences aligned to the minus strand in reference orientation.
#' With both reverse-complement options `FALSE`, `seq` and `qual` remain in that
#' BAM orientation. `revcomp_minus_strand = TRUE` asks Rsamtools to return
#' minus-strand reads in sequenced orientation. `revcomp_plus_strand = TRUE`
#' additionally reverse-complements plus-strand read sequences and reverses
#' their quality strings. Setting both to `TRUE` puts returned sequences in the
#' orientation opposite to the reference plus strand, which can be useful when
#' comparing with a minus-strand reference sequence.
#'
#' @param bam Path to a coordinate-sorted BAM file with an index accessible to
#'   Rsamtools.
#' @param granges A [GenomicRanges::GRanges()] object specifying the ranges to
#'   query. Empty ranges are rejected. Its strand is ignored when selecting
#'   alignments. The default queries
#'   bases 1 through 10,000 of the reference sequence with the highest mapped
#'   read count in the BAM index.
#' @param tags Character vector of BAM tag names passed to
#'   [Rsamtools::ScanBamParam()]. Use `character(0)` to request no tags. The
#'   default requests common Cell Ranger barcode, UMI, and alignment tags;
#'   unavailable tags are omitted from the result and reported in a message.
#' @param scores Logical; if `TRUE`, calculate per-read Phred-quality summaries
#'   from `qual`.
#' @param revcomp_minus_strand Logical; passed as `reverseComplement` to
#'   [Rsamtools::ScanBamParam()]. If `TRUE`, reverse-complement `seq` and reverse
#'   `qual` for alignments on the minus strand.
#' @param revcomp_plus_strand Logical; if `TRUE`, reverse-complement `seq` and
#'   reverse `qual` for alignments on the plus strand after reading the BAM.
#' @param decode_bam_flags run decode_bam_flags and attach columns
#'
#' @return A data frame containing fields returned by
#'   [Rsamtools::scanBamWhat()], including `qname`, `flag`, `rname`, `strand`,
#'   `pos`, `mapq`, `cigar`, `seq`, and `qual`. `genomic_range` is the character
#'   form of the one-based queried-range index. Requested tags present in the BAM are
#'   added as columns. When `scores = TRUE`, four more columns are added:
#'   `readQualNum` (a list of per-base Phred scores), `minQual`, `meanQual`, and
#'   `n_belowQ30` (the number of bases with Phred score below 30). Array-valued
#'   BAM tags are stored as list columns. Missing qualities produce an empty
#'   integer vector in `readQualNum` and `NA` in the three quality summaries.
#' @seealso [get_hla_reads()]
#' @export
#'
#' @examples
#' \dontrun{
#' | value|meaning                                                    |
#' |-----:|:----------------------------------------------------------|
#' |     1|Paired-end or other multi-segment read                     |
#' |     2|Properly aligned pair, according to the aligner            |
#' |     4|Read is unmapped                                           |
#' |     8|Mate is unmapped                                           |
#' |    16|Read is reverse-complemented; minus strand if mapped       |
#' |    32|Mate is reverse-complemented; minus strand if mapped       |
#' |    64|First segment, usually R1                                  |
#' |   128|Last segment, usually R2                                   |
#' |   256|Secondary alignment: an alternative mapping                |
#' |   512|Failed platform/vendor quality checks                      |
#' |  1024|Marked as a PCR or optical duplicate                       |
#' |  2048|Supplementary alignment: another part of a split alignment |
#'
#' # see igsc::decode_bam_flags()
#'
#' bam <- "path/to/possorted_genome_bam.bam"
#' reference_names <- names(Rsamtools::scanBamHeader(bam)[[1]]$targets)
#' chr6 <- intersect(c("chr6", "6"), reference_names)[1]
#' stopifnot(!is.na(chr6))
#'
#' mhc <- GenomicRanges::GRanges(
#'   seqnames = chr6,
#'   ranges = IRanges::IRanges(start = 29000000L, end = 35000000L)
#' )
#'
#' mhc2 <- range <- GenomicRanges::GRanges(seqnames = "chr6",
#'                                         ranges = "28000000..34000000")
#'
#' reads <- get_bam_reads(bam, granges = mhc, tags = c("CB", "UB"))
#' reads <- reads[!is.na(reads$mapq) & reads$mapq >= 20, , drop = FALSE]
#' head(reads[, c("qname", "rname", "pos", "mapq", "seq")])
#' }
get_bam_reads <- function(bam,
                          granges = GenomicRanges::GRanges(seqnames = as.character(Rsamtools::idxstatsBam(bam)[which.max(Rsamtools::idxstatsBam(bam)[["mapped"]]),"seqnames"]),
                                                           ranges = "1..10000"),
                          tags = c("CR", "CB", "CY", "AS", "UR", "UB",
                                   "UY", "HI", "NH", "nM", "RE"),
                          scores = T,
                          revcomp_minus_strand = F,
                          revcomp_plus_strand = F,
                          decode_bam_flags = T) {

  igsc:::.ensure_packages(c("Biostrings", "GenomicRanges", "Rsamtools"))

  if (missing(bam) || bam == "" || !file.exists(bam)) {
    stop("bam not found or missing.")
  }

  if (!methods::is(granges, "GRanges") || length(granges) == 0L) {
    stop("`granges` must be a non-empty GRanges object.", call. = FALSE)
  }

  message("Reading BAM file.")
  params <- Rsamtools::ScanBamParam(
    which = granges,
    what = Rsamtools::scanBamWhat(),
    tag = tags,
    reverseComplement = revcomp_minus_strand
  )
  reads <- Rsamtools::scanBam(bam, param = params)

  # start of a read always refers to the (+)Strand, so for reads on the (-)Strand start is actually the end, (see IGV browser, read details)
  reads <- purrr::map_dfr(reads, function(x) {
    tag_values <- x$tag
    x$tag <- NULL
    no_tags <- names(tag_values)[vapply(tag_values, is.null, logical(1))]
    if (length(no_tags) > 0L) {
      message("These tags were not found: ", paste(no_tags, collapse = ","), ".")
    }
    result <- data.frame(x)
    for (tag_name in names(tag_values)) {
      if (!is.null(tag_values[[tag_name]])) {
        result[[tag_name]] <- tag_values[[tag_name]]
      }
    }
    result
  }, .id = "grange")

  if (revcomp_plus_strand) {
    # start and CIGAR remain the same ?? Not sure. But in ScanBamParam start and CIGAR remain the same when reverseComplement = T.
    reads[which(reads$strand == "+"), "seq"] <- revcomp_dna(reads[which(reads$strand == "+"), "seq"])
    reads[which(reads$strand == "+"), "qual"] <- stringi::stri_reverse(reads[which(reads$strand == "+"), "qual"])
    #lapply(lapply(strsplit(reads[which(reads$strand == "+"), "qual"], ""), rev), paste, collapse = "")
  }

  if (anyDuplicated(reads$qname)) {
    message("Some read names (qname) are duplicated.")
  }

  if (scores) {
    message("Calculating read score.")
    pq <- methods::as(Biostrings::PhredQuality(reads$qual), "IntegerList")
    # Rsamtools represents absent BAM qualities with negative Phred sentinels.
    pq <- lapply(pq, function(x) {
      if (anyNA(x) || any(x < 0L)) integer() else x
    })
    stats <- igsc:::qual_stats_cpp(pq)
    for (i in names(stats)) {
      reads[[i]] <- stats[[i]]
    }
    # reads$readQualNum <- unlist(lapply_fun(pq, paste, collapse = ".", ...))
    # reads$minQual <- unlist(lapply_fun(pq, min, ...))
    # reads$meanQual <- unlist(lapply_fun(pq, mean, ...))
    # reads$n_belowQ30 <- unlist(lapply_fun(pq, function(x) sum(x < 30), ...))
  }

  reads[["seq_duplicate"]] <- as.numeric(duplicated(reads$seq))

  if (decode_bam_flags) {
    reads <- cbind(reads, decode_bam_flags(reads$flag))
  }

  return(reads)
}
