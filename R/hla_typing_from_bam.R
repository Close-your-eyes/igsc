#' Extract candidate HLA reads from an RNA-seq BAM file
#'
#' Extract alignments from the chromosome 6 major histocompatibility complex
#' (MHC), together with chromosome 6 alternative or HLA-named contigs. This
#' provides a compact set of candidate HLA reads suitable for [hla_typing()].
#'
#' The default interval, 28--34 Mb on chromosome 6, deliberately extends beyond
#' the classical HLA genes. It works for both GRCh37 and GRCh38 and avoids
#' relying on exact exon coordinates, which can lose informative reads near
#' splice junctions or reads placed at a related HLA locus. Alternatively,
#' `reference_gtf` can be supplied to query only the annotated boundaries of
#' the genes in `genes`. Explicit `regions` take precedence over both methods.
#'
#' When `reference_genome` is supplied it must be an indexed FASTA file. Its
#' sequence names and lengths are checked against the BAM header. A FASTA is
#' not otherwise required because a BAM already records the reference sequence
#' names and lengths used for alignment.
#'
#' Only alignments represented in the BAM can be recovered. Unmapped reads and
#' reads aligned outside the selected regions are not classified as HLA reads by
#' this function.
#'
#' @param bam Path to a coordinate-sorted BAM file with a corresponding index.
#' @param reference_genome Optional path to the FASTA file used for alignment.
#'   The FASTA must have a `.fai` index and is used to validate reference names
#'   and lengths.
#' @param reference_gtf Optional path to a reference GTF file. When `regions`
#'   is `NULL`, gene records are read with [read_gtf()] and the annotated
#'   boundaries of `genes` are used as the query ranges.
#' @param regions Optional [GenomicRanges::GRanges()] object. Explicit regions
#'   take precedence over `reference_gtf`. When both are `NULL`, the chromosome
#'   6 interval specified by `mhc_start` and `mhc_end` is used and, when
#'   `include_alt_contigs = TRUE`, chromosome 6 alternative and HLA-named
#'   contigs are added.
#' @param genes Character vector of HLA genes whose GTF boundaries should be
#'   queried, for example `c("A", "B", "C")` or `c("HLA-A", "HLA-B")`.
#'   This argument is used only when `reference_gtf` is supplied and `regions`
#'   is `NULL`.
#' @param mhc_start,mhc_end One-based inclusive bounds of the broad MHC interval
#'   used when `regions` is `NULL`.
#' @param include_alt_contigs Logical; include whole BAM reference sequences
#'   whose names identify them as chromosome 6 alternative or HLA contigs.
#' @param tags BAM tags passed to [get_bam_reads()]. The default includes common
#'   Cell Ranger cell-barcode, UMI, and alignment tags.
#' @param scores Logical; calculate base-quality summaries in
#'   [get_bam_reads()].
#' @param min_mapq Minimum mapping quality to retain.
#' @param primary_only Logical; remove secondary and supplementary alignments.
#' @param cell_barcodes Optional character vector of corrected cell barcodes to
#'   retain. This requires the tag named by `cell_barcode_tag`.
#' @param cell_barcode_tag Name of the BAM tag containing corrected cell
#'   barcodes.
#'
#' @return A data frame in the format returned by [get_bam_reads()], with a
#'   unique `readName` column for direct use by [hla_typing()]. The queried
#'   ranges are stored in the `hla_regions` attribute.
#' @export
#'
#' @examples
#' \dontrun{
#' reads_hla <- get_hla_reads(
#'   bam = "outs/possorted_genome_bam.bam",
#'   reference_genome = "reference/fasta/genome.fa",
#'   reference_gtf = "reference/genes/genes.gtf",
#'   genes = "HLA-A"
#' )
#'
#' type_a <- hla_typing(
#'   hla_ref = subset(hla_ref, gene == "A"),
#'   reads = reads_hla,
#'   read_name_col_name = "readName"
#' )
#' }
get_hla_reads <- function(
    bam,
    reference_genome = NULL,
    reference_gtf = NULL,
    regions = NULL,
    genes = c("A", "B", "C"),
    mhc_start = 28000000L,
    mhc_end = 34000000L,
    include_alt_contigs = TRUE,
    tags = c("CR", "CB", "CY", "AS", "UR", "UB", "UY", "HI", "NH",
             "nM", "RE"),
    scores = TRUE,
    min_mapq = 0L,
    primary_only = TRUE,
    cell_barcodes = NULL,
    cell_barcode_tag = "CB") {

  igsc:::.ensure_packages(c("GenomicRanges", "IRanges", "Rsamtools"))

  if (missing(bam) || length(bam) != 1L || is.na(bam) || !nzchar(bam) ||
      !file.exists(bam)) {
    stop("`bam` must be the path to an existing BAM file.", call. = FALSE)
  }
  if (!is.null(reference_genome) &&
      (length(reference_genome) != 1L || is.na(reference_genome) ||
       !nzchar(reference_genome) || !file.exists(reference_genome))) {
    stop("`reference_genome` must be NULL or an existing FASTA file.",
         call. = FALSE)
  }
  if (is.null(regions) && !is.null(reference_gtf) &&
      (length(reference_gtf) != 1L || is.na(reference_gtf) ||
       !nzchar(reference_gtf) || !file.exists(reference_gtf))) {
    stop("`reference_gtf` must be NULL or an existing GTF file.",
         call. = FALSE)
  }
  if (!is.logical(include_alt_contigs) || length(include_alt_contigs) != 1L ||
      is.na(include_alt_contigs)) {
    stop("`include_alt_contigs` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(scores) || length(scores) != 1L || is.na(scores)) {
    stop("`scores` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(primary_only) || length(primary_only) != 1L ||
      is.na(primary_only)) {
    stop("`primary_only` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(min_mapq) != 1L || !is.numeric(min_mapq) || is.na(min_mapq) ||
      !is.finite(min_mapq) || min_mapq < 0 || min_mapq > 255 ||
      min_mapq != floor(min_mapq)) {
    stop("`min_mapq` must be one integer from 0 through 255.", call. = FALSE)
  }
  if (!is.character(tags) || anyNA(tags)) {
    stop("`tags` must be a character vector without missing values.",
         call. = FALSE)
  }
  if (length(cell_barcode_tag) != 1L || !is.character(cell_barcode_tag) ||
      is.na(cell_barcode_tag) || !nzchar(cell_barcode_tag)) {
    stop("`cell_barcode_tag` must be a non-empty character scalar.",
         call. = FALSE)
  }

  bam_stats <- Rsamtools::idxstatsBam(bam)
  bam_stats$seqnames <- as.character(bam_stats$seqnames)
  bam_stats <- bam_stats[!is.na(bam_stats$seqlength) &
                           bam_stats$seqlength > 0L, , drop = FALSE]
  if (nrow(bam_stats) == 0L) {
    stop("The BAM index contains no reference sequences.", call. = FALSE)
  }

  if (!is.null(reference_genome)) {
    validate_hla_reference_genome(reference_genome, bam_stats)
  }

  if (is.null(regions) && !is.null(reference_gtf)) {
    regions <- make_hla_gtf_regions(
      reference_gtf = reference_gtf,
      genes = genes,
      bam_stats = bam_stats
    )
  } else if (is.null(regions)) {
    regions <- make_hla_bam_regions(
      bam_stats = bam_stats,
      mhc_start = mhc_start,
      mhc_end = mhc_end,
      include_alt_contigs = include_alt_contigs
    )
  } else {
    if (!methods::is(regions, "GRanges") || length(regions) == 0L) {
      stop("`regions` must be a non-empty GRanges object.", call. = FALSE)
    }
    missing_seqnames <- setdiff(
      unique(as.character(GenomicRanges::seqnames(regions))),
      bam_stats$seqnames
    )
    if (length(missing_seqnames) > 0L) {
      stop(
        "These region sequence names are absent from the BAM: ",
        paste(missing_seqnames, collapse = ", "),
        call. = FALSE
      )
    }
  }

  message(
    "Extracting candidate HLA reads from ", length(regions), " region",
    if (length(regions) == 1L) "." else "s."
  )
  reads <- get_bam_reads(
    bam = bam,
    granges = regions,
    tags = unique(c(tags, if (!is.null(cell_barcodes)) cell_barcode_tag)),
    scores = scores,
    revcomp_minus_strand = FALSE,
    revcomp_plus_strand = FALSE
  )

  if (nrow(reads) == 0L) {
    reads$readName <- character()
    attr(reads, "hla_regions") <- regions
    return(reads)
  }

  if (primary_only && "flag" %in% names(reads)) {
    secondary_or_supplementary <-
      bitwAnd(as.integer(reads$flag), 256L) != 0L |
      bitwAnd(as.integer(reads$flag), 2048L) != 0L
    reads <- reads[!secondary_or_supplementary, , drop = FALSE]
  }
  if (min_mapq > 0L && "mapq" %in% names(reads)) {
    reads <- reads[!is.na(reads$mapq) & reads$mapq >= min_mapq, , drop = FALSE]
  }

  if (!is.null(cell_barcodes)) {
    if (!cell_barcode_tag %in% names(reads)) {
      stop(
        "The requested cell-barcode tag `", cell_barcode_tag,
        "` was not found in the BAM records.",
        call. = FALSE
      )
    }
    cell_barcodes <- unique(as.character(cell_barcodes))
    reads <- reads[
      !is.na(reads[[cell_barcode_tag]]) &
        reads[[cell_barcode_tag]] %in% cell_barcodes,
      ,
      drop = FALSE
    ]
  }

  # A read overlapping two supplied ranges can be returned twice by scanBam().
  alignment_key_columns <- intersect(
    c("qname", "flag", "rname", "pos", "cigar", "seq"),
    names(reads)
  )
  if (length(alignment_key_columns) > 0L) {
    reads <- reads[
      !duplicated(reads[alignment_key_columns]),
      ,
      drop = FALSE
    ]
  }

  qname <- as.character(reads$qname)
  if ("flag" %in% names(reads)) {
    read_end <- ifelse(
      bitwAnd(as.integer(reads$flag), 64L) != 0L,
      "/1",
      ifelse(bitwAnd(as.integer(reads$flag), 128L) != 0L, "/2", "")
    )
    qname <- paste0(qname, read_end)
  }
  reads$readName <- make.unique(qname)
  rownames(reads) <- NULL
  attr(reads, "hla_regions") <- regions
  reads
}


#' Type HLA genes directly from an RNA-seq BAM file
#'
#' A convenience workflow that first calls [get_hla_reads()] and then calls
#' [hla_typing()] separately for each requested HLA gene.
#'
#' @inheritParams get_hla_reads
#' @param hla_ref HLA allele reference data frame, preferably created by
#'   [hla_df_from_xml()].
#' @param genes Character vector of genes to type, for example `c("A", "B",
#'   "C")`. An optional `HLA-` prefix is ignored. Use `NULL` to type every gene
#'   represented in `hla_ref`. When `reference_gtf` is supplied, the same genes
#'   determine which annotated gene boundaries are queried.
#' @param hla_gene_col_name Column in `hla_ref` containing HLA gene names.
#' @param allele_diff Maximum fold difference used to retain candidate alleles;
#'   passed to [hla_typing()].
#' @param top_n_pairwise_results Number of leading allele pairs to plot; passed
#'   to [hla_typing()].
#' @param hla_seq_col_name Name of the HLA sequence column in `hla_ref`.
#' @param hla_allele_col_name Name of the allele identifier column in
#'   `hla_ref`.
#' @param p_group_col_name Name of the P-group column in `hla_ref`.
#' @param g_group_col_name Name of the G-group column in `hla_ref`.
#' @param lapply_fun Apply function used for allele matching; passed to
#'   [hla_typing()].
#' @param maxmis Maximum mismatches allowed per read match; passed to
#'   [hla_typing()].
#' @param make_reads_distinct Logical; remove duplicate read sequences before
#'   matching; passed to [hla_typing()].
#' @param ... Additional arguments passed to `lapply_fun` by [hla_typing()],
#'   such as `mc.cores` for [parallel::mclapply()].
#'
#' @return A list with `reads`, the candidate HLA read data frame; `typing`, a
#'   named list of [hla_typing()] results; and `regions`, the queried genomic
#'   ranges.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- hla_typing_from_bam(
#'   bam = "outs/possorted_genome_bam.bam",
#'   hla_ref = hla_ref,
#'   reference_gtf = "reference/genes/genes.gtf",
#'   genes = c("A", "B", "C"),
#'   maxmis = 1)
#'
#' # check strandness of HLA genes
#' gtf <- read_gtf(".../refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
#'                 gene_names = "HLA-",
#'                 gene_names_full_match = F,
#'                 features = "gene")$gtf |>
#'   dplyr::select(gene_name, start, end, strand) |>
#'   dplyr::filter(!grepl("AS", gene_name)) |>
#'   dplyr::arrange(gene_name)
#' # genes on (-) Strand: revcomp the reference as reads from BAM file are always
#' # mapped to (+) Strand but reference give sequence from (-) Strand
#' #   |gene_name |    start|      end|strand |
#' #   |:---------|--------:|--------:|:------|
#' #   |HLA-A     | 29941260| 29945884|+      |
#' #   |HLA-B     | 31269491| 31357188|-      |
#' #   |HLA-C     | 31268749| 31272130|-      |
#' #   |HLA-DMA   | 32948613| 32969094|-      |
#' #   |HLA-DMB   | 32934629| 32941028|-      |
#' #   |HLA-DOA   | 33004182| 33009591|-      |
#' #   |HLA-DOB   | 32812763| 32820466|-      |
#' #   |HLA-DPA1  | 33064569| 33080775|-      |
#' #   |HLA-DPB1  | 33075990| 33089696|+      |
#' #   |HLA-DQA1  | 32628179| 32647062|+      |
#' #   |HLA-DQA2  | 32741391| 32747198|+      |
#' #   |HLA-DQB1  | 32659467| 32668383|-      |
#' #   |HLA-DQB2  | 32756098| 32763532|-      |
#' #   |HLA-DRA   | 32439878| 32445046|+      |
#' #   |HLA-DRB1  | 32578769| 32589848|-      |
#' #   |HLA-DRB5  | 32517353| 32530287|-      |
#' #   |HLA-E     | 30489509| 30494194|+      |
#' #   |HLA-F     | 29722775| 29738528|+      |
#' #   |HLA-G     | 29826967| 29831125|+      |
#'
#'
#' library(igsc)
#' hla <- hla_df_from_xml(
#'   "/Volumes/CMS_SSD_2TB/hla.xml.gz",
#'   lapply_fun = parallel::mclapply,
#'   mc.cores = 8)
#'
#' hlaA <- hla |>
#'   dplyr::filter(grepl("^HLA-A", allele)) |>
#'   dplyr::filter(seq_length == max(seq_length), .by = allele_protein)
#'
#' hlaresA <- hla_typing_from_bam(
#'   bam = "/Users/chris/Documents/possorted_genome_bam.bam",
#'   reference_gtf = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
#'   hla_ref = hlaA,
#'   genes = "A",
#'   lapply_fun = parallel::mclapply,
#'   make_reads_distinct = T,
#'   mc.cores = 4,
#'   maxmis = 3,
#'   allele_diff = 4)
#' hlaresA$typing$B$plot_pair_res2
#'
#' ## HLA-B and HLA-C are on (-) Strand but all bam reads are returned as on (+) Strand
#' ## hence revcomp the reference: (-) --> (+)
#' hlaB <- hla |>
#'   dplyr::filter(grepl("^HLA-B", allele)) |>
#'   dplyr::filter(seq_length == max(seq_length), .by = allele_protein) |>
#'   dplyr::mutate(seq_Exon2_3 = revcompDNA(seq_Exon2_3))
#'
#' hlaresB <- hla_typing_from_bam(
#'   bam = "/Users/chris/Documents/possorted_genome_bam.bam",
#'   reference_gtf = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
#'   hla_ref = hlaB,
#'   genes = "B",
#'   lapply_fun = parallel::mclapply,
#'   make_reads_distinct = T,
#'   mc.cores = 4,
#'   maxmis = 3,
#'   allele_diff = 4)
#' hlaresB$typing$B$plot_pair_res2
#'
#' hlaC <- hla |>
#'   dplyr::filter(grepl("^HLA-C", allele)) |>
#'   dplyr::filter(seq_length == max(seq_length), .by = allele_protein) |>
#'   dplyr::mutate(seq_Exon2_3 = revcompDNA(seq_Exon2_3))
#'
#' hlaresC <- hla_typing_from_bam(
#'   bam = "/Users/chris/Documents/possorted_genome_bam.bam",
#'   reference_gtf = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
#'   hla_ref = hlaC,
#'   genes = "C",
#'   lapply_fun = parallel::mclapply,
#'   make_reads_distinct = T,
#'   mc.cores = 4,
#'   maxmis = 3,
#'   allele_diff = 4)
#' hlaresC$typing$C$plot_pair_res2
#'
#' hlaDPA1 <- hla |>
#'   dplyr::filter(grepl("^HLA-DPA1", allele)) |>
#'   dplyr::filter(seq_length == max(seq_length), .by = allele_protein) |>
#'   dplyr::mutate(seq_Exon2_3 = revcompDNA(seq_Exon2_3))
#'
#' hlaresDPA1 <- hla_typing_from_bam(
#'   bam = "/Users/chris/Documents/possorted_genome_bam.bam",
#'   reference_gtf = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
#'   hla_ref = hlaDPA1,
#'   genes = "DPA1",
#'   lapply_fun = parallel::mclapply,
#'   make_reads_distinct = T,
#'   mc.cores = 4,
#'   maxmis = 3,
#'   allele_diff = 2)
#' hlaresDPA1$typing$DPA1$plot_pair_res2
#' }
hla_typing_from_bam <- function(
    bam,
    hla_ref,
    reference_genome = NULL,
    reference_gtf = NULL,
    regions = NULL,
    genes = c("A", "B", "C"),
    hla_gene_col_name = "gene",
    mhc_start = 28000000L,
    mhc_end = 34000000L,
    include_alt_contigs = TRUE,
    tags = c("CR", "CB", "CY", "AS", "UR", "UB", "UY", "HI", "NH",
             "nM", "RE"),
    scores = TRUE,
    min_mapq = 0L,
    primary_only = TRUE,
    cell_barcodes = NULL,
    cell_barcode_tag = "CB",
    allele_diff = 5,
    top_n_pairwise_results = 50,
    hla_seq_col_name = "seq_Exon2_3",
    hla_allele_col_name = "allele",
    p_group_col_name = "p_group",
    g_group_col_name = "g_group",
    lapply_fun = lapply,
    maxmis = 3,
    make_reads_distinct = FALSE,
    ...) {

  if (!is.data.frame(hla_ref) || nrow(hla_ref) == 0L) {
    stop("`hla_ref` must be a non-empty data frame.", call. = FALSE)
  }
  if (length(hla_gene_col_name) != 1L || !is.character(hla_gene_col_name) ||
      is.na(hla_gene_col_name) || !nzchar(hla_gene_col_name) ||
      !hla_gene_col_name %in% names(hla_ref)) {
    stop("`hla_gene_col_name` must name a column in `hla_ref`.",
         call. = FALSE)
  }

  ref_genes <- normalize_hla_gene(hla_ref[[hla_gene_col_name]])
  if (is.null(genes)) {
    genes <- unique(ref_genes[!is.na(ref_genes) & nzchar(ref_genes)])
  } else {
    if (!is.character(genes) || length(genes) == 0L || anyNA(genes)) {
      stop("`genes` must be NULL or a non-empty character vector.",
           call. = FALSE)
    }
    genes <- unique(normalize_hla_gene(genes))
  }

  missing_genes <- setdiff(genes, unique(ref_genes))
  if (length(missing_genes) > 0L) {
    stop(
      "No reference alleles were found for: ",
      paste(paste0("HLA-", missing_genes), collapse = ", "),
      call. = FALSE
    )
  }

  reads <- get_hla_reads(
    bam = bam,
    reference_genome = reference_genome,
    reference_gtf = reference_gtf,
    regions = regions,
    genes = genes,
    mhc_start = mhc_start,
    mhc_end = mhc_end,
    include_alt_contigs = include_alt_contigs,
    tags = tags,
    scores = scores,
    min_mapq = min_mapq,
    primary_only = primary_only,
    cell_barcodes = cell_barcodes,
    cell_barcode_tag = cell_barcode_tag
  )
  if (nrow(reads) == 0L) {
    warning("No candidate HLA reads were found in the selected regions.",
            call. = FALSE)
    return(list(reads = reads, typing = stats::setNames(vector("list", length(genes)), genes),
                regions = attr(reads, "hla_regions")))
  }

  typing <- stats::setNames(vector("list", length(genes)), genes)
  for (gene in genes) {
    message("Typing HLA-", gene, ".")
    gene_ref <- hla_ref[ref_genes == gene & !is.na(ref_genes), , drop = FALSE]
    typing[[gene]] <- hla_typing(
      hla_ref = gene_ref,
      reads = reads,
      allele_diff = allele_diff,
      top_n_pairwise_results = top_n_pairwise_results,
      hla_seq_col_name = hla_seq_col_name,
      read_seq_col_name = "seq",
      hla_allele_col_name = hla_allele_col_name,
      read_name_col_name = "readName",
      p_group_col_name = p_group_col_name,
      g_group_col_name = g_group_col_name,
      lapply_fun = lapply_fun,
      maxmis = maxmis,
      make_reads_distinct = make_reads_distinct,
      rev_comp_minus = FALSE,
      strand_col_name = "strand",
      ...
    )
  }

  list(
    reads = reads,
    typing = typing,
    regions = attr(reads, "hla_regions")
  )
}


normalize_hla_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("^HLA-", "", x)
}


make_hla_gtf_regions <- function(reference_gtf, genes, bam_stats) {
  if (!is.character(genes) || length(genes) == 0L || anyNA(genes) ||
      any(!nzchar(trimws(genes)))) {
    stop(
      "`genes` must be a non-empty character vector when `reference_gtf` is used.",
      call. = FALSE
    )
  }

  genes <- unique(paste0("HLA-", normalize_hla_gene(genes)))
  gtf <- tryCatch(
    read_gtf(
      file_path = reference_gtf,
      features = "gene",
      gene_names = genes
    )[["gtf"]],
    error = function(e) {
      stop(
        "Could not read HLA gene boundaries from `reference_gtf`: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  required_columns <- c("seqname", "start", "end", "strand", "gene_name")
  missing_columns <- setdiff(required_columns, names(gtf))
  if (length(missing_columns) > 0L) {
    stop(
      "The processed GTF is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  gene_name <- toupper(trimws(as.character(gtf$gene_name)))
  keep <- !is.na(gene_name) & gene_name %in% genes
  gtf <- gtf[keep, , drop = FALSE]
  gene_name <- gene_name[keep]

  missing_genes <- setdiff(genes, unique(gene_name))
  if (length(missing_genes) > 0L) {
    stop(
      "No gene boundary was found in `reference_gtf` for: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  seqname <- as.character(gtf$seqname)
  start <- suppressWarnings(as.integer(gtf$start))
  end <- suppressWarnings(as.integer(gtf$end))
  strand <- as.character(gtf$strand)
  valid <- !is.na(seqname) & nzchar(seqname) & !is.na(start) & !is.na(end) &
    start >= 1L & end >= start
  if (!all(valid)) {
    stop("The selected GTF gene records contain invalid genomic bounds.",
         call. = FALSE)
  }

  # A GTF can contain repeated gene records. Collapse records for the same gene,
  # sequence and strand while retaining separate ranges on alternate contigs.
  key <- paste(gene_name, seqname, strand, sep = "\r")
  groups <- split(seq_along(key), key)
  bounds <- lapply(groups, function(i) {
    data.frame(
      gene_name = gene_name[i[[1L]]],
      seqname = seqname[i[[1L]]],
      strand = strand[i[[1L]]],
      start = min(start[i]),
      end = max(end[i]),
      stringsAsFactors = FALSE
    )
  })
  bounds <- do.call(rbind, bounds)
  rownames(bounds) <- NULL

  missing_seqnames <- setdiff(unique(bounds$seqname), bam_stats$seqnames)
  if (length(missing_seqnames) > 0L) {
    stop(
      "These GTF sequence names are absent from the BAM: ",
      paste(missing_seqnames, collapse = ", "),
      ". Use a GTF from the same reference assembly as the BAM.",
      call. = FALSE
    )
  }

  bam_lengths <- stats::setNames(
    as.numeric(bam_stats$seqlength),
    as.character(bam_stats$seqnames)
  )
  beyond_bam <- bounds$end > bam_lengths[bounds$seqname]
  if (any(beyond_bam)) {
    stop(
      "GTF gene boundaries extend beyond the BAM reference length for: ",
      paste(unique(bounds$seqname[beyond_bam]), collapse = ", "),
      call. = FALSE
    )
  }

  GenomicRanges::GRanges(
    seqnames = bounds$seqname,
    ranges = IRanges::IRanges(start = bounds$start, end = bounds$end),
    strand = bounds$strand,
    gene_name = bounds$gene_name
  )
}


make_hla_bam_regions <- function(bam_stats,
                                 mhc_start,
                                 mhc_end,
                                 include_alt_contigs) {
  if (length(mhc_start) != 1L || length(mhc_end) != 1L ||
      !is.numeric(mhc_start) || !is.numeric(mhc_end) ||
      is.na(mhc_start) || is.na(mhc_end) ||
      !is.finite(mhc_start) || !is.finite(mhc_end) ||
      mhc_start < 1L || mhc_end < mhc_start ||
      mhc_start != floor(mhc_start) || mhc_end != floor(mhc_end)) {
    stop(
      "`mhc_start` and `mhc_end` must be valid one-based integer bounds.",
      call. = FALSE
    )
  }

  seqnames <- as.character(bam_stats$seqnames)
  lengths <- as.numeric(bam_stats$seqlength)
  primary <- grepl("^(chr)?6$", seqnames, ignore.case = TRUE) |
    grepl("^NC_000006(\\.[0-9]+)?$", seqnames, ignore.case = TRUE)

  if (!any(primary)) {
    known_chr6_lengths <- c(170805979, 170805980, 171115067, 171115068)
    primary <- lengths %in% known_chr6_lengths
  }
  if (sum(primary) != 1L) {
    stop(
      "Could not identify one primary chromosome 6 sequence in the BAM. ",
      "Supply `regions` explicitly.",
      call. = FALSE
    )
  }

  primary_name <- seqnames[primary]
  primary_length <- lengths[primary]
  if (mhc_start > primary_length) {
    stop("`mhc_start` lies beyond the chromosome 6 sequence length.",
         call. = FALSE)
  }

  region_names <- primary_name
  starts <- as.integer(mhc_start)
  ends <- as.integer(min(mhc_end, primary_length))

  if (include_alt_contigs) {
    alt <- !primary & (
      grepl("HLA", seqnames, ignore.case = TRUE) |
      grepl("^(chr)?6[_-].*(alt|hap|patch|fix)", seqnames, ignore.case = TRUE) |
      grepl("^(chr)?6_", seqnames, ignore.case = TRUE)
    )
    if (any(alt)) {
      region_names <- c(region_names, seqnames[alt])
      starts <- c(starts, rep.int(1L, sum(alt)))
      ends <- c(ends, as.integer(lengths[alt]))
    }
  }

  GenomicRanges::GRanges(
    seqnames = region_names,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )
}


validate_hla_reference_genome <- function(reference_genome, bam_stats) {
  fai <- paste0(reference_genome, ".fai")
  if (!file.exists(fai)) {
    stop(
      "The reference FASTA must be indexed; file not found: ", fai,
      call. = FALSE
    )
  }

  fasta_index <- Rsamtools::scanFaIndex(reference_genome)
  fasta_names <- as.character(GenomicRanges::seqnames(fasta_index))
  common <- intersect(as.character(bam_stats$seqnames), fasta_names)
  if (length(common) == 0L) {
    stop("The reference FASTA and BAM have no sequence names in common.",
         call. = FALSE)
  }

  bam_length <- stats::setNames(
    as.numeric(bam_stats$seqlength),
    as.character(bam_stats$seqnames)
  )
  fasta_length <- stats::setNames(
    as.numeric(IRanges::width(fasta_index)),
    fasta_names
  )
  mismatch <- common[bam_length[common] != fasta_length[common]]
  if (length(mismatch) > 0L) {
    stop(
      "The reference FASTA and BAM report different lengths for: ",
      paste(utils::head(mismatch, 5L), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
