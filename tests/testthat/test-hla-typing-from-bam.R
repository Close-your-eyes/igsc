make_hla_test_bam <- function(directory, barcode_tag = "CB", chr6 = "chr6") {
  sam <- file.path(directory, "hla-test.sam")
  quality <- paste(rep("I", 20L), collapse = "")
  record <- function(qname, flag, rname, pos, mapq, barcode, base = "A") {
    sequence <- paste(rep(base, 20L), collapse = "")
    paste(
      qname, flag, rname, pos, mapq, "20M", "*", 0, 0,
      sequence, quality, paste0(barcode_tag, ":Z:", barcode),
      sep = "\t"
    )
  }

  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      "@SQ\tSN:chr1\tLN:1000",
      paste0("@SQ\tSN:", chr6, "\tLN:1000"),
      "@SQ\tSN:chr6_GL000250v2_alt\tLN:1000",
      record("outside", 0, "chr1", 100, 60, "CELL-1"),
      record("hla_primary", 0, chr6, 100, 60, "CELL-1"),
      record("hla_secondary", 256, chr6, 150, 60, "CELL-1"),
      record("hla_low_mapq", 0, chr6, 200, 5, "CELL-2", base = "C"),
      record("hla_alt", 0, "chr6_GL000250v2_alt", 10, 60, "CELL-1")
    ),
    sam
  )

  destination <- file.path(directory, "hla-test")
  suppressMessages(
    Rsamtools::asBam(sam, destination, indexDestination = TRUE)
  )
  paste0(destination, ".bam")
}


make_hla_test_gtf <- function(directory) {
  gtf <- file.path(directory, "hla-test.gtf")
  writeLines(
    c(
      paste(
        "chr6", "test", "gene", 80, 130, ".", "+", ".",
        'gene_id "HLA-A"; gene_name "HLA-A";',
        sep = "\t"
      ),
      paste(
        "chr6", "test", "gene", 180, 230, ".", "-", ".",
        'gene_id "HLA-B"; gene_name "HLA-B";',
        sep = "\t"
      ),
      paste(
        "chr6", "test", "gene", 300, 350, ".", "+", ".",
        'gene_id "OTHER"; gene_name "OTHER";',
        sep = "\t"
      )
    ),
    gtf
  )
  gtf
}


test_that("default regions contain the MHC and chromosome 6 alt contigs", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  bam_stats <- data.frame(
    seqnames = c("chr1", "chr6", "chr6_GL000250v2_alt"),
    seqlength = c(248956422, 170805979, 4672374)
  )

  regions <- igsc:::make_hla_bam_regions(
    bam_stats,
    mhc_start = 28000000L,
    mhc_end = 34000000L,
    include_alt_contigs = TRUE
  )

  expect_setequal(
    as.character(GenomicRanges::seqnames(regions)),
    c("chr6", "chr6_GL000250v2_alt")
  )
  expect_identical(IRanges::start(regions), c(28000000L, 1L))
  expect_identical(IRanges::end(regions), c(34000000L, 4672374L))
})


test_that("candidate extraction filters alignments and includes alt contigs", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  regions <- GenomicRanges::GRanges(
    c("chr6", "chr6_GL000250v2_alt"),
    IRanges::IRanges(c(1, 1), c(1000, 1000))
  )

  reads <- suppressMessages(
    get_hla_reads(
      bam,
      regions = regions,
      scores = FALSE,
      min_mapq = 10
    )
  )

  expect_setequal(reads$qname, c("hla_primary", "hla_alt"))
  expect_true(all(c("qname", "seq", "readName", "CB") %in% names(reads)))
  expect_identical(anyDuplicated(reads$readName), 0L)
  expect_s4_class(attr(reads, "hla_regions"), "GRanges")
  expect_setequal(
    as.character(GenomicRanges::seqnames(attr(reads, "hla_regions"))),
    c("chr6", "chr6_GL000250v2_alt")
  )
})


test_that("cell barcodes and custom regions filter candidate reads", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  region <- GenomicRanges::GRanges(
    "chr6",
    IRanges::IRanges(1, 300)
  )

  reads <- suppressMessages(
    get_hla_reads(
      bam,
      regions = region,
      scores = FALSE,
      cell_barcodes = "CELL-2"
    )
  )

  expect_identical(reads$qname, "hla_low_mapq")
  expect_identical(reads$CB, "CELL-2")
})


test_that("reference GTF gene boundaries restrict candidate reads", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("vroom")
  skip_if_not_installed("Gmisc")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  gtf <- make_hla_test_gtf(directory)

  reads <- suppressMessages(
    get_hla_reads(
      bam,
      reference_gtf = gtf,
      genes = "A",
      scores = FALSE
    )
  )

  expect_identical(reads$qname, "hla_primary")
  regions <- attr(reads, "hla_regions")
  expect_s4_class(regions, "GRanges")
  expect_identical(as.character(GenomicRanges::seqnames(regions)), "chr6")
  expect_identical(IRanges::start(regions), 80L)
  expect_identical(IRanges::end(regions), 130L)
  expect_identical(as.character(regions$gene_name), "HLA-A")
})


test_that("reference GTF reports missing HLA genes", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("vroom")
  skip_if_not_installed("Gmisc")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  gtf <- make_hla_test_gtf(directory)

  expect_error(
    get_hla_reads(
      bam,
      reference_gtf = gtf,
      genes = "C",
      scores = FALSE
    ),
    "Could not read HLA gene boundaries.*No rows left in gtf"
  )
})


test_that("explicit regions take precedence over a reference GTF", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  region <- GenomicRanges::GRanges(
    "chr6_GL000250v2_alt",
    IRanges::IRanges(1, 30)
  )

  reads <- suppressMessages(
    get_hla_reads(
      bam,
      reference_gtf = file.path(directory, "does-not-exist.gtf"),
      regions = region,
      genes = "A",
      scores = FALSE
    )
  )

  expect_identical(reads$qname, "hla_alt")
  expect_identical(attr(reads, "hla_regions"), region)
})


test_that("an indexed reference FASTA is validated against the BAM", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  fasta <- file.path(directory, "genome.fa")
  sequence <- paste(rep("A", 1000L), collapse = "")
  writeLines(
    c(
      ">chr1", sequence,
      ">chr6", sequence,
      ">chr6_GL000250v2_alt", sequence
    ),
    fasta
  )
  suppressMessages(Rsamtools::indexFa(fasta))
  region <- GenomicRanges::GRanges(
    "chr6",
    IRanges::IRanges(1, 300)
  )

  reads <- suppressMessages(
    get_hla_reads(
      bam,
      reference_genome = fasta,
      regions = region,
      scores = FALSE
    )
  )

  expect_setequal(reads$qname, c("hla_primary", "hla_low_mapq"))
})


test_that("BAM workflow extracts reads and types each requested gene", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("brathering")

  directory <- tempfile("hla-test-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  region <- GenomicRanges::GRanges(
    "chr6",
    IRanges::IRanges(1, 300)
  )
  sequence_a <- paste(rep("A", 20L), collapse = "")
  sequence_c <- paste(rep("C", 20L), collapse = "")
  hla_ref <- data.frame(
    gene = c("A", "A"),
    allele = c("HLA-A*01:01", "HLA-A*02:01"),
    seq_Exon2_3 = c(sequence_a, sequence_c),
    p_group = c("A*01:01P", "A*02:01P"),
    g_group = c("A*01:01G", "A*02:01G"),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(
    hla_typing_from_bam(
      bam = bam,
      hla_ref = hla_ref,
      regions = region,
      genes = "HLA-A",
      tags = character(),
      scores = FALSE,
      top_n_pairwise_results = 5
    )
  )

  expect_named(result, c("reads", "typing", "regions"))
  expect_named(result$typing, "A")
  expect_s3_class(result$typing$A$pair_res1_df, "data.frame")
  expect_gt(nrow(result$typing$A$pair_res1_df), 0L)
})


test_that("custom allele columns are honored when checking one-gene input", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("brathering")

  sequence <- paste(rep("A", 20L), collapse = "")
  hla_ref <- data.frame(
    allele_id = c("HLA-A*01:01", "HLA-A*02:01"),
    sequence = c(sequence, sequence),
    p = c("A*01:01P", "A*02:01P"),
    g = c("A*01:01G", "A*02:01G"),
    stringsAsFactors = FALSE
  )
  reads <- data.frame(
    read_id = "read-1",
    sequence = sequence,
    strand = "+",
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(
    hla_typing(
      hla_ref = hla_ref,
      reads = reads,
      hla_seq_col_name = "sequence",
      read_seq_col_name = "sequence",
      hla_allele_col_name = "allele_id",
      read_name_col_name = "read_id",
      p_group_col_name = "p",
      g_group_col_name = "g",
      top_n_pairwise_results = 5
    )
  )

  expect_s3_class(result$pair_res1_df, "data.frame")
})


test_that("multiple BAMs retain cross-file reads and source metadata", {
  skip_if_not_installed("Rsamtools")
  directories <- c(tempfile("hla-run1-"), tempfile("hla-run2-"))
  invisible(lapply(directories, dir.create))
  bam <- vapply(directories, make_hla_test_bam, character(1))
  region <- GenomicRanges::GRanges("chr6", IRanges::IRanges(1, 300))

  reads <- suppressMessages(get_hla_reads(
    bam, regions = region, scores = FALSE, min_mapq = 10
  ))

  expect_equal(nrow(reads), 2L)
  expect_identical(reads$qname, rep("hla_primary", 2L))
  expect_identical(reads$bam_file, unname(bam))
  expect_identical(reads$readName, c("1:hla_primary", "2:hla_primary"))
  expect_identical(anyDuplicated(reads$readName), 0L)
  regions <- attr(reads, "hla_regions")
  expect_s4_class(regions, "GRangesList")
  expect_identical(names(regions), unname(bam))
  expect_identical(regions[[1L]], region)
  expect_identical(regions[[2L]], region)
})


test_that("pooling fills unavailable BAM tags with missing values", {
  skip_if_not_installed("Rsamtools")
  directories <- c(tempfile("hla-run1-"), tempfile("hla-run2-"))
  invisible(lapply(directories, dir.create))
  bam <- c(make_hla_test_bam(directories[[1L]]),
           make_hla_test_bam(directories[[2L]], barcode_tag = "ZZ"))
  region <- GenomicRanges::GRanges("chr6", IRanges::IRanges(1, 130))

  reads <- suppressMessages(get_hla_reads(
    bam, regions = region, tags = c("CB", "ZZ"), scores = FALSE
  ))

  expect_identical(reads$CB, c("CELL-1", NA_character_))
  expect_identical(reads$ZZ, c(NA_character_, "CELL-1"))
})


test_that("default query regions are determined separately for each BAM", {
  skip_if_not_installed("Rsamtools")
  directories <- c(tempfile("hla-chr6-"), tempfile("hla-6-"))
  invisible(lapply(directories, dir.create))
  bam <- c(make_hla_test_bam(directories[[1L]]),
           make_hla_test_bam(directories[[2L]], chr6 = "6"))

  reads <- suppressMessages(get_hla_reads(
    bam, mhc_start = 1, mhc_end = 300, include_alt_contigs = FALSE,
    scores = FALSE
  ))

  expect_equal(nrow(reads), 4L)
  regions <- attr(reads, "hla_regions")
  expect_identical(as.character(GenomicRanges::seqnames(regions[[1L]])), "chr6")
  expect_identical(as.character(GenomicRanges::seqnames(regions[[2L]])), "6")
})


test_that("empty BAM queries can be pooled with nonempty or empty queries", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-empty-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  empty_sam <- file.path(directory, "empty.sam")
  writeLines(readLines(file.path(directory, "hla-test.sam"), n = 4L), empty_sam)
  empty_bam <- suppressMessages(Rsamtools::asBam(
    empty_sam, file.path(directory, "empty"), indexDestination = TRUE
  ))
  region <- GenomicRanges::GRanges("chr6", IRanges::IRanges(1, 130))

  reads <- suppressMessages(get_hla_reads(
    c(empty_bam, bam), regions = region, scores = FALSE
  ))
  expect_equal(nrow(reads), 1L)
  expect_identical(reads$readName, "2:hla_primary")
  expect_identical(reads$bam_file, bam)
  expect_length(attr(reads, "hla_regions"), 2L)

  empty <- suppressMessages(get_hla_reads(
    c(empty_bam, empty_bam), regions = region, scores = FALSE
  ))
  expect_equal(nrow(empty), 0L)
  expect_identical(empty$readName, character())
  expect_identical(empty$bam_file, character())
  expect_length(attr(empty, "hla_regions"), 2L)
})


test_that("BAM path vectors are validated before extraction", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-validation-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)

  for (invalid in list(character(), NA_character_, "", 1L,
                       c(bam, NA_character_), c(bam, ""),
                       c(bam, file.path(directory, "missing.bam")))) {
    expect_error(get_hla_reads(invalid), "paths to existing BAM files")
  }
  expect_error(get_hla_reads(), "paths to existing BAM files")
})


test_that("the typing workflow passes pooled reads to each requested gene", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("hla-pooled-typing-")
  dir.create(directory)
  bam <- make_hla_test_bam(directory)
  region <- GenomicRanges::GRanges("chr6", IRanges::IRanges(1, 130))
  workflow <- hla_typing_from_bam
  workflow_env <- new.env(parent = environment(workflow))
  workflow_env$hla_typing <- function(reads, hla_ref, ...) {
    list(n_reads = nrow(reads), gene = unique(hla_ref$gene))
  }
  environment(workflow) <- workflow_env

  result <- suppressMessages(workflow(
    bam = c(bam, bam),
    hla_ref = data.frame(gene = c("A", "B")),
    genes = c("A", "B"), regions = region, scores = FALSE
  ))

  expect_equal(nrow(result$reads), 2L)
  expect_named(result$typing, c("A", "B"))
  expect_identical(result$typing$A, list(n_reads = 2L, gene = "A"))
  expect_identical(result$typing$B, list(n_reads = 2L, gene = "B"))
  expect_identical(result$regions, attr(result$reads, "hla_regions"))
})
