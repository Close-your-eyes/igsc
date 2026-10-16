make_bam_reads_fixture <- function(directory) {
  sam <- file.path(directory, "reads.sam")
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:2000",
    "read1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tAACC\tIJKL\tCB:Z:CELL1\tXY:B:i,1,2",
    "read2\t16\tchr1\t200\t60\t4M\t*\t0\t0\tGGTT\t*\tXY:B:i,3,4,5"
  ), sam)
  Rsamtools::asBam(sam, file.path(directory, "reads"), indexDestination = TRUE)
}

test_that("BAM extraction handles absent, scalar, and array tags", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  directory <- tempfile("bam-reads-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  bam <- make_bam_reads_fixture(directory)
  region <- GenomicRanges::GRanges("chr1:1-1000")
  extract <- function(tags, ...) suppressMessages(
    get_bam_reads(bam, granges = region, tags = tags, scores = FALSE, ...)
  )

  absent <- extract("UB")
  expect_identical(absent$qname, c("read1", "read2"))
  expect_false("UB" %in% names(absent))
  expect_equal(extract(character()), absent)

  tagged <- extract(c("CB", "XY", "UB"))
  expect_identical(tagged$qname, absent$qname)
  expect_identical(tagged$CB, c("CELL1", NA_character_))
  expect_type(tagged$XY, "list")
  expect_equal(tagged$XY, list(c(1L, 2L), c(3L, 4L, 5L)))
  expect_false("UB" %in% names(tagged))

  reversed <- extract("XY", revcomp_plus_strand = TRUE)
  expect_identical(reversed$seq, c("GGTT", "GGTT"))
  expect_identical(reversed$qual[1], "LKJI")

  regions <- GenomicRanges::GRanges(c("chr1:90-110", "chr1:190-210", "chr1:900-910"))
  multi <- suppressMessages(get_bam_reads(
    bam, granges = regions, tags = c("CB", "XY"), scores = FALSE
  ))
  expect_identical(multi$qname, absent$qname)
  expect_identical(multi$CB, tagged$CB)
  expect_equal(multi$XY, tagged$XY)
})

test_that("missing BAM qualities yield unavailable quality summaries", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  directory <- tempfile("bam-quality-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  bam <- make_bam_reads_fixture(directory)
  reads <- suppressMessages(get_bam_reads(
    bam, granges = GenomicRanges::GRanges("chr1:1-1000"), tags = character()
  ))
  expect_identical(reads$readQualNum[[1]], 40:43)
  expect_identical(reads$readQualNum[[2]], integer())
  expect_equal(reads$minQual, c(40L, NA_integer_))
  expect_equal(reads$meanQual, c(41.5, NA_real_))
  expect_identical(reads$n_belowQ30, c(0L, NA_integer_))

  empty <- suppressMessages(get_bam_reads(
    bam, granges = GenomicRanges::GRanges("chr1:900-1000"), tags = "XY"
  ))
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("qname", "minQual", "meanQual", "n_belowQ30") %in% names(empty)))
})

test_that("empty BAM ranges are rejected before reading alignments", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  # An existing non-BAM file verifies that validation precedes scanBam.
  bam <- tempfile(fileext = ".bam")
  file.create(bam)
  on.exit(unlink(bam), add = TRUE)
  expect_error(
    get_bam_reads(bam, granges = GenomicRanges::GRanges()),
    "must be a non-empty GRanges object", fixed = TRUE
  )
})
