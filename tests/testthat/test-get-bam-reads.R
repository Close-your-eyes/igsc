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

test_that("BAM subset output validates paths and protects existing files", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  directory <- tempfile("bam-subset-paths-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  bam <- make_bam_reads_fixture(directory)
  region <- GenomicRanges::GRanges("chr1:90-110")
  subset <- function(path) suppressMessages(get_bam_reads(
    bam, granges = region, return_bam_subset = path
  ))
  for (path in list(TRUE, 1, NA_character_, character(),
                   c("a.bam", "b.bam"), "", "  ")) {
    expect_error(subset(path), "must be one non-empty, non-NA path", fixed = TRUE)
  }
  expect_error(subset(file.path(directory, "subset.sam")), "must end in .bam",
               fixed = TRUE)
  expect_error(subset(file.path(directory, "missing", "subset.bam")),
               "output directory must exist and be writable", fixed = TRUE)
  expect_error(subset(bam), "must differ from the input BAM", fixed = TRUE)
  expect_error(subset(file.path(directory, ".", basename(bam))),
               "must differ from the input BAM", fixed = TRUE)

  destination <- file.path(directory, "subset.bam")
  for (existing in c(destination, paste0(destination, c(".bai", ".csi")),
                     file.path(directory, c("subset.bai", "subset.csi")))) {
    writeLines("keep me", existing)
    expect_error(subset(destination), "refusing to overwrite", fixed = TRUE)
    expect_identical(readLines(existing), "keep me")
    unlink(existing)
  }
  dir.create(destination)
  expect_error(subset(destination), "refusing to overwrite", fixed = TRUE)
  expect_equal(Rsamtools::countBam(bam)$records, 2)
})

test_that("BAM subset output validates ranges and requires an input index", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  directory <- tempfile("bam-subset-input-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  bam <- make_bam_reads_fixture(directory)
  destination <- file.path(directory, "subset.bam")
  subset <- function(region) suppressMessages(get_bam_reads(
    bam, granges = region, return_bam_subset = destination
  ))
  expect_error(subset(GenomicRanges::GRanges("1:90-110")),
               "absent from the BAM header: 1", fixed = TRUE)
  expect_error(subset(GenomicRanges::GRanges(c("chr1:90-110", "chr2:1-10"))),
               "absent from the BAM header: chr2", fixed = TRUE)
  expect_error(subset(GenomicRanges::GRanges("chr1:0-10")),
               "positive-width ranges starting at base 1 or later", fixed = TRUE)
  expect_error(subset(GenomicRanges::GRanges("chr1", IRanges::IRanges(100, width = 0))),
               "positive-width ranges starting at base 1 or later", fixed = TRUE)
  unlink(paste0(bam, ".bai"))
  expect_error(subset(GenomicRanges::GRanges("chr1:90-110")),
               "input BAM needs a readable index", fixed = TRUE)
  expect_false(file.exists(destination))
  expect_false(file.exists(paste0(destination, ".bai")))
})

test_that("BAM subset output retains complete overlapping alignments and tags", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  directory <- tempfile("bam-subset-success-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  bam <- make_bam_reads_fixture(directory)
  destination <- file.path(directory, "subset.bam")
  # A single-base overlap must keep the whole read and its tags unchanged.
  result <- suppressMessages(get_bam_reads(
    bam, granges = GenomicRanges::GRanges("chr1:201-201"),
    tags = character(), revcomp_minus_strand = TRUE,
    revcomp_plus_strand = TRUE, return_bam_subset = destination
  ))
  expect_identical(result, destination)
  expect_true(file.exists(destination))
  expect_true(file.exists(paste0(destination, ".bai")))
  expect_equal(Rsamtools::countBam(destination)$records, 1)
  param <- Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(), tag = "XY")
  actual <- Rsamtools::scanBam(destination, param = param)[[1]]
  expected <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(
    which = GenomicRanges::GRanges("chr1:201-201"),
    what = Rsamtools::scanBamWhat(), tag = "XY"
  ))[[1]]
  expect_identical(actual, expected)
  expect_identical(actual$qname, "read2")

  empty_destination <- file.path(directory, "empty.bam")
  suppressMessages(get_bam_reads(
    bam, granges = GenomicRanges::GRanges("chr1:900-1000"),
    return_bam_subset = empty_destination
  ))
  expect_true(file.exists(paste0(empty_destination, ".bai")))
  expect_equal(Rsamtools::countBam(empty_destination)$records, 0)
})
