gtf_test_record <- function(seqname, gene_id, gene_name, feature = "gene",
                            extra = "") {
  paste(
    seqname, "test", feature, 1, 10, ".", "+", ".",
    paste0(
      'gene_id "', gene_id, '"; gene_name "', gene_name, '";', extra
    ),
    sep = "\t"
  )
}


test_that("sequence filtering is exact with interleaved records and spaced paths", {
  directory <- tempfile("gtf test ")
  dir.create(directory)
  path <- file.path(directory, "genes with spaces.gtf")
  records <- c(
    gtf_test_record("chr1", "g1", "A"),
    gtf_test_record("chr2", "g2", "B"),
    gtf_test_record("chr1", "g3", "C"),
    gtf_test_record("chr1[alt]", "g4", "D")
  )
  writeLines(records, path)

  result <- read_gtf(path, seqnames = "chr1", process_attr_col = FALSE)$gtf
  expect_identical(as.character(result$seqname), c("chr1", "chr1"))
  expect_equal(nrow(result), 2L)

  exact <- read_gtf(path, seqnames = "chr1[alt]",
                    process_attr_col = FALSE)$gtf
  expect_identical(as.character(exact$seqname), "chr1[alt]")

  compressed_path <- paste0(path, ".gz")
  connection <- gzfile(compressed_path, open = "wt")
  writeLines(records, connection)
  close(connection)
  compressed <- read_gtf(compressed_path, seqnames = "chr1",
                         process_attr_col = FALSE)$gtf
  expect_identical(as.character(compressed$seqname), c("chr1", "chr1"))
})


test_that("gene-name filtering uses the gene_name attribute literally", {
  path <- tempfile(fileext = ".gtf")
  writeLines(
    c(
      gtf_test_record("chr6", "g1", "HLA-A"),
      gtf_test_record("chr6", "g2", "HLAXA"),
      gtf_test_record("chr6", "HLA-A", "OTHER")
    ),
    path
  )

  expect_error(
    read_gtf(path, gene_names = "HLA.A", process_attr_col = FALSE),
    "No rows left in gtf"
  )
})


test_that("literal and regex gene-name filters do not match other attributes", {
  path <- tempfile(fileext = ".gtf")
  writeLines(
    c(
      gtf_test_record("chr6", "g1", "HLA-A"),
      gtf_test_record("chr6", "g2", "HLAXA"),
      gtf_test_record("chr6", "HLA-A", "OTHER")
    ),
    path
  )

  literal <- read_gtf(path, gene_names = "HLA-A",
                      process_attr_col = FALSE)$gtf
  expect_equal(nrow(literal), 1L)
  expect_match(literal$attribute, 'gene_id "g1"', fixed = TRUE)

  regex <- read_gtf(path, gene_names = "HLA.A",
                    gene_names_full_match = FALSE,
                    process_attr_col = FALSE)$gtf
  expect_equal(nrow(regex), 2L)
  expect_false(any(grepl('gene_id "HLA-A"', regex$attribute, fixed = TRUE)))
})


test_that("annotated gene names remain stable across contigs", {
  path <- tempfile(fileext = ".gtf")
  writeLines(
    c(
      gtf_test_record("chr6", "main", "HLA-A"),
      gtf_test_record("chr6_alt", "alt", "HLA-A")
    ),
    path
  )

  result <- read_gtf(path, features = "gene", gene_names = "HLA-A")$gtf
  expect_identical(as.character(result$gene_name), c("HLA-A", "HLA-A"))

  unique_names <- read_gtf(
    path,
    features = "gene",
    gene_names = "HLA-A",
    process_attr_col_args_repl = list(unique_gene_names = TRUE)
  )$gtf
  expect_setequal(unique_names$gene_name, c("HLA-A--1", "HLA-A--2"))
})


test_that("fill_na handles gene-only annotations", {
  path <- tempfile(fileext = ".gtf")
  writeLines(gtf_test_record("chr1", "g1", "A"), path)

  result <- read_gtf(
    path,
    process_attr_col_args_repl = list(fill_na = TRUE)
  )$gtf
  expect_identical(result$transcript_id, "g1")
  expect_identical(result$transcript_name, "A")
})


test_that("fill_na preserves exon-number and exon-ID parsing", {
  path <- tempfile(fileext = ".gtf")
  writeLines(
    c(
      gtf_test_record(
        "chr1", "g1", "A", feature = "exon",
        extra = ' exon_number "1; exon_id e1";'
      ),
      gtf_test_record("chr1", "g1", "A")
    ),
    path
  )

  result <- read_gtf(
    path,
    process_attr_col_args_repl = list(fill_na = TRUE)
  )$gtf
  expect_equal(result$exon_number, c(1, 0))
  expect_identical(result$exon_id, c("e1", "0"))
})


test_that("NULL overrides remain explicit in attribute-processing arguments", {
  path <- tempfile(fileext = ".gtf")
  writeLines(gtf_test_record("chr1", "g1", "A"), path)

  result <- read_gtf(
    path,
    process_attr_col_args = list(
      attr_keep = "gene_id",
      gene_name_force = "gene_id"
    ),
    process_attr_col_args_repl = list(gene_name_force = NULL)
  )$gtf
  expect_false("gene_name" %in% names(result))
})
