check_fasta_bounds <- function(directory) {
  records <- c(">chr1", rep("ACGT", 100), ">chr2", "TGCA", "AAAA")
  for (filename in c("genome.fa", "genome ' $literal; [1] `literal`.fa")) {
    path <- file.path(directory, filename)
    writeLines(records, path)
    compressed_path <- paste0(path, ".gz")
    connection <- gzfile(compressed_path, "wt")
    writeLines(records, connection)
    close(connection)

    for (input in c(path, compressed_path)) {
      bounds <- get_fasta_seq_bounds(input)
      expect_equal(bounds$start_line, c(1, 102))
      expect_equal(bounds$end_line, c(101, 104))
      expect_equal(bounds$n_seq_lines, c(100, 2))
      expect_equal(bounds$bp_approx, c(400, 8))
      expect_identical(as.character(bounds$seqname), c("chr1", "chr2"))

      headers <- get_fasta_names(input, prep_chr = FALSE)
      expect_equal(headers$start_line, c(1, 102))
      expect_identical(headers$fastaname, c("chr1", "chr2"))
    }
  }
}

test_that("FASTA boundaries handle gzip input and literal shell characters", {
  skip_if_not(nzchar(Sys.which("rg")))
  skip_if_not(nzchar(Sys.which("gunzip")))
  directory <- tempfile("fasta bounds ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))

  check_fasta_bounds(directory)
})

test_that("FASTA boundaries fall back to grep without ripgrep", {
  skip_on_os("windows")
  commands <- Sys.which(c("grep", "gunzip", "gzip"))
  skip_if_not(all(nzchar(commands)))
  directory <- tempfile("fasta fallback ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  bin <- file.path(directory, "bin")
  dir.create(bin)
  linked <- file.symlink(unname(commands), file.path(bin, names(commands)))
  skip_if_not(all(linked))
  original_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = original_path), add = TRUE)
  Sys.setenv(PATH = bin)
  expect_identical(unname(Sys.which("rg")), "")

  check_fasta_bounds(directory)
})
