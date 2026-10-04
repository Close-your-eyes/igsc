.mafft_test_executable <- function() {
  skip_if_not_installed("Biostrings")
  executable <- tryCatch(.find_mafft_executable(), error = function(e) NULL)
  if (is.null(executable)) skip("MAFFT executable is not installed")
  executable
}
.mafft_test_content <- function(x, y) {
  expect_identical(class(y), class(x))
  expect_identical(names(y), names(x))
  expect_equal(length(unique(nchar(as.character(y)))), 1)
  expect_identical(unname(gsub("-", "", as.character(y), fixed = TRUE)),
                   unname(gsub("-", "", as.character(x), fixed = TRUE)))
}

test_that("DNA aligns with MacPorts MAFFT and keeps names and order", {
  executable <- .mafft_test_executable()
  x <- Biostrings::DNAStringSet(stats::setNames(
    c("AAAACCCCGGGG", "AAAATCCCCGGGG", "AAAACCCGGGG"),
    c("name with spaces", "duplicate", "duplicate")))
  expect_no_message(y <- align_seqs_mafft(x, processors = 1, verbose = FALSE, mafft = executable))
  .mafft_test_content(x, y)
  expect_true(any(grepl("-", as.character(y), fixed = TRUE)))
  .mafft_test_content(x, align_seqs_mafft(x, verbose = FALSE, mafft_args = "--reorder"))
})

test_that("AA input forces protein alignment even for nucleotide-like residues", {
  executable <- .mafft_test_executable()
  for (sequences in list(c("MKVLW", "MKVLAW", "MKLW"), c("ACGT", "ACGGT", "AGT"),
                        c("MKVLW*", "MKVLAW*", "MKLW*"))) {
    x <- Biostrings::AAStringSet(sequences)
    y <- align_seqs_mafft(x, verbose = FALSE, mafft = executable)
    .mafft_test_content(x, y)
  }
})

test_that("RNA sequences retain RNA type and uracil", {
  executable <- .mafft_test_executable()
  x <- Biostrings::RNAStringSet(c(first = "AAAACCCCUUUU", second = "AAAACCCUUUU"))
  y <- align_seqs_mafft(x, verbose = FALSE, mafft = executable)
  .mafft_test_content(x, y)
})

test_that("strategies and native scoring arguments work", {
  executable <- .mafft_test_executable()
  x <- Biostrings::DNAStringSet(c(a = "AAACCCGGG", b = "AAAACCCGGG", c = "AAACCGGG"))
  for (strategy in c("linsi", "ginsi", "einsi", "fftns2")) {
    .mafft_test_content(x, align_seqs_mafft(x, strategy = strategy, verbose = FALSE,
                                          mafft = executable, mafft_args = c("--op", "2")))
  }
  .mafft_test_content(x, align_seqs_mafft(x, processors = 2, verbose = FALSE, mafft = executable))
})

test_that("existing gaps and absent names are handled consistently", {
  executable <- .mafft_test_executable()
  x <- Biostrings::DNAStringSet(c("AA--CC", "AACCC"))
  .mafft_test_content(x, align_seqs_mafft(x, verbose = FALSE, mafft = executable))
  single <- Biostrings::AAStringSet(c("only one" = "M-KL"))
  expect_identical(align_seqs_mafft(single, verbose = FALSE, mafft = "/missing/mafft"),
                   Biostrings::AAStringSet(c("only one" = "MKL")))
})

test_that("invalid inputs and unsupported DECIPHER arguments fail clearly", {
  skip_if_not_installed("Biostrings")
  x <- Biostrings::DNAStringSet(c("AA", "AC"))
  expect_error(align_seqs_mafft(c("AA", "AC")), "myXStringSet must be")
  expect_error(align_seqs_mafft(Biostrings::DNAStringSet()), "must not be empty")
  expect_error(align_seqs_mafft(Biostrings::DNAStringSet(c("--", "AA"))), "non-gap")
  expect_error(align_seqs_mafft(x, iterations = 2), "Unsupported argument.*iterations")
  expect_error(align_seqs_mafft(x, gapOpening = -18), "Unsupported argument.*gapOpening")
  for (bad in list(0, -1, 1.5, NA, Inf, c(1, 2), "2")) {
    expect_error(align_seqs_mafft(x, processors = bad), "processors")
  }
  expect_error(align_seqs_mafft(x, verbose = NA), "verbose")
  expect_error(align_seqs_mafft(x, mafft = "/missing/mafft"), "executable not found")
  expect_error(align_seqs_mafft(x, mafft_args = NA_character_), "mafft_args")
  for (flag in c("--clustalout", "--thread", "--addfragments", "--adjustdirection", "--amino")) {
    expect_error(align_seqs_mafft(x, mafft_args = flag), "Unsupported mafft_args")
  }
})

test_that("external failures include diagnostics and clean temporary files", {
  skip_on_os("windows")
  skip_if_not_installed("Biostrings")
  folder <- tempfile("fake mafft ")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  executable <- file.path(folder, "mafft executable")
  writeLines(c("#!/bin/sh", "echo 'test MAFFT diagnostic' >&2", "exit 17"), executable)
  Sys.chmod(executable, "0755")
  before <- list.files(tempdir(), pattern = "^igsc-mafft-")
  x <- Biostrings::DNAStringSet(c("AA", "AC"))
  expect_error(align_seqs_mafft(x, mafft = executable, verbose = FALSE),
               "exit status 17.*\\n.*test MAFFT diagnostic")
  expect_identical(list.files(tempdir(), pattern = "^igsc-mafft-"), before)
  writeLines(c("#!/bin/sh", "for arg do input=$arg; done", 'cat "$input"'), executable)
  expect_identical(align_seqs_mafft(x, mafft = executable, verbose = FALSE), x)
  expect_identical(list.files(tempdir(), pattern = "^igsc-mafft-"), before)
  writeLines(c("#!/bin/sh", "previous=''; found=0", 
               'for arg do if [ "$previous" = "--thread" ] && [ "$arg" = "-1" ]; then found=1; fi; previous=$arg; input=$arg; done',
               '[ "$found" = 1 ] || exit 19', 'cat "$input"'), executable)
  expect_identical(align_seqs_mafft(x, processors = NULL, mafft = executable, verbose = FALSE), x)
  writeLines(c("#!/bin/sh", "printf '>wrong_id\\nAAAA\\n'"), executable)
  expect_error(align_seqs_mafft(x, mafft = executable, verbose = FALSE), "expected sequence identifiers")
  expect_identical(list.files(tempdir(), pattern = "^igsc-mafft-"), before)
})
