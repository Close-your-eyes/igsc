test_that("a missing strand is generated in antiparallel display order", {
  skip_if_not_installed("Biostrings")
  expected <- c("top     5' ATGC 3'", "           ||||",
                "bottom  3' TACG 5'")
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGC", col_out = FALSE
  )), expected)
  expect_identical(capture.output(dna_double_strand_print(
    bottom = "TACG", col_out = FALSE
  )), expected)

  capture.output(value <- withVisible(dna_double_strand_print(
    top = " aR-y.N ", col_out = FALSE
  )))
  expect_false(value$visible)
  expect_identical(value$value, c(top = " aR-y.N ", bottom = " TY-R.N "))
})

test_that("supplied strands preserve case, gaps, mismatches, and overhangs", {
  expect_identical(capture.output(dna_double_strand_print(
    top = "  ATGCaa", bottom = "GGTACG", col_out = FALSE
  )), c("top     5'   ATGCaa 3'", "             ||||  ",
        "bottom  3' GGTACG   5'"))

  expect_identical(capture.output(dna_double_strand_print(
    top = "ACGN-.", bottom = "TAAN-.", col_out = FALSE
  )), c("top     5' ACGN-. 3'", "           |     ",
        "bottom  3' TAAN-. 5'"))

  capture.output(value <- dna_double_strand_print(
    top = "AcG", bottom = " t", col_out = FALSE
  ))
  expect_identical(value, c(top = "AcG", bottom = " t"))
})

test_that("wrapping retains overhang columns and labels on both strands", {
  expect_identical(capture.output(dna_double_strand_print(
    top = "  ATGCaa", bottom = "GGTACG", linewidth = 4, col_out = FALSE
  )), c("top     5'   AT 3'", "             ||", "bottom  3' GGTA 5'",
        "", "top     5' GCaa 3'", "           ||  ", "bottom  3' CG   5'"))

  expect_identical(capture.output(dna_double_strand_print(
    top = "A", bottom = "TAG", linewidth = 2, col_out = FALSE
  )), c("top     5' A  3'", "           | ", "bottom  3' TA 5'",
        "", "top     5'   3'", "            ", "bottom  3' G 5'"))
})

test_that("right alignment pads either strand before wrapping", {
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGCAA", bottom = "CGTT", align = "right", linewidth = 4,
    col_out = FALSE
  )), c("top     5' ATGC 3'", "             ||", "bottom  3'   CG 5'",
        "", "top     5' AA 3'", "           ||", "bottom  3' TT 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "GC", bottom = "TACG", align = "right", col_out = FALSE
  )), c("top     5'   GC 3'", "             ||", "bottom  3' TACG 5'"))
  capture.output(value <- dna_double_strand_print(
    top = " GC", bottom = "TACG", align = "right", col_out = FALSE
  ))
  expect_identical(value, c(top = " GC", bottom = "TACG"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGC", bottom = "TA", align = "left", col_out = FALSE
  )), capture.output(dna_double_strand_print(
    top = "ATGC", bottom = "TA", col_out = FALSE
  )))
  expect_error(dna_double_strand_print(top = "AT", align = "center"), "arg")
})

test_that("output limits select aligned columns before wrapping", {
  skip_if_not_installed("Biostrings")
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGCAA", first_n = 3, col_out = FALSE
  )), c("top     5' ATG 3'", "           |||", "bottom  3' TAC 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    bottom = "TACGTT", last_n = 3, col_out = FALSE
  )), c("top     5' CAA 3'", "           |||", "bottom  3' GTT 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGCAA", bottom = "CGTT", align = "right", first_n = 3,
    col_out = FALSE
  )), c("top     5' ATG 3'", "             |", "bottom  3'   C 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGCAA", bottom = "CGTT", align = "right", last_n = 3,
    linewidth = 2, col_out = FALSE
  )), c("top     5' CA 3'", "           ||", "bottom  3' GT 5'",
        "", "top     5' A 3'", "           |", "bottom  3' T 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGCAA", bottom = "TACG", last_n = 3, col_out = FALSE
  )), c("top     5' CAA 3'", "           |  ", "bottom  3' G   5'"))
  expected <- capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", col_out = FALSE
  ))
  expect_identical(capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", first_n = 10, col_out = FALSE
  )), expected)
  expect_identical(capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", last_n = 10, col_out = FALSE
  )), expected)
  capture.output(value <- dna_double_strand_print(
    top = "ATGCAA", bottom = "CGTT", last_n = 1, col_out = FALSE
  ))
  expect_identical(value, c(top = "ATGCAA", bottom = "CGTT"))
  expect_error(dna_double_strand_print(top = "AT", first_n = 1, last_n = 1),
               "at most one")
  for (n in list(0, -1, 1.5, Inf, NA_real_, numeric(), c(2, 3), "2", TRUE)) {
    expect_error(dna_double_strand_print(top = "AT", first_n = n), "first_n")
    expect_error(dna_double_strand_print(top = "AT", last_n = n), "last_n")
  }
})

test_that("end labels and pipes can be disabled independently", {
  expect_identical(capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", print_ends = FALSE, col_out = FALSE
  )), c("top     AT", "        ||", "bottom  TA"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", print_pipes = FALSE, col_out = FALSE
  )), c("top     5' AT 3'", "bottom  3' TA 5'"))
  expect_identical(capture.output(dna_double_strand_print(
    top = "AT", bottom = "TA", print_ends = FALSE, print_pipes = FALSE,
    col_out = FALSE
  )), c("top     AT", "bottom  TA"))
})

test_that("DNAString inputs and custom row labels are supported", {
  skip_if_not_installed("Biostrings")
  expect_identical(capture.output(dna_double_strand_print(
    top = Biostrings::DNAStringSet(c(forward = "AT")),
    bottom = c(reverse = "TA"), col_out = FALSE
  )), c("forward  5' AT 3'", "            ||", "reverse  3' TA 5'"))
  capture.output(value <- dna_double_strand_print(
    bottom = Biostrings::DNAString("TA"), col_out = FALSE
  ))
  expect_identical(value, c(top = "AT", bottom = "TA"))
})

test_that("file output agrees with console output and preserves connections", {
  path <- tempfile(fileext = ".txt")
  other <- file(tempfile(), open = "wt")
  on.exit(close(other), add = TRUE)
  on.exit(unlink(path), add = TRUE)
  expected <- capture.output(dna_double_strand_print(
    top = "ATGC", bottom = "TA", linewidth = 3, col_out = FALSE
  ))
  sink_count <- sink.number()
  expect_identical(capture.output(dna_double_strand_print(
    top = "ATGC", bottom = "TA", linewidth = 3, col_out = FALSE,
    out_file = path
  )), character())
  expect_identical(readLines(path), expected)
  expect_true(isOpen(other))
  expect_identical(sink.number(), sink_count)
})

test_that("coloring leaves the displayed bases and spacing intact", {
  skip_if_not_installed("crayon")
  skip_if_not_installed("RColorBrewer")
  previous <- options(cli.num_colors = 256, crayon.enabled = TRUE, crayon.colors = 256)
  on.exit(options(previous), add = TRUE)
  plain <- capture.output(dna_double_strand_print(
    top = "  AT-Gn", bottom = "GGTANC", col_out = FALSE
  ))
  colored <- capture.output(dna_double_strand_print(
    top = "  AT-Gn", bottom = "GGTANC", col_out = TRUE
  ))
  expect_identical(crayon::strip_style(colored), plain)
  expect_true(any(crayon::has_style(colored)))
})

test_that("invalid inputs fail with useful messages", {
  expect_error(dna_double_strand_print(), "at least one")
  for (strand in list("", NA_character_, character(), c("AT", "GC"),
                      "AX", "AU", "A\nT", "---", "   ")) {
    expect_error(dna_double_strand_print(top = strand), "nonempty DNA")
  }
  expect_error(dna_double_strand_print(top = 12), "DNA string")
  for (width in list(0, -1, 1.5, Inf, NA_real_, numeric(), c(2, 3), "2")) {
    expect_error(dna_double_strand_print(top = "AT", linewidth = width),
                 "positive integer")
  }
  expect_error(dna_double_strand_print(top = "AT", print_ends = NA), "print_ends")
  expect_error(dna_double_strand_print(top = "AT", print_pipes = 1), "print_pipes")
  expect_error(dna_double_strand_print(top = "AT", col_out = logical()), "col_out")
  expect_error(dna_double_strand_print(top = "AT", out_file = NA), "out_file")
})
