make_test_tcr_protein <- function() {
  structure(
    list(
      sequence = "ACDE",
      length_aa = 4L,
      chain = "TRB",
      cdr3 = "CDE",
      segments = c(V = "TRBV1*01", D = NA_character_,
                   J = "TRBJ1-1*01", C = "TRBC1*01")
    ),
    class = "tcr_protein"
  )
}


test_that("tcr_protein print method dispatches", {
  protein <- make_test_tcr_protein()

  expect_s3_class(protein, "tcr_protein")
  expect_output(print(protein), "Full TRB TCR protein: 4 aa", fixed = TRUE)
  expect_output(print(protein), "ACDE", fixed = TRUE)
  expect_identical(
    getS3method("print", "tcr_protein", optional = TRUE),
    igsc:::print.tcr_protein
  )
})


test_that("as_fasta is exported and dispatches for tcr_protein", {
  protein <- make_test_tcr_protein()

  expect_identical(
    as_fasta(protein, header = "example", width = 2L),
    ">example\nAC\nDE\n"
  )
  expect_match(as_fasta(protein), "^>TCR_TRB\\|")
  expect_identical(
    getS3method("as_fasta", "tcr_protein", optional = TRUE),
    igsc:::as_fasta.tcr_protein
  )
})
