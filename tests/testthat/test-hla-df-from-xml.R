make_strand_hla_xml <- function(genes) {
  path <- tempfile(fileext = ".xml")
  records <- vapply(genes, function(gene) paste0(
    '<allele name="HLA-', gene, '*01:01">',
    '<hla_g_group name="None"/><hla_p_group name="None"/>',
    '<sequence><nucsequence>AACCGGTTAGCA</nucsequence>',
    '<feature name="Exon 2"><coordinates start="2" end="4"/></feature>',
    '<feature name="Exon 3"><coordinates start="8" end="10"/></feature>',
    '<feature name="Intron 2"><coordinates start="5" end="7"/></feature>',
    '<feature name="5UTR"><coordinates start="1" end="1"/></feature>',
    '<feature name="Translation"><coordinates start="1" end="12" frame="1"/></feature>',
    '</sequence></allele>'
  ), character(1))
  writeLines(c('<hla>', records, '</hla>'), path)
  path
}

test_that("XML orientation is retained by default and reported", {
  path <- make_strand_hla_xml(c("A", "B"))
  on.exit(unlink(path))
  expect_message(result <- hla_df_from_xml(path),
                 "rev_comp_minus = FALSE.*not reverse complemented")
  expect_identical(result$seq, rep("AACCGGTTAGCA", 2L))
  expect_identical(result$seq_Exon2_3, rep("ACCTAG", 2L))
  expect_identical(result$Exon2, rep("2_4", 2L))
})

test_that("only the supplied minus-strand genes are reverse complemented", {
  minus <- c("B", "C", "DMA", "DMB", "DOA", "DOB", "DPA1", "DQB1",
             "DQB2", "DRB1", "DRB5")
  plus <- c("A", "DPB1", "DQA1", "DQA2", "DRA", "E", "F", "G", "DRB10")
  path <- make_strand_hla_xml(c(minus, plus))
  on.exit(unlink(path))
  expect_message(result <- hla_df_from_xml(path, rev_comp_minus = TRUE),
                 "rev_comp_minus = TRUE.*11 minus-strand HLA alleles")
  negative <- seq_along(minus)
  positive <- length(minus) + seq_along(plus)
  expect_identical(result$seq[negative], rep("TGCTAACCGGTT", length(minus)))
  expect_identical(result$seq_Exon2_3[negative], rep("CTAGGT", length(minus)))
  expect_identical(result$seq[positive], rep("AACCGGTTAGCA", length(plus)))
  expect_identical(result$seq_Exon2_3[positive], rep("ACCTAG", length(plus)))
  expect_identical(result$Exon2[negative], rep("9_11", length(minus)))
  expect_identical(result$Exon3[negative], rep("3_5", length(minus)))
  expect_identical(result$Intron2[negative], rep("6_8", length(minus)))
  expect_identical(result$`5UTR`[negative], rep("12_12", length(minus)))
  expect_identical(result$Exon2[positive], rep("2_4", length(plus)))
  expect_identical(extract_cds(result[negative, ], exons = c(3, 2)),
                   result$seq_Exon2_3[negative])
  expect_identical(result$seq_length, rep(12L, nrow(result)))
  expect_identical(result$seq_length_Exon2_3, rep(6L, nrow(result)))
})

test_that("reverse complement works with filtering and boundary removal", {
  path <- make_strand_hla_xml(c("A", "B", "B"))
  on.exit(unlink(path))
  result <- suppressMessages(hla_df_from_xml(
    path, rev_comp_minus = TRUE, unique_exon2_3 = TRUE,
    rm_intron_exon_bounds = TRUE
  ))
  expect_equal(nrow(result), 2L)
  expect_identical(result$seq_Exon2_3, c("ACCTAG", "CTAGGT"))
  expect_false(any(grepl("^Intron|^Exon|UTR", names(result))))
})

test_that("the reverse complement option is validated and reports no matches", {
  path <- make_strand_hla_xml("A")
  on.exit(unlink(path))
  for (invalid in list(NA, NULL, c(TRUE, FALSE), "TRUE", 1L)) {
    expect_error(hla_df_from_xml(path, rev_comp_minus = invalid),
                 "`rev_comp_minus` must be TRUE or FALSE.", fixed = TRUE)
  }
  expect_message(result <- hla_df_from_xml(path, rev_comp_minus = TRUE),
                 "0 minus-strand HLA alleles")
  expect_identical(result$seq, "AACCGGTTAGCA")
})
