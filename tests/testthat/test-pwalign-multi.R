.pwalign_test_dependencies <- function() {
  for (package in c("Biostrings", "brathering", "pwalign", "ggplot2", "colrr",
                    "Peptides", "RColorBrewer", "crayon", "ggrepel", "scales")) {
    testthat::skip_if_not_installed(package)
  }
}

.pwalign_test_run <- function(patterns, subject = "AAAACCCCGGGGTTTT", ...) {
  suppressMessages(pwalign_multi(
    subject = stats::setNames(subject, "ref"), patterns = patterns,
    seq_type = "NT", verbose = FALSE,
    pair_aln_args = list(gapOpening = 2, gapExtension = 1), ...
  ))
}

.pwalign_test_reference <- function(result, subject) {
  df <- result$data_wide
  expect_identical(paste0(df$ref[!is.na(df$subject.position)], collapse = ""), subject)
  expect_equal(df$subject.position[!is.na(df$subject.position)], seq_len(nchar(subject)))
  expect_true(all(df$ref[is.na(df$subject.position)] == "-"))
  expect_identical(df$position, seq_len(nrow(df)))
}

test_that("ordinary alignments and pattern limits work in a clean session", {
  .pwalign_test_dependencies()
  expect_no_warning(x <- .pwalign_test_run(c(p1 = "AAAA", p2 = "CCCC")))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_equal(which(!is.na(x$data_wide$p2)), 5:8)
  expect_no_warning(ggplot2::ggplot_build(x$plot))
})

test_that("disjoint indels retain the full subject and downstream positions", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(p1 = "AAAATCCCC", p2 = "TTTT"))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_identical(paste0(x$data_wide$ref, collapse = ""), "AAAA-CCCCGGGGTTTT")
  expect_identical(x$data_wide$p1[1:9], strsplit("AAAATCCCC", "")[[1]])
  expect_equal(which(!is.na(x$data_wide$p2)), 14:17)
})

test_that("unsupported overlaps fail regardless of pattern order", {
  .pwalign_test_dependencies()
  patterns <- c(gapped = "AAAATCCCC", spanning = "AACCCC")
  for (p in list(patterns, rev(patterns))) {
    expect_error(.pwalign_test_run(p), "Overlapping indel.*gapped")
  }
})

test_that("indel repair uses pattern coordinates and recomputes alignments", {
  .pwalign_test_dependencies()
  subject <- "GGGGAAAACCCCTTTT"
  x <- .pwalign_test_run(c(gapped = "AAAATCCCC", spanning = "AACCCC"),
                         subject, fix_subject_indels = TRUE)
  .pwalign_test_reference(x, subject)
  expect_identical(unname(as.character(x$pairwise_alignments@pattern@unaligned)),
                   c("AAAA", "AACCCC"))
  expect_false(anyNA(x$data_wide$subject.position))
  expect_equal(which(!is.na(x$data_wide$gapped)), 5:8)
})

test_that("identical shared gaps occupy one set of reference columns", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(p1 = "AAAATCCCC", p2 = "AAAATCCCC"))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_equal(sum(is.na(x$data_wide$subject.position)), 1)
  expect_identical(x$data_wide$p1, x$data_wide$p2)
})

test_that("a gap outside another overlapping alignment is supported", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(p1 = "AAAATCCCC", p2 = "CCCC"))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_equal(which(!is.na(x$data_wide$p2)), 6:9)
})

test_that("gap positions account for earlier gaps in the same alignment", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(p1 = "AAAATCCCCTGGGG", p2 = "TTTT"))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_identical(paste0(x$data_wide$ref, collapse = ""), "AAAA-CCCC-GGGGTTTT")
  expect_equal(which(!is.na(x$data_wide$p2)), 15:18)
  y <- .pwalign_test_run(c(p1 = "AAAATCCCCTGGGG", p2 = "AAAATCCCCTGGGG"))
  expect_equal(sum(is.na(y$data_wide$subject.position)), 2)
  expect_identical(y$data_wide$p1, y$data_wide$p2)
})

test_that("removing inducing patterns clears their gap information", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(p1 = "AAAATCCCC", p2 = "TTTT"),
                         rm_indel_inducing_pattern = TRUE)
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_false(anyNA(x$data_wide$subject.position))
  expect_equal(which(!is.na(x$data_wide$p2)), 13:16)
  expect_identical(names(x$pattern_indel_inducing), "p1")
  expect_error(.pwalign_test_run(c(p1 = "AAAATCCCC"),
                                rm_indel_inducing_pattern = TRUE),
               "No patterns left")
})

test_that("ambiguous subjects disable mismatch filtering without removing patterns", {
  .pwalign_test_dependencies()
  subject <- Biostrings::DNAStringSet(c(ref = "AAAANCCCC"))
  patterns <- Biostrings::DNAStringSet(c(p1 = "AAAA", p2 = "AANC"))
  expect_message(x <- check_for_invalid_chars(subject, patterns, 1), "set to NA")
  expect_identical(x$patterns, patterns)
  expect_true(is.na(x$max_mismatch))
  expect_null(x$patterns_invalid)
  expect_null(x$pattern_mismatching_return)
  expect_no_warning(.pwalign_test_run(c(p1 = "AAAA", p2 = "CCCC"),
                                      "AAAANCCCC", max_mismatch = 1))
})

test_that("single alignments preserve sequence content across alignment modes", {
  .pwalign_test_dependencies()
  subject <- "GGGGAAAACCCCTTTT"
  for (mode in c("global-local", "global", "local", "overlap", "local-global")) {
    x <- .pwalign_test_run(c(p1 = "AAAATCCCC"), subject, type = mode)
    .pwalign_test_reference(x, subject)
    expect_identical(paste0(stats::na.omit(x$data_wide$p1), collapse = ""),
                     unname(as.character(x$pairwise_alignments@pattern)))
    expect_no_warning(ggplot2::ggplot_build(x$plot))
  }
})

test_that("sanitized names remain distinct and original labels survive", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c("a-b" = "AAAA", "a.b" = "CCCC", "a.b.1" = "GGGG"))
  expect_true(all(c("a-b", "a.b", "a.b.1") %in% as.character(x$data$seq.name)))
  expect_identical(unname(x$pairwise_alignments@pattern@unaligned@ranges@NAMES),
                   c("a-b", "a.b", "a.b.1"))
})

test_that("singleton groups pass and within-group overlaps fail", {
  .pwalign_test_dependencies()
  expect_no_error(.pwalign_test_run(list(g1 = c("AAAA", "CCCC"), g2 = "GGGG")))
  expect_no_error(.pwalign_test_run(list(g1 = c("AAAA", "TTTA"), g2 = c("CCCC", "GGGG")),
                                    max_mismatch = 0))
  expect_error(.pwalign_test_run(list(g1 = c("AAAAC", "CCCC"))),
               "groups have patterns with overlapping")
})

test_that("terminal overhangs follow the pairwise alignment ranges", {
  .pwalign_test_dependencies()
  for (pattern in c("TTAAAACCCC", "AAAACCCCTT", "TTAAAACCCCTT")) {
    x <- .pwalign_test_run(c(p1 = pattern), "AAAACCCC", type = "global")
    .pwalign_test_reference(x, "AAAACCCC")
    expect_identical(paste0(stats::na.omit(x$data_wide$p1), collapse = ""),
                     unname(as.character(x$pairwise_alignments@pattern)))
    expect_equal(which(!is.na(x$data_wide$p1)), 1:8)
  }
})

test_that("different insertion lengths at the same position require repair", {
  .pwalign_test_dependencies()
  patterns <- c(p1 = "AAAATCCCC", p2 = "AAAATTCCCC")
  expect_error(.pwalign_test_run(patterns), "Overlapping indel")
  x <- .pwalign_test_run(patterns, fix_subject_indels = TRUE)
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_false(anyNA(x$data_wide$subject.position))
})

test_that("cached gaps stay with their patterns after alignment sorting", {
  .pwalign_test_dependencies()
  x <- .pwalign_test_run(c(late = "GGGGATTTT", early = "AAAATCCCC"))
  .pwalign_test_reference(x, "AAAACCCCGGGGTTTT")
  expect_identical(paste0(x$data_wide$ref, collapse = ""), "AAAA-CCCCGGGG-TTTT")
  expect_equal(which(!is.na(x$data_wide$early)), 1:9)
  expect_equal(which(!is.na(x$data_wide$late)), 10:18)
  expect_identical(paste0(stats::na.omit(x$data_wide$late), collapse = ""), "GGGGATTTT")
  y <- .pwalign_test_run(c(late = "GGGGATTTT", spanning = "AACCCC", early = "AAAATCCCC"),
                         fix_subject_indels = TRUE)
  .pwalign_test_reference(y, "AAAACCCCGGGGTTTT")
  expect_equal(which(is.na(y$data_wide$subject.position)), 13L)
  expect_equal(which(!is.na(y$data_wide$early)), 1:4)
  expect_equal(which(!is.na(y$data_wide$late)), 9:17)
})

test_that("group overlap validation rejects containment and tied starts", {
  .pwalign_test_dependencies()
  expect_error(.pwalign_test_run(list(g = c("AAAACCCCGGGG", "CCCC", "GGGG"))),
               "groups have patterns with overlapping.*g")
  expect_error(.pwalign_test_run(list(g = c("AAAA", "AAAACCCC"))),
               "groups have patterns with overlapping.*g")
  expect_no_error(.pwalign_test_run(list(g = c("GGGG", "AAAA", "CCCC"))))
})

test_that("mismatch filtering preserves first-seen length and within-group order", {
  .pwalign_test_dependencies()
  subject <- Biostrings::DNAStringSet(c(ref = "AAAACCCCGGGGTTTT"))
  patterns <- Biostrings::DNAStringSet(c(
    long_exact = "AAAACCCC", short_exact = "CCCCGG", medium_exact = "GGGTTTT",
    short_mismatch = "CCACGG", long_mismatch = "AAATCCCC"
  ))
  invisible(capture.output(x <- suppressMessages(check_for_invalid_chars(subject, patterns, 1))))
  expect_identical(names(x$pattern_mismatching_return$max_mis_0),
                   c("long_exact", "short_exact", "medium_exact"))
  expected <- c("long_exact", "long_mismatch", "short_exact", "short_mismatch", "medium_exact")
  expect_identical(names(x$pattern_mismatching_return$max_mis_1), expected)
  expect_identical(x$patterns, patterns[expected])
})
