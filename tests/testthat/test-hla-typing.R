test_that("print_single_res reports matches by strand and returns counts", {
  skip_if_not_installed("Matrix")

  single_res <- Matrix::sparseMatrix(
    i = c(1L, 1L, 3L),
    j = c(1L, 2L, 2L),
    x = 1L,
    dims = c(3L, 2L)
  )
  reads <- data.frame(strand = c("+", NA_character_, "-"))
  messages <- character()

  result <- withCallingHandlers(
    igsc:::print_single_res(single_res, reads),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )

  expect_identical(result$reads_w_min_one_match, c(TRUE, FALSE, TRUE))
  expect_equal(result$expl_reads_per_allele, c(1, 2))
  expect_equal(result$reads_w_min_one_match_sum, 2)
  expect_equal(result$reads_w_no_match_sum, 1)
  expect_match(paste(messages, collapse = "\n"), "Strand \\(\\+\\):")
  expect_match(paste(messages, collapse = "\n"), "Strand \\(<NA>\\):")
  expect_match(paste(messages, collapse = "\n"), "Strand \\(-\\):")
  expect_match(paste(messages, collapse = "\n"), "Total:")
})

test_that("print_single_res returns NULL when no reads match", {
  skip_if_not_installed("Matrix")

  single_res <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = integer(),
    dims = c(2L, 2L)
  )
  reads <- data.frame(qname = c("r1", "r2"))

  expect_null(suppressMessages(igsc:::print_single_res(single_res, reads)))
})

test_that("pairwise_matching returns allele summaries and pair counts", {
  skip_if_not_installed("Matrix")

  allele_names <- c("HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01")
  top_single_res <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L, 3L, 1L, 3L),
    j = c(1L, 1L, 2L, 2L, 3L, 3L),
    x = 1L,
    dims = c(3L, 3L),
    dimnames = list(c("r1", "r2", "r3"), allele_names)
  )
  attr(top_single_res, "read_weights") <- c(2L, 3L, 5L)
  hla_ref <- data.frame(allele = allele_names)

  result <- suppressMessages(igsc:::pairwise_matching(
    top_single_res,
    hla_ref,
    hla_allele_col_name = "allele",
    lapply_fun = lapply,
    arg_list = list()
  ))

  expect_setequal(result$groups$allele, allele_names)
  expect_setequal(result$top_sin_res_df$allele, allele_names)
  expect_equal(result$top_sin_res_df$expl_reads, c(5, 8, 7))
  expect_equal(nrow(result$col.combs), 3)
  expected_counts <- t(apply(result$col.combs, 1L, function(cols) {
    first <- as.logical(top_single_res[, cols[1L]])
    second <- as.logical(top_single_res[, cols[2L]])
    weights <- attr(top_single_res, "read_weights")
    c(sum(weights[xor(first, second)]), sum(weights[first & second]))
  }))
  expect_equal(result$pairwise_results, expected_counts)
})

test_that("single_matching matches unique sequences and records duplicate weights", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("brathering")

  reads <- data.frame(
    qname = paste0("r", seq_len(9L)),
    seq = c("AAC", "GGG", "AAC", "CCC", "GGG",
            "AACT", "AACT", "TGT", "TGT")
  )
  hla_ref <- data.frame(
    allele = c("HLA-A*01:01", "HLA-A*02:01"),
    seq_Exon2_3 = c("AACTTT", "GGGCCC")
  )
  matched_indices <- integer()
  recording_apply <- function(X, FUN, ...) {
    matched_indices <<- unlist(X, use.names = FALSE)
    lapply(X, FUN, ...)
  }

  result <- igsc:::single_matching(
    reads = reads,
    lapply_fun = recording_apply,
    hla_ref = hla_ref,
    maxmis = 0L,
    hla_seq_col_name = "seq_Exon2_3",
    hla_allele_col_name = "allele",
    read_name_col_name = "qname",
    read_seq_col_name = "seq"
  )

  expected_unique <- matrix(
    0,
    nrow = length(unique(reads$seq)),
    ncol = nrow(hla_ref),
    dimnames = list(c("r1", "r2", "r4", "r6", "r8"), hla_ref$allele)
  )
  expected_unique[c(1L, 4L), 1L] <- 1
  expected_unique[c(2L, 3L), 2L] <- 1

  expect_s4_class(result, "dgCMatrix")
  expect_equal(as.matrix(result), expected_unique)
  expect_equal(attr(result, "read_weights"), c(2L, 2L, 1L, 2L, 2L))
  expect_equal(attr(result, "read_to_unique"), c(1L, 2L, 1L, 3L, 2L,
                                                  4L, 4L, 5L, 5L))
  stats <- suppressMessages(igsc:::print_single_res(result, reads))
  expect_equal(unname(stats$expl_reads_per_allele), c(4, 3))
  expect_equal(unname(stats$reads_w_min_one_match),
               c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_setequal(matched_indices, seq_along(unique(reads$seq)))
  expect_length(matched_indices, length(unique(reads$seq)))

  distinct_reads <- reads[!duplicated(reads$seq), , drop = FALSE]
  distinct_result <- igsc:::single_matching(
    reads = distinct_reads,
    lapply_fun = lapply,
    hla_ref = hla_ref,
    maxmis = 0L,
    hla_seq_col_name = "seq_Exon2_3",
    hla_allele_col_name = "allele",
    read_name_col_name = "qname",
    read_seq_col_name = "seq"
  )
  expect_equal(
    as.matrix(distinct_result),
    expected_unique
  )
  expect_equal(attr(distinct_result, "read_weights"), rep.int(1L, 5L))
})

test_that("heuristic pairwise matching shortlists pairs and agrees when unbounded", {
  skip_if_not_installed("Matrix")

  allele_names <- sprintf("HLA-A*%02d:01", seq_len(6L))
  read_sets <- list(1:6, 5:10, 1:3, 7:9, 2:4, 8:10)
  hit_rows <- unlist(read_sets, use.names = FALSE)
  hit_cols <- rep.int(seq_along(read_sets), lengths(read_sets))
  top_single_res <- Matrix::sparseMatrix(
    i = hit_rows,
    j = hit_cols,
    x = 1L,
    dims = c(10L, 6L),
    dimnames = list(paste0("r", seq_len(10L)), allele_names)
  )
  hla_ref <- data.frame(allele = allele_names)

  exhaustive <- suppressMessages(igsc:::pairwise_matching(
    top_single_res, hla_ref, "allele", lapply, list()
  ))
  full_shortlist <- suppressMessages(igsc:::pairwise_matching_heuristic(
    top_single_res, hla_ref, "allele", lapply, list(),
    anchor_n = 6L, partner_n = 5L
  ))
  short_shortlist <- suppressMessages(igsc:::pairwise_matching_heuristic(
    top_single_res, hla_ref, "allele", lapply, list(),
    anchor_n = 1L, partner_n = 2L
  ))

  pair_keys <- function(result) {
    paste(result$col.combs[, 1], result$col.combs[, 2], sep = ":")
  }
  expect_setequal(pair_keys(full_shortlist), pair_keys(exhaustive))
  expect_equal(
    full_shortlist$pairwise_results[
      match(pair_keys(exhaustive), pair_keys(full_shortlist)),
      ,
      drop = FALSE
    ],
    exhaustive$pairwise_results
  )
  expect_equal(nrow(short_shortlist$col.combs), 2)
  expect_equal(short_shortlist$col.combs, rbind(c(1L, 2L), c(1L, 4L)))
  expect_equal(
    short_shortlist$pairwise_results,
    igsc:::countOccurrencesSparseCpp(top_single_res, short_shortlist$col.combs)
  )
})

test_that("hla_typing dispatches between exhaustive and heuristic pairwise search", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("brathering")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  alleles <- sprintf("HLA-A*%02d:01", seq_len(4L))
  hla_ref <- data.frame(
    allele = alleles,
    seq_Exon2_3 = c("AAACCC", "GGGTTT", "CCCGGG", "TTTAAA"),
    p_group = paste0("A*", sprintf("%02d", seq_len(4L)), "P"),
    g_group = paste0("A*", sprintf("%02d", seq_len(4L)), "G")
  )
  reads <- data.frame(
    qname = paste0("r", seq_len(8L)),
    seq = rep(c("AAA", "GGG", "CCC", "TTT"), 2L)
  )

  exhaustive <- suppressMessages(hla_typing(hla_ref, reads))
  heuristic <- suppressMessages(hla_typing(
    hla_ref,
    reads,
    pairwise_method = "heuristic",
    pairwise_anchor_n = 1L,
    pairwise_partner_n = 1L
  ))

  expect_equal(nrow(exhaustive$pair_res_df), 6)
  expect_equal(nrow(heuristic$pair_res_df), 1)
  expect_equal(
    heuristic$top_sin_res_df$expl_reads,
    exhaustive$top_sin_res_df$expl_reads
  )
})
