.aln_test_dependencies <- function() {
  for (p in c("Biostrings", "pwalign", "brathering", "colrr", "ggplot2", "Peptides",
              "RColorBrewer", "crayon", "ggrepel", "scales", "collapse")) skip_if_not_installed(p)
}
.aln_test_data <- function() {
  data.frame(position = rep(1:10, 2), seq.name = rep(c("zref", "apat"), each = 10),
             seq = c(rep("A", 10), rep(NA, 3), "C", "C", rep(NA, 5)))
}
.aln_test_plot <- function(aln = .aln_test_data(), ...) {
  suppressMessages(aln_plot(aln, subject_name = "zref", aln_type = "NT", verbose = FALSE, ...))
}

test_that("AA presets color the expected residues", {
  .aln_test_dependencies()
  x <- data.frame(position = 1:4, seq.name = "ref", seq = c("A", "K", "E", "F"))
  for (scheme in c("Chemistry_AA", "Shapely_AA", "Zappo_AA", "Taylor_AA")) {
    expect_no_error(p <- aln_plot(x, aln_type = "AA", subject_name = "ref", tile_fill = scheme))
    expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill,
                     unname(igsc:::scheme_AA[x$seq, scheme]))
  }
  expect_no_error(aln_plot(x, aln_type = "AA", subject_name = "ref"))
  x <- rbind(x, transform(x, seq.name = "pat"))
  expect_warning(p <- aln_plot(x, aln_type = "AA", subject_name = "ref", ref = "ref"),
                 "change_nonref and change_ref not functioning currently")
  expect_identical(x$seq, rep(c("A", "K", "E", "F"), 2))
  expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill,
                   rep(unname(igsc:::scheme_AA[c("A", "K", "E", "F"), "Chemistry_AA"]), 2))
})

test_that("pairwise objects and lists convert using their sequence type", {
  .aln_test_dependencies()
  pa <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(c(pat = "AAAA")),
                                  Biostrings::DNAStringSet(c(ref = "AAAACCCC")), type = "global-local")
  for (input in list(pa, list(pa))) {
    expect_no_error(p <- aln_plot(input, verbose = FALSE))
    expect_setequal(as.character(p$data$seq.name), c("ref", "pat"))
    expect_no_warning(ggplot2::ggplot_build(p))
  }
  pa2 <- pwalign::pairwiseAlignment(Biostrings::AAStringSet(c(pat = "AKEF")),
                                   Biostrings::AAStringSet(c(ref = "AKEF")))
  expect_no_error(ggplot2::ggplot_build(aln_plot(pa2, verbose = FALSE)))
  expect_error(aln_plot(list(), verbose = FALSE), "non-empty list")
})

test_that("multiple pairwise alignments retain all pattern rows", {
  .aln_test_dependencies()
  pa <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(c(p1 = "AAAA", p2 = "CCCC")),
                                  Biostrings::DNAStringSet(c(ref = "AAAACCCC")), type = "global-local")
  for (input in list(pa, list(pa[1], pa[2]))) {
    p <- aln_plot(input, verbose = FALSE)
    expect_setequal(as.character(p$data$seq.name), c("p1", "p2", "ref"))
    expect_equal(sum(p$data$seq.name == "ref"), 8)
  }
})

test_that("circular shifting preserves coordinates including no-op shifts", {
  expect_identical(shifted_pos(1:10, start_pos = 1, verbose = FALSE), 1:10)
  expect_identical(shifted_pos(1:10, start_pos = 4, verbose = FALSE), c(8:10, 1:7))
  expect_identical(shifted_pos(5L, start_pos = 5, verbose = FALSE), 5L)
  expect_identical(shifted_pos(1:10, n = 10, verbose = FALSE), 1:10)
  expect_error(shifted_pos(1:10, start_pos = 99, verbose = FALSE), "not found")
})

test_that("shifted tiles and original-coordinate labels agree", {
  .aln_test_dependencies()
  for (shift in list(1, 4, "+1", "+3")) {
    p <- .aln_test_plot(pos_shift = shift, x_breaks = 1:10)
    expect_equal(sort(unique(p$data$position)), 1:10)
    expect_false(anyNA(p$data$position))
    axis <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x
    labels <- as.numeric(axis$get_labels())
    positions <- p$data$position[1:10]
    expect_equal(labels[match(positions, axis$get_breaks())], 1:10)
  }
  a <- .aln_test_plot(pos_shift = "+3")
  b <- .aln_test_plot(pos_shift = 9)
  expect_identical(a$data$position, b$data$position)
  p <- .aln_test_plot(pos_shift = 4, pos_shift_adjust_axis = FALSE, x_breaks = 1:10)
  expect_equal(as.numeric(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()), 1:10)
})

test_that("shifting maps offset and sparse positions by value", {
  .aln_test_dependencies()
  x <- .aln_test_data(); x$position <- x$position + 100
  p <- .aln_test_plot(x, pos_shift = 105)
  expect_equal(sort(unique(p$data$position)), 101:110)
  expect_equal(p$data$position[5], 101)
  x <- x[x$position %in% c(101, 105, 110), ]
  p <- .aln_test_plot(x, pos_shift = 105)
  expect_equal(sort(unique(p$data$position)), c(101, 105, 110))
  expect_false(anyNA(p$data$position))
  for (bad in list(NA, Inf, c(1, 2), "+bad", "+0", "+11")) {
    expect_error(.aln_test_plot(pos_shift = bad), "pos_shift")
  }
})

test_that("length labels use sequence names even for reordered factors", {
  .aln_test_dependencies()
  x <- .aln_test_data(); x$seq.name <- factor(x$seq.name, levels = c("apat", "zref"))
  p <- .aln_test_plot(x, add_length_suffix = TRUE)
  axis <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y
  expect_identical(as.character(axis$get_labels()), c("apat\n2 nt", "zref\n10 nt"))
  pa <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(c(apat = "CC")),
                                  Biostrings::DNAStringSet(c(zref = "AAACCAAAAA")), type = "global-local")
  p <- .aln_test_plot(x, add_length_suffix = TRUE, pairwise_alignment = pa)
  expect_identical(as.character(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()),
                   c("apat\n2 nt", "zref\n10 nt"))
})

test_that("focus ignores NA padding and handles empty patterns", {
  .aln_test_dependencies()
  p <- .aln_test_plot(focus = 1)
  expect_equal(range(p$data$position), c(3, 6))
  x <- tibble::as_tibble(.aln_test_data())
  expect_equal(range(.aln_test_plot(x, focus = 0)$data$position), c(4, 5))
  x$seq[x$seq.name == "apat"] <- NA
  expect_no_warning(p <- .aln_test_plot(x, focus = 1))
  expect_equal(range(p$data$position), c(1, 10))
  for (bad in list(-1, NA, Inf, c(1, 2))) expect_error(.aln_test_plot(focus = bad), "focus")
})

test_that("custom column names support grouping, ordering, focus, and limits", {
  .aln_test_dependencies()
  x <- rbind(.aln_test_data(), data.frame(position = 8:9, seq.name = "bpat", seq = "G"))
  names(x) <- c("pos", "id", "base")
  p <- .aln_test_plot(x, name_col = "id", pos_col = "pos", seq_col = "base",
                      group_on_yaxis = TRUE, min_gap = 1)
  expect_true("group" %in% names(p$data))
  expect_equal(as.character(p$data$group[p$data$id == "apat"])[1],
               as.character(p$data$group[p$data$id == "bpat"])[1])
  p <- .aln_test_plot(x, name_col = "id", pos_col = "pos", seq_col = "base", y_order = "decreasing")
  expect_identical(levels(p$data$id), c("zref", "bpat", "apat"))
  p <- .aln_test_plot(x, name_col = "id", pos_col = "pos", seq_col = "base", focus = 0)
  expect_equal(range(p$data$pos), c(4, 9))
  pa <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(c(apat = "CC")),
                                  Biostrings::DNAStringSet(c(zref = "AAACCAAAAA")), type = "global-local")
  expect_no_warning(ggplot2::ggplot_build(.aln_test_plot(
    x, name_col = "id", pos_col = "pos", seq_col = "base", pairwise_alignment = pa, pattern_lim_size = 2)))
  p <- aln_plot(pa, name_col = "id", pos_col = "pos", seq_col = "base", verbose = FALSE)
  expect_true(all(c("id", "pos", "base") %in% names(p$data)))
})

test_that("numeric labels sort numerically without changing their spelling", {
  .aln_test_dependencies()
  x <- data.frame(position = 1:4, seq.name = c("10", "02", "1", "2"), seq = "A")
  p <- aln_plot(x, subject_name = "1", aln_type = "NT", order_numeric_seq_names = TRUE)
  expect_identical(levels(p$data$seq.name), c("1", "02", "2", "10"))
  expect_false(anyNA(p$data$seq.name))
})
