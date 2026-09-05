#' Infer HLA type from RNA-seq reads
#'
#' Align reads from RNA-seq or single-cell RNA-seq data to known HLA reference
#' alleles, count matching reads, and rank candidate allele pairs. The plots in
#' the returned object support visual assessment of the inferred HLA type.
#'
#' Supply reference alleles for one HLA gene at a time, for example HLA-A,
#' HLA-B, or HLA-C, and run the function separately for each gene. Reads are
#' matched against every reference allele (or against a selected region such as
#' exons 2 and 3). Candidate alleles are retained according to
#' `allele_diff`, and all retained allele pairs are scored by the number of
#' reads they explain uniquely, redundantly, and in total. Pairwise results are
#' summarized by P group because typing at greater resolution may be uncertain.
#'
#' Matching is performed first with `hla_seq_colName`. If
#' `hla_seq_colName2` is present in `hla_ref`, a second round uses that sequence
#' column and the complete reference alleles associated with the leading P
#' groups from the first round. If the column is absent, the second round is
#' skipped with a message.
#'
#' @param hla_ref A data frame of HLA reference alleles, preferably created by
#'   [hla_df_from_xml()]. It must contain the sequence, allele, P-group, and
#'   G-group columns specified below.
#' @param reads A data frame of reads, preferably created by [get_bam_reads()].
#'   It must contain the read sequence and read name columns specified below.
#' @param allele_diff A numeric value greater than 1 giving the maximum fold
#'   difference in single-allele read counts used to retain alleles for pairwise
#'   matching. For example, with a maximum count of 20 and `allele_diff = 5`,
#'   alleles with at least 4 matches are retained.
#' @param top_n_pairwise_results A positive integer giving the number of leading
#'   pairwise results to include in the rank plot.
#' @param hla_seq_colName A character scalar naming the reference sequence
#'   column used for the first matching round.
#' @param hla_seq_colName2 A character scalar naming the reference sequence
#'   column used for the optional second matching round. If it is absent from
#'   `hla_ref`, the second round is skipped.
#' @param read_seq_colName A character scalar naming the sequence column in
#'   `reads`.
#' @param hla_allele_colName A character scalar naming the allele column in
#'   `hla_ref`.
#' @param read_name_colName A character scalar naming the read identifier column
#'   in `reads`.
#' @param p_group_colName A character scalar naming the P-group column in
#'   `hla_ref`.
#' @param g_group_colName A character scalar naming the G-group column in
#'   `hla_ref`.
#' @param lapply_fun A function, or the name of a function, used to apply the
#'   matching operation. Suggested values are [base::lapply()],
#'   `pbapply::pblapply`, and [parallel::mclapply()].
#' @param maxmis A non-negative integer giving the maximum number of mismatches
#'   allowed per read match.
#' @param complete_seq_columns A character vector naming columns used to identify
#'   reference alleles with complete sequence annotation for the second round.
#'   Columns not present in `hla_ref` are ignored.
#' @param make_reads_distinct A logical value indicating whether duplicated read
#'   sequences should be removed before matching.
#' @param ... Additional arguments passed to `lapply_fun`, such as `mc.cores`
#'   when [parallel::mclapply()] is used.
#'
#' @import Matrix
#'
#' @return A named list containing first- and, when available, second-round
#'   results. Each non-`NULL` round contains:
#'   \describe{
#'     \item{top_sin_res_df}{A data frame of retained alleles and their
#'       explained-read counts.}
#'     \item{top_sin_res_mat}{The read-by-allele match matrix for retained
#'       alleles.}
#'     \item{pair_res_df}{A data frame of all scored allele pairs.}
#'     \item{pair_res1_df}{Pairwise results reduced to the best result per
#'       P-group pair.}
#'     \item{pair_res2_df}{The leading pairwise results used in the rank plot.}
#'     \item{plot_sin_res}{A `ggplot` of single-allele explained-read counts.}
#'     \item{plot_pair_res1}{A `ggplot` overview of pairwise results.}
#'     \item{plot_pair_res2}{A `patchwork` object showing the leading pairwise
#'       results.}
#'   }
#'   The function returns `NULL` if no read matches any reference allele.
#' @export
#'
#' @examples
#' \dontrun{
#' hla_a_results <- hla_typing(
#'   hla_ref = subset(hla_reference, grepl("HLA-A", allele)),
#'   reads = bam_reads
#' )
#' }
hla_typing <- function(hla_ref,
                       reads,
                       allele_diff = 5,
                       top_n_pairwise_results = 50,
                       hla_seq_colName = "seq_Exon2_3",
                       hla_seq_colName2 = "seq",
                       read_seq_colName = "seq",
                       hla_allele_colName = "allele",
                       read_name_colName = "readName",
                       p_group_colName = "p_group",
                       g_group_colName = "g_group",
                       lapply_fun = lapply,
                       maxmis = 0,
                       complete_seq_columns = c("5'UTR", paste0("Exon", 1:8), "3'UTR"),
                       make_reads_distinct = F,
                       ...) {

  # check for unique read names
  # check for columns in data frames

  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }

  lapply_fun <- match.fun(lapply_fun)
  arg_list <- list(...)

  if (missing(hla_ref)) {
    stop("hla_ref missing.")
  }
  if (missing(reads)) {
    stop("reads missing.")
  }
  if (allele_diff < 1 || allele_diff == 1) {
    stop("allele_diff has to be larger than 1.")
  }
  if (top_n_pairwise_results < 1) {
    stop("top_n_pairwise_results has to be at least 1.")
  }

  if (!read_seq_colName %in% names(reads)) {
    stop(read_seq_colName, " not found in names of reads.")
  }
  if (!read_name_colName %in% names(reads)) {
    stop(read_name_colName, " not found in names of reads.")
  }
  if (!hla_seq_colName %in% names(hla_ref)) {
    stop(hla_seq_colName, " not found in names of hla_ref")
  }
  if (!hla_seq_colName2 %in% names(hla_ref)) {
    message(hla_seq_colName2, " not found in names of hla_ref. Will skip second round of read alignment.")
    hla_seq_colName2 <- NULL
  }
  if (!p_group_colName %in% names(hla_ref)) {
    stop(p_group_colName, " not found in names of hla_ref")
  }
  if (!g_group_colName %in% names(hla_ref)) {
    stop(g_group_colName, " not found in names of hla_ref")
  }
  if (!hla_allele_colName %in% names(hla_ref)) {
    stop(hla_allele_colName, " not found in names of hla_ref")
  }

  if (anyDuplicated(reads[,read_name_colName,drop=T])) {
    message("Duplicated ", read_name_colName, " found. Will make them unique.")
    reads[,read_name_colName] <- make.unique(reads[,read_name_colName,drop=T])
  }

  if (length(invalid <- which(grepl("[^ACTGU]", reads[,read_seq_colName,drop=T]))) > 0) {
    message(length(invalid), " sequences with non-DNA or non-RNA characters detected. Those are excluded.")
    reads <- reads[-invalid,]
    if (nrow(reads) == 0) {
      stop("No reads left after filtering for ones with valid DNA/RNA characters only. Please fix the sequences.")
    }
  }

  if (sum(any(grepl("HLA-A", hla_ref[[hla_allele_colName]])),
          any(grepl("HLA-B", hla_ref[[hla_allele_colName]])),
          any(grepl("HLA-C", hla_ref[[hla_allele_colName]]))) > 1) {
    message("More than one of 'HLA-A', 'HLA-B' or 'HLA-C' detected in column '", hla_allele_colName, "' of hla_ref. Did you forget to filter for one reference allele only.")
  }

  # remove columns with all NA, e.g. for HLA-B there is no exon8 --> remove the column
  hla_ref <- hla_ref[,colSums(is.na(hla_ref)) < nrow(hla_ref)]

  complete_alleles <-
    hla_ref |>
    dplyr::select(dplyr::all_of(names(hla_ref)[which(names(hla_ref) %in% complete_seq_columns)]), allele)
  complete_alleles <- complete_alleles |>
    dplyr::filter(stats::complete.cases(complete_alleles))|>
    dplyr::pull(!!rlang::sym(hla_allele_colName))

  hla_ref <- hla_ref |>
    dplyr::mutate(complete_seq = allele %in% complete_alleles) |>
    dplyr::distinct(!!rlang::sym(hla_seq_colName), .keep_all = TRUE)


  if (make_reads_distinct) {
    n_before <- nrow(reads)
    reads <- dplyr::distinct(reads, !!rlang::sym(read_seq_colName), .keep_all = T)
    n_after <- nrow(reads)
    if (n_after < n_before) {
      message(n_after, " of ", n_before, " (", round(n_after/(n_after+n_before)*100), " %) reads are unique. Will do matching only with those.")
    } else {
      message("No duplicated reads found.")
    }
  }

  first_round_results <- run_read_matching_and_report_results(hla_ref = hla_ref,
                                                              reads = reads,
                                                              allele_diff = allele_diff,
                                                              top_n_pairwise_results = top_n_pairwise_results,
                                                              hla_seq_colName = hla_seq_colName,
                                                              read_seq_colName = read_seq_colName,
                                                              hla_allele_colName = hla_allele_colName,
                                                              read_name_colName = read_name_colName,
                                                              p_group_colName = p_group_colName,
                                                              g_group_colName = g_group_colName,
                                                              lapply_fun = lapply_fun,
                                                              maxmis = maxmis,
                                                              arg_list = arg_list,
                                                              ...)
  if (is.null(first_round_results)) {
    return(NULL)
  }

  top_alleles <- colnames(first_round_results[["top_sin_res_mat"]])

  hla_ref_top <- dplyr::filter(hla_ref, !!rlang::sym(hla_allele_colName) %in% unique(top_alleles))

  hla_ref_top_complete <- hla_ref |>
    dplyr::filter(!!rlang::sym(p_group_colName) %in% unique(hla_ref_top[[p_group_colName]])) |>
    dplyr::filter(complete_seq)

  second_round_results <- NULL
  if (!is.null(hla_seq_colName2)) {
    message("Second round on sequences from ", hla_seq_colName2, " column.")
    second_round_results <- run_read_matching_and_report_results(hla_ref = hla_ref_top_complete,
                                                                 reads = reads,
                                                                 allele_diff = 1000, # arbitrary high number to make it uneffective; I did not want a second round of filtering
                                                                 top_n_pairwise_results = top_n_pairwise_results,
                                                                 hla_seq_colName = hla_seq_colName2, # diff to first round
                                                                 read_seq_colName = read_seq_colName,
                                                                 hla_allele_colName = hla_allele_colName,
                                                                 read_name_colName = read_name_colName,
                                                                 p_group_colName = p_group_colName,
                                                                 g_group_colName = g_group_colName,
                                                                 lapply_fun = lapply_fun,
                                                                 maxmis = maxmis,
                                                                 arg_list = arg_list,
                                                                 ...)
  }


  return(stats::setNames(list(first_round_results, second_round_results), c(hla_seq_colName, hla_seq_colName2)))

}


run_read_matching_and_report_results <- function(hla_ref,
                                                 reads,
                                                 allele_diff = 5,
                                                 top_n_pairwise_results = 50,
                                                 hla_seq_colName = "seq_Exon2_3",
                                                 read_seq_colName = "seq",
                                                 hla_allele_colName = "allele",
                                                 read_name_colName = "readName",
                                                 p_group_colName = "p_group",
                                                 g_group_colName = "g_group",
                                                 lapply_fun = lapply,
                                                 maxmis = 0,
                                                 arg_list = arg_list,
                                                 ...) {

  requireNamespace(Matrix) # required for sparseMatrix below; saves memory
  message("Calculating single matches.")
  if ("strand" %in% names(reads)) {
    single_res <- lapply(split(reads, as.character(reads$strand)),
                         single_matching,
                         lapply_fun = lapply_fun,
                         hla_ref = hla_ref,
                         maxmis = maxmis,
                         hla_seq_colName = hla_seq_colName,
                         hla_allele_colName = hla_allele_colName,
                         read_name_colName = read_name_colName,
                         read_seq_colName = read_seq_colName,
                         ...)
    for (i in names(single_res)) {
      message("Strand (", i, "):")
      reads_w_no_match_sum <- sum(Matrix::rowSums(single_res[[i]]) == 0)
      reads_w_min_one_match_sum <- sum(Matrix::rowSums(single_res[[i]]) > 0)
      message("  ", reads_w_min_one_match_sum, " reads with at least one match/hit, (", round(reads_w_min_one_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
      message("  ", reads_w_no_match_sum, " reads with no match/hit, (", round(reads_w_no_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
    }
    single_res <- do.call(rbind, single_res)
  } else {
    single_res <- lapply(reads,
                         single_matching,
                         lapply_fun = lapply_fun,
                         hla_ref = hla_ref,
                         maxmis = maxmis,
                         hla_seq_colName = hla_seq_colName,
                         hla_allele_colName = hla_allele_colName,
                         ...)
  }

  reads_w_min_one_match <- Matrix::rowSums(single_res) > 0
  expl_reads_per_allele <- Matrix::colSums(single_res)

  reads_w_no_match_sum <- sum(Matrix::rowSums(single_res) == 0)
  reads_w_min_one_match_sum <- sum(reads_w_min_one_match)

  message("Total:")
  if (reads_w_min_one_match_sum == 0) {
    message("No matches/hits determined.")
    return(NULL)
  }
  message("  ", reads_w_min_one_match_sum, " reads with at least one match/hit, (", round(reads_w_min_one_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
  message("  ", reads_w_no_match_sum, " reads with no match/hit, (", round(reads_w_no_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")

  # as.matrix here as this will speed up iteration over col.combs below!
  top_single_res <- as.matrix(single_res[which(reads_w_min_one_match), which(expl_reads_per_allele >= max(expl_reads_per_allele)/allele_diff)])
  top_sin_res_df <-
    data.frame(expl_reads = Matrix::colSums(top_single_res)) |>
    tibble::rownames_to_column(hla_allele_colName) |>
    dplyr::left_join(hla_ref, by = hla_allele_colName) |>
    dplyr::mutate(rank = dplyr::row_number(-expl_reads))

  #col.combs <- utils::combn(1:ncol(top_single_res), 2, simplify = F)
  # TODO: catch error when only 1 read matches
  col.combs <- t(utils::combn(1:ncol(top_single_res), 2, simplify = T))
  message("Calculating pairwise matches. Combinations: ", nrow(col.combs), ".")

  # doing this in R was too slow. other packages did not have the functionality
  # so, written in c++
  if (identical(lapply_fun, parallel::mclapply) && "mc.cores" %in% names(arg_list)) {
    # Split the matrix into chunks for multithreading
    pairwise_results <- lapply_fun(brathering::split_mat(col.combs, n_chunks = arg_list[["mc.cores"]], byrow = T), function(x) {
      igsc:::countOccurrencesInCpp(top_single_res, x)
    }, ...)
    pairwise_results <- do.call(rbind, pairwise_results)
  } else {
    pairwise_results <- igsc:::countOccurrencesInCpp(top_single_res, col.combs)
  }

  pair_res_df <-
    data.frame(uni_expl_reads = pairwise_results[,1], # sapply(pairwise_results, "[", 1),
               dbl_expl_reads = pairwise_results[,2], # sapply(pairwise_results, "[", 2),
               allele1 = colnames(top_single_res)[col.combs[,1]],
               allele2 = colnames(top_single_res)[col.combs[,2]]) |>
    dplyr::mutate(tot_expl_reads = uni_expl_reads + dbl_expl_reads, .after = dbl_expl_reads) |>
    dplyr::mutate(uni_expl_reads_rank = dplyr::dense_rank(-uni_expl_reads),
                  tot_expl_reads_rank = dplyr::dense_rank(-tot_expl_reads)) |>
    dplyr::mutate(rank_sum = (tot_expl_reads_rank + uni_expl_reads_rank)/2)
  pair_res_df <- pair_res_df |>
    dplyr::mutate(allele_comb = scexpr:::orderAndConcatenateStrings(as.matrix(pair_res_df[c("allele1", "allele2")])), .after = "allele2")|>
    #dplyr::mutate(rank_int = dplyr::dense_rank(base::interaction(-tot_expl_reads, -uni_expl_reads, lex.order = TRUE)))|>
    dplyr::group_by(allele_comb)|>
    #dplyr::filter(rank_int == min(rank_int))|>
    dplyr::filter(rank_sum == min(rank_sum))|>
    dplyr::ungroup()|>
    dplyr::distinct(allele_comb, .keep_all = T)|>
    dplyr::left_join(hla_ref|> dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele1" = hla_allele_colName))|>
    dplyr::rename("p_group1" = p_group_colName, "g_group1" = g_group_colName)|>
    dplyr::left_join(hla_ref|> dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele2" = hla_allele_colName))|>
    dplyr::rename("p_group2" = p_group_colName, "g_group2" = g_group_colName)|>
    dplyr::left_join(top_sin_res_df|> dplyr::select(dplyr::all_of(hla_allele_colName), expl_reads), by = c("allele1" = hla_allele_colName))|>
    dplyr::rename("allele1_expl_reads" = expl_reads)|>
    dplyr::left_join(top_sin_res_df|> dplyr::select(dplyr::all_of(hla_allele_colName), expl_reads), by = c("allele2" = hla_allele_colName))|>
    dplyr::rename("allele2_expl_reads" = expl_reads)|>
    dplyr::mutate(expl_reads_overall = !!reads_w_min_one_match_sum)|>
    dplyr::mutate(non_expl_reads_overall = !!reads_w_no_match_sum)|>
    dplyr::mutate(allele12_expl_read_diff = abs(allele1_expl_reads - allele2_expl_reads))|>
    dplyr::mutate(frac_match_reads_expl = tot_expl_reads/expl_reads_overall)|>
    dplyr::mutate(frac_all_reads_expl = tot_expl_reads/(expl_reads_overall + non_expl_reads_overall))|>
    dplyr::mutate(p_group12 = paste0(p_group1, "_", p_group2))|>
    dplyr::mutate(g_group12 = paste0(g_group1, "_", g_group2))|>
    dplyr::mutate(allele_group1 = stringr::str_extract(allele1, "[:alpha:]\\*[:digit:]{2}"))|>
    dplyr::mutate(allele_group2 = stringr::str_extract(allele2, "[:alpha:]\\*[:digit:]{2}"))|>
    dplyr::mutate(allele_group12 = paste0(allele_group1, "_", allele_group2))

  ## plotting
  top_pair_res_plot1 <-
    #top_pair_res_df
    pair_res_df|>
    dplyr::group_by(p_group12)|>
    dplyr::slice_min(order_by = rank_sum, n = 1)|>
    dplyr::ungroup()

  allele_group12_medians <-
    top_pair_res_plot1|>
    dplyr::group_by(allele_group12)|>
    dplyr::summarise(median_frac_match_reads_expl = stats::median(frac_match_reads_expl))|>
    dplyr::arrange(dplyr::desc(median_frac_match_reads_expl))

  overview_plot <- ggplot2::ggplot(top_pair_res_plot1, ggplot2::aes(x = brathering::reorder_within(allele_group2, frac_match_reads_expl, allele_group1), y = frac_match_reads_expl)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    brathering::scale_x_reordered() +
    ggplot2::labs(x = "allele_group2") +
    ggplot2::geom_hline(yintercept = allele_group12_medians[1,2,drop=T], color = "tomato2") +
    ggplot2::geom_hline(yintercept = max(top_pair_res_plot1$frac_match_reads_expl), color = "forestgreen") +
    ggplot2::facet_wrap(ggplot2::vars(allele_group1), nrow = 1, scales = "free_x")

  top_pair_res_plot2 <-
    top_pair_res_plot1|>
    dplyr::mutate(rank_plot = dplyr::row_number(-tot_expl_reads))|>
    dplyr::slice_min(order_by = rank_plot, n = top_n_pairwise_results)|>
    dplyr::arrange(rank_plot)|>
    dplyr::mutate(plot.color = as.factor(ifelse(rank_plot %% 2 != 0, 1, 2)))

  sin_plot <- ggplot2::ggplot(top_sin_res_df, ggplot2::aes(x = brathering::reorder_within(!!rlang::sym(hla_allele_colName), expl_reads, allele_group), y = expl_reads)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("allele") +
    ggplot2::ylab("n explained reads") +
    ggplot2::theme_bw() +
    brathering::scale_x_reordered() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::facet_wrap(ggplot2::vars(allele_group), scales = "free_x")

  rank.plot.p1 <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(rank_plot), y = stats::reorder(p_group1, rank_plot), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 1") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.plot.p2 <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(rank_plot), y = stats::reorder(p_group2, rank_plot), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 2") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.read.plot <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(rank_plot), y = tot_expl_reads, fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("total\nexplained\nreads") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10))

  expl_reads_overall <- as.numeric(levels(as.factor(top_pair_res_plot2$expl_reads_overall)))
  max_tot_reads <- max(top_pair_res_plot2$tot_expl_reads)
  min_tot_reads <- min(top_pair_res_plot2$tot_expl_reads)
  max_min_diff <- max_tot_reads - min_tot_reads
  if (max_min_diff > 10) {
    breaks <- seq(max_tot_reads,
                  min_tot_reads,
                  by = -max(c(10, brathering::floor2((max_tot_reads - min_tot_reads)/3, 10))))
  } else if (max_min_diff > 5) {
    breaks <- seq(max_tot_reads, min_tot_reads, by = -5)
  } else {
    breaks <- seq(max_tot_reads, min_tot_reads, by = -1)
  }
  rank.read.plot <-
    rank.read.plot +
    ggplot2::scale_y_continuous(breaks = breaks,
                                sec.axis = ggplot2::sec_axis(~ . / expl_reads_overall * 100, name = "% of matching\nreads that\nmatched min. one\nref. allele"))

  height.1 <- nlevels(as.factor(top_pair_res_plot2$p_group1))
  height.2 <- nlevels(as.factor(top_pair_res_plot2$p_group2))
  if (height.1 / height.2 < 0.1) {
    height.1 <- 10
    height.2 <- 100 - height.1
  }
  if (height.2 / height.1 < 0.1) {
    height.2 <- 10
    height.1<- 100 - height.2
  }
  sum.1.2 <- height.1 + height.2
  height.3 <- sum.1.2*0.15
  height.1 = (sum.1.2 - height.3)*height.1/sum.1.2
  height.2 = (sum.1.2 - height.3)*height.2/sum.1.2
  total <- height.1 + height.2 + height.3
  height.1 <- height.1/total
  height.2 <- height.2/total
  height.3 <- height.3/total

  #pair_plot <- cowplot::plot_grid(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, align = "v", rel_heights = c(height.1,height.2,height.3)) # check how to replace with patchwork
  pair_plot <- patchwork::wrap_plots(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, heights = c(height.1,height.2,height.3))

  return(list(top_sin_res_df = top_sin_res_df,
              top_sin_res_mat = top_single_res,
              pair_res_df = pair_res_df,
              pair_res1_df = top_pair_res_plot1,
              pair_res2_df = top_pair_res_plot2,
              plot_sin_res = sin_plot,
              plot_pair_res1 = overview_plot,
              plot_pair_res2 = pair_plot))
}



single_matching <- function(reads,
                            lapply_fun,
                            hla_ref,
                            maxmis,
                            hla_seq_colName,
                            hla_allele_colName,
                            read_name_colName,
                            read_seq_colName,
                            ...) {

  reads <- stats::setNames(reads[,read_seq_colName,drop=T], reads[,read_name_colName,drop=T])
  single_res <- do.call(rbind, lapply_fun(split(reads, ceiling(seq_along(reads)/1e3)), function(rrr) {

    # vwhichPDict allows for max.mismatch, but vcountPDict does not
    hit_matrix <- Biostrings::vwhichPDict(subject = Biostrings::DNAStringSet(hla_ref[,hla_seq_colName,drop=T]),
                                          pdict = Biostrings::PDict(rrr, max.mismatch = maxmis),
                                          max.mismatch = maxmis)
    hit_matrix <- lapply(hit_matrix, function(z) replace(rep(0, length(rrr)), z, 1))
    hit_matrix <- methods::as(do.call(cbind, hit_matrix), "sparseMatrix")
    colnames(hit_matrix) <- hla_ref[,hla_allele_colName,drop=T]
    rownames(hit_matrix) <- names(rrr)
    return(hit_matrix)
  }, ...))
  return(single_res)
}

