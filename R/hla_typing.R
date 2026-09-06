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
#' @param hla_seq_col_name A character scalar naming the reference sequence
#'   column in `hla_ref`.
#' @param read_seq_col_name A character scalar naming the sequence column in
#'   `reads`.
#' @param hla_allele_col_name A character scalar naming the allele column in
#'   `hla_ref`.
#' @param read_name_col_name A character scalar naming the read identifier column
#'   in `reads`.
#' @param p_group_col_name A character scalar naming the P-group column in
#'   `hla_ref`.
#' @param g_group_col_name A character scalar naming the G-group column in
#'   `hla_ref`.
#' @param lapply_fun A function, or the name of a function, used to apply the
#'   matching operation. Suggested values are [base::lapply()],
#'   `pbapply::pblapply`, and [parallel::mclapply()].
#' @param maxmis A non-negative integer giving the maximum number of mismatches
#'   allowed per read match.
#' @param make_reads_distinct A logical value indicating whether duplicated read
#'   sequences should be removed before matching.
#' @param rev_comp_minus A logical value indicating whether sequences whose
#'   strand equals `minus_strand_value` should be reverse-complemented before
#'   matching. This is useful when minus-strand reads are stored in their
#'   sequenced orientation rather than reference orientation.
#' @param strand_col_name A character scalar naming the strand column in
#'   `reads`. The column is required when `rev_comp_minus = TRUE` and is used
#'   for strand-specific match summaries when present.
#' @param minus_strand_value The value in `strand_col_name` identifying
#'   minus-strand reads.
#' @param ... Additional arguments passed to `lapply_fun`, such as `mc.cores`
#'   when [parallel::mclapply()] is used.
#'
#' @import Matrix
#'
#' @return `NULL` if no read matches a reference allele; otherwise, a list
#'   containing:
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
#' # get hla refs
#' hla_ref <- hla_df_from_xml("/Volumes/CMS_SSD_2TB/hla.xml.gz",
#'                            lapply_fun = parallel::mclapply, mc.cores = 8)
#' # create synthetic reads
#' reads_hla <- simulate_hla_reads(hla_ref, snps_per_read = c(1),
#'                                 gene = "A", minus_strand_prob = 0)
#' # run typing algorithm
#' type <- hla_typing(hla_ref = hla_ref |> dplyr::filter(gene == "A"),
#'                    reads = reads_hla, maxmis = 1)
#' }
hla_typing <- function(hla_ref,
                       reads,
                       allele_diff = 5,
                       top_n_pairwise_results = 50,
                       hla_seq_col_name = "seq_Exon2_3",
                       read_seq_col_name = "seq",
                       hla_allele_col_name = "allele",
                       read_name_col_name = "readName",
                       p_group_col_name = "p_group",
                       g_group_col_name = "g_group",
                       lapply_fun = lapply,
                       maxmis = 0,
                       make_reads_distinct = FALSE,
                       rev_comp_minus = FALSE,
                       strand_col_name = "strand",
                       minus_strand_value = "-",
                       ...) {
  required_packages <- c(
    "Biostrings", "Matrix", "brathering", "dplyr", "ggplot2",
    "patchwork", "rlang", "stringr", "tibble"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop("Required packages are not installed: ",
         paste(missing_packages, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.data.frame(hla_ref) || nrow(hla_ref) == 0L) {
    stop("hla_ref must be a non-empty data frame.")
  }
  if (!is.data.frame(reads) || nrow(reads) == 0L) {
    stop("reads must be a non-empty data frame.")
  }
  if (length(allele_diff) != 1L || !is.numeric(allele_diff) ||
      is.na(allele_diff) || !is.finite(allele_diff) || allele_diff <= 1) {
    stop("allele_diff must be one finite numeric value greater than 1.")
  }
  if (length(top_n_pairwise_results) != 1L ||
      !is.numeric(top_n_pairwise_results) ||
      is.na(top_n_pairwise_results) ||
      !is.finite(top_n_pairwise_results) ||
      top_n_pairwise_results < 1L ||
      top_n_pairwise_results != floor(top_n_pairwise_results)) {
    stop("top_n_pairwise_results must be a positive integer.")
  }
  if (length(maxmis) != 1L || !is.numeric(maxmis) || is.na(maxmis) ||
      !is.finite(maxmis) || maxmis < 0L || maxmis != floor(maxmis)) {
    stop("maxmis must be a non-negative integer.")
  }
  if (length(make_reads_distinct) != 1L || is.na(make_reads_distinct) ||
      !is.logical(make_reads_distinct)) {
    stop("make_reads_distinct must be TRUE or FALSE.")
  }
  if (length(rev_comp_minus) != 1L || is.na(rev_comp_minus) ||
      !is.logical(rev_comp_minus)) {
    stop("rev_comp_minus must be TRUE or FALSE.")
  }

  column_arguments <- list(
    hla_seq_col_name,
    hla_allele_col_name,
    p_group_col_name,
    g_group_col_name,
    read_seq_col_name,
    read_name_col_name,
    strand_col_name
  )
  valid_column_arguments <- vapply(
    column_arguments,
    function(x) {
      is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    },
    logical(1)
  )
  if (!all(valid_column_arguments)) {
    stop("Column-name arguments must be non-empty character scalars.")
  }
  if (length(minus_strand_value) != 1L || is.na(minus_strand_value)) {
    stop("minus_strand_value must be one non-missing value.")
  }

  missing_hla_columns <- setdiff(
    c(hla_seq_col_name, hla_allele_col_name,
      p_group_col_name, g_group_col_name),
    names(hla_ref)
  )
  missing_read_columns <- setdiff(
    c(read_seq_col_name, read_name_col_name),
    names(reads)
  )
  if (length(missing_hla_columns) > 0L) {
    stop("Missing hla_ref columns: ", paste(missing_hla_columns, collapse = ", "))
  }
  if (length(missing_read_columns) > 0L) {
    stop("Missing reads columns: ", paste(missing_read_columns, collapse = ", "))
  }
  if (rev_comp_minus && !strand_col_name %in% names(reads)) {
    stop("rev_comp_minus = TRUE requires the strand column '",
         strand_col_name, "'.")
  }

  lapply_fun <- match.fun(lapply_fun)
  arg_list <- list(...)

  reads[[read_name_col_name]] <- as.character(reads[[read_name_col_name]])
  if (anyNA(reads[[read_name_col_name]]) ||
      any(reads[[read_name_col_name]] == "")) {
    stop("Read identifiers must not be missing or empty.")
  }
  if (anyDuplicated(reads[[read_name_col_name]])) {
    message("Duplicated ", read_name_col_name, " found. Will make them unique.")
    reads[[read_name_col_name]] <- make.unique(reads[[read_name_col_name]])
  }

  reads[[read_seq_col_name]] <- toupper(as.character(reads[[read_seq_col_name]]))
  reads[[read_seq_col_name]] <- chartr("U", "T", reads[[read_seq_col_name]])
  invalid_reads <- is.na(reads[[read_seq_col_name]]) |
    !grepl("^[ACGT]+$", reads[[read_seq_col_name]])
  if (any(invalid_reads)) {
    message(sum(invalid_reads),
            " reads with missing or invalid sequences were excluded.")
    reads <- reads[!invalid_reads, , drop = FALSE]
    if (nrow(reads) == 0) {
      stop("No reads remain after sequence validation.")
    }
  }

  if (rev_comp_minus) {
    is_minus <- !is.na(reads[[strand_col_name]]) &
      reads[[strand_col_name]] == minus_strand_value
    if (any(is_minus)) {
      reads[[read_seq_col_name]][is_minus] <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(reads[[read_seq_col_name]][is_minus])
        )
      )
      message(sum(is_minus), " minus-strand reads were reverse-complemented.")
    }
  }

  hla_ref[[hla_allele_col_name]] <- as.character(hla_ref[[hla_allele_col_name]])
  hla_ref[[hla_seq_col_name]] <- toupper(as.character(hla_ref[[hla_seq_col_name]]))
  hla_ref[[hla_seq_col_name]] <- chartr("U", "T", hla_ref[[hla_seq_col_name]])
  invalid_refs <- is.na(hla_ref[[hla_allele_col_name]]) |
    hla_ref[[hla_allele_col_name]] == "" |
    is.na(hla_ref[[hla_seq_col_name]]) |
    !grepl("^[ACGT]+$", hla_ref[[hla_seq_col_name]])
  if (any(invalid_refs)) {
    message(sum(invalid_refs),
            " reference rows with missing or invalid sequences were excluded.")
    hla_ref <- hla_ref[!invalid_refs, , drop = FALSE]
  }
  if (nrow(hla_ref) < 2L) {
    stop("At least two valid reference alleles are required.")
  }
  hla_ref <- hla_ref[
    !duplicated(hla_ref[c(hla_allele_col_name, hla_seq_col_name)]),
    ,
    drop = FALSE
  ]
  if (anyDuplicated(hla_ref[[hla_allele_col_name]])) {
    stop("hla_ref contains duplicated allele identifiers with different sequences.")
  }


  allele_genes <- sapply(strsplit(sapply(strsplit(hla_ref$allele, "-"), "[", 2), "\\*"), "[", 1)
  genes <- table(allele_genes)
  if (length(genes) > 1L) {
    stop("hla_ref must contain one HLA gene; detected: ",
         paste(names(genes), genes, sep = "=", collapse = ", "), ".")
  }


  if (make_reads_distinct) {
    n_before <- nrow(reads)
    reads <- dplyr::distinct(
      reads,
      !!rlang::sym(read_seq_col_name),
      .keep_all = TRUE
    )
    n_after <- nrow(reads)
    if (n_after < n_before) {
      message(n_after, " of ", n_before, " (",
              round(n_after / n_before * 100),
              " %) reads are unique. Matching will use only those reads.")
    } else {
      message("No duplicated reads found.")
    }
  }

  first_round_results <- run_read_matching_and_report_results(hla_ref = hla_ref,
                                                              reads = reads,
                                                              allele_diff = allele_diff,
                                                              top_n_pairwise_results = top_n_pairwise_results,
                                                              hla_seq_col_name = hla_seq_col_name,
                                                              read_seq_col_name = read_seq_col_name,
                                                              hla_allele_col_name = hla_allele_col_name,
                                                              read_name_col_name = read_name_col_name,
                                                              p_group_col_name = p_group_col_name,
                                                              g_group_col_name = g_group_col_name,
                                                              lapply_fun = lapply_fun,
                                                              maxmis = maxmis,
                                                              arg_list = arg_list,
                                                              strand_col_name = strand_col_name,
                                                              ...)
  if (is.null(first_round_results)) {
    return(NULL)
  }
  first_round_results

}


run_read_matching_and_report_results <- function(hla_ref,
                                                 reads,
                                                 allele_diff = 5,
                                                 top_n_pairwise_results = 50,
                                                 hla_seq_col_name = "seq_Exon2_3",
                                                 read_seq_col_name = "seq",
                                                 hla_allele_col_name = "allele",
                                                 read_name_col_name = "readName",
                                                 p_group_col_name = "p_group",
                                                 g_group_col_name = "g_group",
                                                 lapply_fun = lapply,
                                                 maxmis = 0,
                                                 arg_list = list(),
                                                 strand_col_name = "strand",
                                                 ...) {

  message("Calculating single matches.")
  single_res <- single_matching(
    reads = reads,
    lapply_fun = lapply_fun,
    hla_ref = hla_ref,
    maxmis = maxmis,
    hla_seq_col_name = hla_seq_col_name,
    hla_allele_col_name = hla_allele_col_name,
    read_name_col_name = read_name_col_name,
    read_seq_col_name = read_seq_col_name,
    ...
  )

  if (strand_col_name %in% names(reads)) {
    strand <- as.character(reads[[strand_col_name]])
    strand[is.na(strand)] <- "<NA>"
    for (strand_value in unique(strand)) {
      strand_rows <- strand == strand_value
      strand_res <- single_res[strand_rows, , drop = FALSE]
      strand_matched <- sum(Matrix::rowSums(strand_res) > 0)
      strand_unmatched <- nrow(strand_res) - strand_matched
      message("Strand (", strand_value, "):")
      message("  ", strand_matched, " reads with at least one match/hit (",
              round(strand_matched / nrow(strand_res) * 100), " %)")
      message("  ", strand_unmatched, " reads with no match/hit (",
              round(strand_unmatched / nrow(strand_res) * 100), " %)")
    }
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
  total_reads <- reads_w_min_one_match_sum + reads_w_no_match_sum
  message("  ", reads_w_min_one_match_sum,
          " reads with at least one match/hit (",
          round(reads_w_min_one_match_sum / total_reads * 100), " %)")
  message("  ", reads_w_no_match_sum, " reads with no match/hit (",
          round(reads_w_no_match_sum / total_reads * 100), " %)")

  retained_alleles <- which(
    expl_reads_per_allele >= max(expl_reads_per_allele) / allele_diff
  )
  if (length(retained_alleles) < 2L) {
    stop(
      "Fewer than two alleles passed the allele_diff filter; ",
      "pairwise matching cannot be performed."
    )
  }

  # A dense matrix speeds up the compiled pairwise operation. Preserve matrix
  # dimensions when only one read matched.
  top_single_res <- as.matrix(
    single_res[reads_w_min_one_match, retained_alleles, drop = FALSE]
  )
  top_sin_res_df <-
    data.frame(expl_reads = Matrix::colSums(top_single_res)) |>
    tibble::rownames_to_column(hla_allele_col_name) |>
    dplyr::left_join(hla_ref, by = hla_allele_col_name) |>
    dplyr::mutate(rank = dplyr::dense_rank(-expl_reads))
  top_sin_res_df$allele_group <- stringr::str_extract(
    top_sin_res_df[[hla_allele_col_name]],
    "[[:alpha:]]+\\*[[:digit:]]{2}"
  )

  col.combs <- t(utils::combn(seq_len(ncol(top_single_res)), 2L))
  message("Calculating pairwise matches. Combinations: ", nrow(col.combs), ".")

  # doing this in R was too slow. other packages did not have the functionality
  # so, written in c++
  if (identical(lapply_fun, parallel::mclapply) && "mc.cores" %in% names(arg_list)) {
    # Split the matrix into chunks for multithreading
    pairwise_results <- lapply_fun(brathering::split_mat(col.combs, n_chunks = arg_list[["mc.cores"]], byrow = TRUE), function(x) {
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
    dplyr::mutate(allele_comb = igsc:::orderAndConcatenateStrings(as.matrix(pair_res_df[c("allele1", "allele2")])), .after = "allele2")|>
    #dplyr::mutate(rank_int = dplyr::dense_rank(base::interaction(-tot_expl_reads, -uni_expl_reads, lex.order = TRUE)))|>
    dplyr::group_by(allele_comb)|>
    #dplyr::filter(rank_int == min(rank_int))|>
    dplyr::filter(rank_sum == min(rank_sum))|>
    dplyr::ungroup()|>
    dplyr::distinct(allele_comb, .keep_all = TRUE)|>
    dplyr::left_join(hla_ref|> dplyr::select(dplyr::all_of(c(hla_allele_col_name, p_group_col_name, g_group_col_name))), by = c("allele1" = hla_allele_col_name))|>
    dplyr::rename(
      p_group1 = dplyr::all_of(p_group_col_name),
      g_group1 = dplyr::all_of(g_group_col_name)
    ) |>
    dplyr::left_join(hla_ref|> dplyr::select(dplyr::all_of(c(hla_allele_col_name, p_group_col_name, g_group_col_name))), by = c("allele2" = hla_allele_col_name))|>
    dplyr::rename(
      p_group2 = dplyr::all_of(p_group_col_name),
      g_group2 = dplyr::all_of(g_group_col_name)
    ) |>
    dplyr::left_join(top_sin_res_df|> dplyr::select(dplyr::all_of(hla_allele_col_name), expl_reads), by = c("allele1" = hla_allele_col_name))|>
    dplyr::rename("allele1_expl_reads" = expl_reads)|>
    dplyr::left_join(top_sin_res_df|> dplyr::select(dplyr::all_of(hla_allele_col_name), expl_reads), by = c("allele2" = hla_allele_col_name))|>
    dplyr::rename("allele2_expl_reads" = expl_reads)|>
    dplyr::mutate(expl_reads_overall = !!reads_w_min_one_match_sum)|>
    dplyr::mutate(non_expl_reads_overall = !!reads_w_no_match_sum)|>
    dplyr::mutate(allele12_expl_read_diff = abs(allele1_expl_reads - allele2_expl_reads))|>
    dplyr::mutate(frac_match_reads_expl = tot_expl_reads/expl_reads_overall)|>
    dplyr::mutate(frac_all_reads_expl = tot_expl_reads/(expl_reads_overall + non_expl_reads_overall))|>
    dplyr::mutate(allele_group1 = stringr::str_extract(allele1, "[:alpha:]\\*[:digit:]{2}"))|>
    dplyr::mutate(allele_group2 = stringr::str_extract(allele2, "[:alpha:]\\*[:digit:]{2}"))

  # Treat pair labels as unordered so X/Y and Y/X are grouped together even
  # when the reference-table order interleaves their alleles.
  pair_res_df$p_group12 <- igsc:::orderAndConcatenateStrings(
    as.matrix(pair_res_df[c("p_group1", "p_group2")])
  )
  pair_res_df$g_group12 <- igsc:::orderAndConcatenateStrings(
    as.matrix(pair_res_df[c("g_group1", "g_group2")])
  )
  pair_res_df$allele_group12 <- igsc:::orderAndConcatenateStrings(
    as.matrix(pair_res_df[c("allele_group1", "allele_group2")])
  )

  ## plotting
  top_pair_res_plot1 <-
    #top_pair_res_df
    pair_res_df|>
    dplyr::group_by(p_group12)|>
    dplyr::slice_min(order_by = rank_sum, n = 1, with_ties = FALSE)|>
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
    ggplot2::geom_hline(yintercept = allele_group12_medians[1, 2, drop = TRUE], color = "tomato2") +
    ggplot2::geom_hline(yintercept = max(top_pair_res_plot1$frac_match_reads_expl), color = "forestgreen") +
    ggplot2::facet_wrap(ggplot2::vars(allele_group1), nrow = 1, scales = "free_x")

  top_pair_res_plot2 <-
    top_pair_res_plot1|>
    dplyr::arrange(rank_sum, dplyr::desc(tot_expl_reads),
                   dplyr::desc(uni_expl_reads)) |>
    dplyr::mutate(rank_plot = dplyr::row_number()) |>
    dplyr::slice_head(n = top_n_pairwise_results) |>
    dplyr::mutate(plot.color = as.factor(ifelse(rank_plot %% 2 != 0, 1, 2)))

  sin_plot <- ggplot2::ggplot(top_sin_res_df, ggplot2::aes(x = brathering::reorder_within(!!rlang::sym(hla_allele_col_name), expl_reads, allele_group), y = expl_reads)) +
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

  expl_reads_overall <- unique(top_pair_res_plot2$expl_reads_overall)
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
                            hla_seq_col_name,
                            hla_allele_col_name,
                            read_name_col_name,
                            read_seq_col_name,
                            ...) {

  read_sequences <- stats::setNames(
    reads[[read_seq_col_name]],
    reads[[read_name_col_name]]
  )

  # PDict requires equal-width patterns. Split first by read width and then
  # into bounded chunks to support variable-length read data without excessive
  # peak memory use.
  width_groups <- split(read_sequences, nchar(read_sequences))
  read_chunks <- unlist(
    lapply(
      width_groups,
      function(x) split(x, ceiling(seq_along(x) / 1000L))
    ),
    recursive = FALSE,
    use.names = FALSE
  )

  single_res <- lapply_fun(read_chunks, function(read_chunk) {
    hit_matrix <- Biostrings::vwhichPDict(
      subject = Biostrings::DNAStringSet(hla_ref[[hla_seq_col_name]]),
      pdict = Biostrings::PDict(read_chunk, max.mismatch = maxmis),
      max.mismatch = maxmis
    )
    hit_matrix <- lapply(
      hit_matrix,
      function(hit_indices) {
        replace(integer(length(read_chunk)), hit_indices, 1L)
      }
    )
    hit_matrix <- methods::as(do.call(cbind, hit_matrix), "sparseMatrix")
    colnames(hit_matrix) <- hla_ref[[hla_allele_col_name]]
    rownames(hit_matrix) <- names(read_chunk)
    hit_matrix
  }, ...)

  single_res <- do.call(rbind, single_res)
  single_res[match(names(read_sequences), rownames(single_res)), , drop = FALSE]
}



#' Simulate reads from two HLA-A alleles
#'
#' Generate fixed-length synthetic DNA reads by sampling substrings from two
#' HLA-A reference alleles. Reads are sampled according to the requested allele
#' proportions and can contain a fixed number of randomly positioned
#' single-nucleotide substitutions. The resulting data frame can be supplied
#' directly to [hla_typing()].
#'
#' Reference RNA sequences are converted from `U` to `T`. Reference rows with
#' missing sequences, characters other than `A`, `C`, `G`, and `T`, or
#' sequences shorter than `read_length` are excluded. When `allele_names` is
#' `NULL`, the first two suitable HLA-A alleles are used.
#'
#' SNP positions are sampled without replacement within each read. Each
#' selected nucleotide is replaced by one of the other three DNA bases, so the
#' number of substitutions in every read is exactly `snps_per_read`.
#'
#' Add decoy reads by running function twice and combine different allele_names.
#'
#' @param hla_ref A data frame containing HLA reference alleles and their
#'   sequences.
#' @param n A positive integer giving the number of reads to generate.
#' @param read_length A positive integer giving the read length in nucleotides.
#' @param allele_names A character vector of HLA allele
#'   names, or `NULL` to use random alleles in `hla_ref`.
#' @param proportions A numeric vector giving the relative
#'   sampling probabilities of `allele_names`. Values must be non-negative and
#'   have a positive sum; they are normalized internally.
#' @param snps_per_read A non-negative vector of integers giving allowed number
#'   of single-nucleotide substitutions inserted into reads. It cannot exceed
#'   `read_length`. The default of zero produces exact reference matches.
#' @param seq_col A character scalar naming the reference sequence column in
#'   `hla_ref`.
#' @param allele_col A character scalar naming the allele column in `hla_ref`.
#' @param seed An integer seed passed to [base::set.seed()] for reproducible
#'   sampling.
#' @param gene which HLA gene to sample from
#' @param minus_strand_prob probability of reads on minus strand. these will be
#'   reverse-complemented
#'
#' @return A data frame with `n` rows and the following columns:
#'   \describe{
#'     \item{readName}{A unique synthetic read identifier.}
#'     \item{seq}{The simulated DNA sequence, including inserted SNPs.}
#'     \item{strand}{The strand label `"+"`, included for compatibility with
#'       [hla_typing()].}
#'     \item{true_allele}{The reference allele from which the read was sampled.}
#'     \item{reference_start}{The one-based start position in the reference
#'       sequence.}
#'     \item{n_snps}{The number of substitutions inserted into the read.}
#'     \item{snp_positions}{A list-column containing the one-based positions of
#'       inserted substitutions within each read.}
#'   }
#'
#' @seealso [hla_typing()]
#' @export
#'
#' @examples
#' \dontrun{
#' reads <- simulate_hla_a_reads(
#'   hla_ref,
#'   n = 1000,
#'   read_length = 100,
#'   snps_per_read = 2
#' )
#'
#' table(reads$true_allele)
#' result <- hla_typing(hla_ref, reads, maxmis = 2)
#' }
simulate_hla_reads <- function(
    hla_ref,
    gene = c("A", "B", "C", "E", "F", "G", "DRA", "DRB1", "DRB3", "DRB4",
             "DRB5", "DPA1", "DPB1", "DQA1", "DQB1", "DMA", "DMB", "DOA",
             "DOB"),
    n = 1000L,
    read_length = 100L,
    allele_names = NULL,
    proportions = c(0.55, 0.45),
    snps_per_read = 0L,
    minus_strand_prob = 0,
    seq_col = "seq_Exon2_3",
    allele_col = "allele",
    seed = 42L
) {

  if (length(n) != 1L || is.na(n) || n < 1L || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (length(read_length) != 1L || is.na(read_length) ||
      read_length < 1L || read_length != as.integer(read_length)) {
    stop("read_length must be a positive integer.")
  }
  if (!is.numeric(snps_per_read) ||
      length(snps_per_read) == 0L ||
      anyNA(snps_per_read) ||
      any(!is.finite(snps_per_read)) ||
      any(snps_per_read < 0L) ||
      any(snps_per_read != floor(snps_per_read)) ||
      any(snps_per_read > read_length)) {
    stop(
      "snps_per_read must contain integers between zero and read_length."
    )
  }
  if (length(seq_col) != 1L || !seq_col %in% names(hla_ref)) {
    stop("seq_col must name a column in hla_ref.")
  }
  if (length(allele_col) != 1L || !allele_col %in% names(hla_ref)) {
    stop("allele_col must name a column in hla_ref.")
  }
  if (length(minus_strand_prob) != 1L ||
      is.na(minus_strand_prob) ||
      minus_strand_prob < 0 ||
      minus_strand_prob > 1) {
    stop("minus_strand_prob must be between zero and one.")
  }

  set.seed(seed)

  gene <- rlang::arg_match(gene)

  # get valid refs
  ref <- hla_ref[
    grepl(paste0("^(HLA-)?", gene, "\\*"), hla_ref[[allele_col]]) & !is.na(hla_ref[[seq_col]]),
    c(allele_col, seq_col),
    drop = FALSE]

  ref[[seq_col]] <- toupper(as.character(ref[[seq_col]]))
  ref[[seq_col]] <- chartr("U", "T", ref[[seq_col]])
  ref[[seq_col]] <- gsub("\\s+", "", ref[[seq_col]])

  ref <- ref[
    grepl("^[ACGT]+$", ref[[seq_col]]) & nchar(ref[[seq_col]]) >= read_length,
    ,
    drop = FALSE]
  ref <- ref[!duplicated(ref[[allele_col]]), , drop = FALSE]
  if (nrow(ref) < 2L) {
    stop(
      "At least two suitable HLA-", gene,
      " reference alleles are required."
    )
  }

  if (is.null(allele_names)) {
    allele_names <- sample(ref[[allele_col]], 2)
  } else {
    missing_alleles <- setdiff(allele_names, ref[[allele_col]])
    if (length(missing_alleles) > 0L) {
      stop(
        "Alleles not found or sequences too short: ",
        paste(missing_alleles, collapse = ", ")
      )
    }
  }
  if (anyDuplicated(allele_names)) {
    stop("allele_names must not contain duplicates.")
  }
  if (length(allele_names) != length(proportions)) {
    stop("allele_names must have same length as proportions.")
  }


  if (anyNA(proportions) ||
      any(!is.finite(proportions)) || any(proportions < 0) ||
      sum(proportions) == 0) {
    stop(
      "proportions must contain one finite, non-negative value per allele ",
      "and have a positive sum."
    )
  }
  proportions <- proportions / sum(proportions)

  truth <- sample(
    allele_names,
    size = n,
    replace = TRUE,
    prob = proportions
  )

  reference_sequences <- stats::setNames(
    ref[[seq_col]][match(allele_names, ref[[allele_col]])],
    allele_names
  )

  snp_positions <- vector("list", n)
  alternatives <- list(
    A = c("C", "G", "T"),
    C = c("A", "G", "T"),
    G = c("A", "C", "T"),
    T = c("A", "C", "G")
  )

  strands <- sample(
    c("+", "-"),
    size = n,
    replace = TRUE,
    prob = c(1 - minus_strand_prob, minus_strand_prob)
  )

  reference <- reference_sequences[truth]
  starts <- purrr::map_int(nchar(reference) - read_length + 1L, ~sample.int(.x ,size = 1L))
  sequence <- stringr::str_sub(reference, starts, starts + read_length - 1L)
  sequence[which(strands == "-")] <- igsc::revcompDNA(sequence[which(strands == "-")])

  # insert SNPs
  # Insert SNPs into each individual read.
  for (i in seq_len(n)) {
    n_snps <- sample(snps_per_read, size = 1L)

    if (n_snps > 0L) {
      positions <- sort(
        sample.int(
          read_length,
          size = n_snps,
          replace = FALSE
        )
      )

      bases <- strsplit(sequence[i], "", fixed = TRUE)[[1]]

      for (position in positions) {
        bases[position] <- sample(
          alternatives[[bases[position]]],
          size = 1L
        )
      }

      sequence[i] <- paste0(bases, collapse = "")
      snp_positions[[i]] <- positions
    } else {
      snp_positions[[i]] <- integer()
    }
  }

  sequences <- sequence

  return(data.frame(
    readName = sprintf(paste0("synthetic_HLA_", gene, "_%04d"), seq_len(n)),
    seq = sequences,
    strand = strands,
    true_allele = truth,
    reference_start = starts,
    n_snps = lengths(snp_positions),
    snp_positions = I(snp_positions),
    stringsAsFactors = FALSE
  ))
}

