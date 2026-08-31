#' Summarize and visualize TCRen structure scores
#'
#' Join a TCRen feature table to a TCRen score table, rank candidate
#' TCR-pMHC structures, and construct four diagnostic plots covering docking
#' geometry, contact energetics, interface quality, and footprint topology.
#'
#' @param features Either a data frame containing TCRen features or a path to a
#'   delimited feature file readable by [vroom::vroom()]. The table must contain
#'   one row per `complex.id`.
#' @param scores Either a data frame containing TCRen scores or a path to a
#'   delimited score file readable by [vroom::vroom()]. The table must contain
#'   one row per `complex.id` and columns `Q` and `T`. `P_native` is optional.
#' @param top_n Positive integer giving the number of highest-ranked structures
#'   to label and return in `top_candidates`. Defaults to `10`.
#' @param rank_by Character vector specifying the preferred ranking-score order.
#'   The first column present in the joined data with at least one finite value
#'   is used. By default, `P_native` is preferred, followed by `Q` and `T`. This
#'   fallback is useful when TCRen skips `P_native` for a small cohort.
#' @param label_top Logical; if `TRUE`, label the top candidates on every plot.
#' @param point_alpha Numeric point opacity between zero and one.
#'
#' @details
#' `F_good` is defined as `-(F_tcr_pep + F_tcr_mhc)`, so larger values represent
#' more favorable TCRen contact scores. `footprint_contiguity` is
#' `-fp_b0_frac_r7`, so larger values represent a less fragmented footprint.
#'
#' `P_native` is cohort-relative and should only be compared among structures
#' generated for the same biological question and with a consistent modeling
#' protocol. The function does not treat a missing `P_native` column as an
#' error; it falls back according to `rank_by`.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{`data`}{A named list containing `features`, `scores`, `combined`,
#'     `ranking`, and `top_candidates` data frames.}
#'   \item{`plots`}{A named list containing the ggplot objects
#'     `geometry_vs_topology`, `geometry_vs_contact_energetics`,
#'     `burial_vs_hydrogen_bonds`, and
#'     `footprint_contiguity_vs_diversity`.}
#'   \item{`ranking_score`}{The score column selected from `rank_by`.}
#'   \item{`top_n`}{The effective number of top candidates returned.}
#' }
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' result <- plot_tcren_structures(
#'   features = "features.tsv",
#'   scores = "scores.tsv",
#'   top_n = 10
#' )
#'
#' result$data$ranking
#' result$plots$geometry_vs_topology
#' result$plots$geometry_vs_contact_energetics
#' }
plot_tcren_structures <- function(
    features,
    scores,
    top_n = 10L,
    rank_by = c("P_native", "Q", "T"),
    label_top = TRUE,
    point_alpha = 0.85) {

  read_input <- function(x, object_name) {
    if (is.character(x) && length(x) == 1L && !is.na(x)) {
      if (!file.exists(x)) {
        stop(object_name, " file does not exist: ", x, call. = FALSE)
      }
      return(vroom::vroom(x, show_col_types = FALSE, progress = FALSE))
    }

    if (inherits(x, "data.frame")) {
      return(tibble::as_tibble(x))
    }

    stop(
      object_name,
      " must be a data frame or a single file path.",
      call. = FALSE
    )
  }

  require_columns <- function(x, required, object_name) {
    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
      stop(
        object_name,
        " is missing required columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }

  check_unique_ids <- function(x, object_name) {
    duplicate_ids <- x |>
      dplyr::count(.data$complex.id, name = "n") |>
      dplyr::filter(.data$n > 1L)

    if (nrow(duplicate_ids) > 0L) {
      stop(
        object_name,
        " contains duplicated complex.id values: ",
        paste(utils::head(duplicate_ids$complex.id, 5L), collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("top_n must be one positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (!is.character(rank_by) || length(rank_by) < 1L) {
    stop("rank_by must be a non-empty character vector.", call. = FALSE)
  }

  if (length(label_top) != 1L || is.na(label_top)) {
    stop("label_top must be TRUE or FALSE.", call. = FALSE)
  }

  if (length(point_alpha) != 1L || is.na(point_alpha) ||
      point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be between zero and one.", call. = FALSE)
  }

  features_tbl <- read_input(features, "features")
  scores_tbl <- read_input(scores, "scores")

  require_columns(features_tbl, "complex.id", "features")
  require_columns(scores_tbl, c("complex.id", "Q", "T"), "scores")
  check_unique_ids(features_tbl, "features")
  check_unique_ids(scores_tbl, "scores")

  feature_columns <- c(
    "F_tcr_pep", "F_tcr_mhc", "n_clashes", "clash_score", "burial",
    "n_hbond", "n_pep_contacted", "chain_balance", "fp_b0_frac_r7",
    "D2_pep24", "H_cell", "n_pep_contacts"
  )
  require_columns(features_tbl, feature_columns, "features")

  score_columns <- c("Q", "T", "P_native")

  # Prefer score-table values if an earlier feature table already contains any
  # score columns, avoiding .x/.y suffixes after the join.
  dat <- features_tbl |>
    dplyr::select(-dplyr::any_of(score_columns)) |>
    dplyr::left_join(
      scores_tbl |>
        dplyr::select(.data$complex.id, dplyr::any_of(score_columns)),
      by = "complex.id"
    )

  unmatched <- dat |>
    dplyr::filter(is.na(.data$Q) & is.na(.data$T)) |>
    dplyr::pull(.data$complex.id)
  if (length(unmatched) > 0L) {
    stop(
      "No matching score row was found for complex.id: ",
      paste(utils::head(unmatched, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  usable_score <- function(column) {
    column %in% names(dat) &&
      is.numeric(dat[[column]]) &&
      any(is.finite(dat[[column]]))
  }

  selected_score <- rank_by[vapply(rank_by, usable_score, logical(1))][1]
  if (is.na(selected_score)) {
    stop(
      "None of the rank_by columns is present with finite numeric values: ",
      paste(rank_by, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- dat |>
    dplyr::mutate(
      F_good = -(.data$F_tcr_pep + .data$F_tcr_mhc),
      clash_free = .data$n_clashes == 0,
      footprint_contiguity = -.data$fp_b0_frac_r7,
      ranking_score = .data[[selected_score]],
      candidate_rank = dplyr::min_rank(dplyr::desc(.data$ranking_score)),
      score_percentile = dplyr::percent_rank(.data$ranking_score)
    )

  ranking <- dat |>
    dplyr::arrange(
      dplyr::desc(.data$ranking_score),
      dplyr::desc(.data$Q),
      dplyr::desc(.data$T),
      .data$n_clashes,
      .data$clash_score
    ) |>
    dplyr::select(
      .data$candidate_rank,
      .data$complex.id,
      ranking_score = .data$ranking_score,
      dplyr::any_of(c("P_native", "Q", "T")),
      .data$burial,
      .data$n_hbond,
      .data$n_pep_contacted,
      .data$chain_balance,
      .data$F_good,
      .data$n_clashes,
      .data$clash_score,
      .data$score_percentile
    )

  effective_top_n <- min(top_n, nrow(dat))
  top_candidates <- dat |>
    dplyr::slice_min(
      order_by = .data$candidate_rank,
      n = effective_top_n,
      with_ties = FALSE
    )

  clash_shape_scale <- ggplot2::scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 4),
    labels = c(`TRUE` = "No clashes", `FALSE` = "Clashes")
  )

  geometry_vs_topology <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data$Q, y = .data$T)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data[[selected_score]],
        size = .data$burial,
        shape = .data$clash_free
      ),
      alpha = point_alpha
    ) +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    clash_shape_scale +
    ggplot2::labs(
      x = "Q: interface geometry",
      y = "T: footprint topology",
      colour = selected_score,
      size = "Buried area",
      shape = NULL
    ) +
    ggplot2::theme_classic()

  geometry_vs_contact_energetics <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data$Q, y = .data$F_good)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data[[selected_score]],
        shape = .data$clash_free,
        size = .data$n_pep_contacted
      ),
      alpha = point_alpha
    ) +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    clash_shape_scale +
    ggplot2::labs(
      x = "Q: interface geometry",
      y = "-(F_tcr_pep + F_tcr_mhc): favorable contact score",
      colour = selected_score,
      size = "Peptide residues contacted",
      shape = NULL
    ) +
    ggplot2::theme_classic()

  burial_vs_hydrogen_bonds <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data$burial, y = .data$n_hbond)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data$chain_balance,
        size = .data$n_pep_contacted,
        shape = .data$clash_free
      ),
      alpha = point_alpha
    ) +
    ggplot2::scale_colour_viridis_c() +
    clash_shape_scale +
    ggplot2::labs(
      x = "Buried interface area (Å²)",
      y = "TCR-peptide hydrogen bonds",
      colour = "Chain balance",
      size = "Peptide residues contacted",
      shape = NULL
    ) +
    ggplot2::theme_classic()

  footprint_contiguity_vs_diversity <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$footprint_contiguity,
      y = .data$D2_pep24
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data$H_cell,
        size = .data$n_pep_contacts
      ),
      alpha = point_alpha
    ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(
      x = "-fp_b0_frac_r7: higher means more contiguous",
      y = "D2_pep24: footprint diversity",
      colour = "Footprint entropy",
      size = "Peptide contacts"
    ) +
    ggplot2::theme_classic()

  if (isTRUE(label_top)) {
    geometry_vs_topology <- geometry_vs_topology +
      ggrepel::geom_text_repel(
        data = top_candidates,
        ggplot2::aes(x = .data$Q, y = .data$T, label = .data$complex.id),
        inherit.aes = FALSE,
        size = 3
      )

    geometry_vs_contact_energetics <- geometry_vs_contact_energetics +
      ggrepel::geom_text_repel(
        data = top_candidates,
        ggplot2::aes(
          x = .data$Q,
          y = .data$F_good,
          label = .data$complex.id
        ),
        inherit.aes = FALSE,
        size = 3
      )

    burial_vs_hydrogen_bonds <- burial_vs_hydrogen_bonds +
      ggrepel::geom_text_repel(
        data = top_candidates,
        ggplot2::aes(
          x = .data$burial,
          y = .data$n_hbond,
          label = .data$complex.id
        ),
        inherit.aes = FALSE,
        size = 3
      )

    footprint_contiguity_vs_diversity <-
      footprint_contiguity_vs_diversity +
      ggrepel::geom_text_repel(
        data = top_candidates,
        ggplot2::aes(
          x = .data$footprint_contiguity,
          y = .data$D2_pep24,
          label = .data$complex.id
        ),
        inherit.aes = FALSE,
        size = 3
      )
  }

  list(
    data = list(
      features = features_tbl,
      scores = scores_tbl,
      combined = dat,
      ranking = ranking,
      top_candidates = top_candidates
    ),
    plots = list(
      geometry_vs_topology = geometry_vs_topology,
      geometry_vs_contact_energetics = geometry_vs_contact_energetics,
      burial_vs_hydrogen_bonds = burial_vs_hydrogen_bonds,
      footprint_contiguity_vs_diversity =
        footprint_contiguity_vs_diversity
    ),
    ranking_score = selected_score,
    top_n = effective_top_n
  )
}


#' Triage and rank TCRen candidate complexes
#'
#' Apply conservative structural-quality filters to a joined TCRen feature and
#' score table, record specific red flags and cautions for every complex, and
#' rank the candidates that survive. The function deliberately separates clear
#' exclusion criteria from softer evidence that warrants inspection.
#'
#' @param dat A data frame containing joined TCRen features and scores, with one
#'   row per `complex.id`.
#' @param thresholds Named list overriding any of the conservative defaults:
#'   `severe_n_clashes = 10`, `severe_clash_score = 5`,
#'   `min_n_contacts_tp = 5`, `min_n_pep_contacted = 3`,
#'   `min_burial = 500`, `max_burial = 3000`, `low_Q = -2`,
#'   `low_T = 0.10`, `low_P_native = 0.10`, `low_p_real = 0.10`,
#'   `low_chain_balance = 0.05`, and
#'   `fragmented_fp_b0_frac_r7 = 0.50`. Distances and burial use the units in
#'   the TCRen output. These defaults are conservative triage heuristics, not
#'   experimentally validated binding cutoffs.
#' @param rank_weights Named numeric vector giving relative weights for
#'   within-cohort percentile ranks. Defaults to `Q = 0.35`,
#'   `P_native = 0.25`, `T = 0.15`, `p_real = 0.15`, and `F_good = 0.10`.
#'   Missing score columns are omitted and the remaining weights are
#'   renormalized separately for every row.
#' @param caution_penalty Non-negative penalty subtracted from the percentile
#'   promise score for each caution. Defaults to `0.05`.
#' @param drop_excluded Logical. If `FALSE` (the default), return every complex
#'   so exclusions remain auditable. If `TRUE`, return only rows for which
#'   `keep_candidate` is `TRUE`.
#'
#' @details
#' Clear red flags are restricted to missing or invalid critical values, severe
#' steric clashes, a simultaneously sparse TCR-peptide contact interface, or
#' very low buried interface area. Low `Q`, low topology or recognition scores,
#' mild clashes, chain imbalance, absent hydrogen bonds, unusual burial, and a
#' fragmented footprint are recorded as cautions rather than automatic
#' exclusions.
#'
#' Ranking is performed only among retained candidates. `promise_score` is the
#' weighted mean of available within-cohort percentiles, and
#' `adjusted_promise_score` subtracts `caution_penalty` for each caution. This is
#' a transparent prioritization device, not a binding probability. In
#' particular, `P_native` and `T` remain cohort-relative.
#'
#' @return A tibble containing all original columns plus derived energies,
#'   individual flag columns, `red_flags`, `cautions`, `keep_candidate`,
#'   `promise_score`, `adjusted_promise_score`, `candidate_rank`, and the verbose
#'   `selection_rationale`. Excluded structures have `candidate_rank = NA`.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' triage <- triage_tcren_candidates(joined_features_and_scores)
#'
#' triage |>
#'   dplyr::select(
#'     complex.id, keep_candidate, candidate_rank,
#'     red_flags, cautions, selection_rationale
#'   )
#'
#' shortlist <- triage_tcren_candidates(
#'   joined_features_and_scores,
#'   drop_excluded = TRUE
#' )
#' }
triage_tcren_candidates <- function(
    dat,
    thresholds = list(),
    rank_weights = c(
      Q = 0.35,
      P_native = 0.25,
      T = 0.15,
      p_real = 0.15,
      F_good = 0.10
    ),
    caution_penalty = 0.05,
    drop_excluded = FALSE) {

  if (!inherits(dat, "data.frame")) {
    stop("dat must be a data frame.", call. = FALSE)
  }

  default_thresholds <- list(
    severe_n_clashes = 10,
    severe_clash_score = 5,
    min_n_contacts_tp = 5,
    min_n_pep_contacted = 3,
    min_burial = 500,
    max_burial = 3000,
    low_Q = -2,
    low_T = 0.10,
    low_P_native = 0.10,
    low_p_real = 0.10,
    low_chain_balance = 0.05,
    fragmented_fp_b0_frac_r7 = 0.50
  )

  if (!is.list(thresholds) || is.null(names(thresholds))) {
    if (length(thresholds) > 0L) {
      stop("thresholds must be a named list.", call. = FALSE)
    }
  }

  unknown_thresholds <- setdiff(names(thresholds), names(default_thresholds))
  if (length(unknown_thresholds) > 0L) {
    stop(
      "Unknown threshold names: ",
      paste(unknown_thresholds, collapse = ", "),
      call. = FALSE
    )
  }
  cut <- utils::modifyList(default_thresholds, thresholds)

  cut_values <- unlist(cut, use.names = TRUE)
  if (!all(is.finite(cut_values))) {
    stop("Every threshold must be a finite number.", call. = FALSE)
  }

  if (!is.numeric(rank_weights) || length(rank_weights) < 1L ||
      is.null(names(rank_weights)) || any(names(rank_weights) == "") ||
      any(!is.finite(rank_weights)) || any(rank_weights < 0) ||
      sum(rank_weights) <= 0) {
    stop(
      "rank_weights must be a named, non-negative numeric vector with positive sum.",
      call. = FALSE
    )
  }

  if (length(caution_penalty) != 1L || !is.finite(caution_penalty) ||
      caution_penalty < 0) {
    stop("caution_penalty must be one non-negative number.", call. = FALSE)
  }

  if (length(drop_excluded) != 1L || is.na(drop_excluded)) {
    stop("drop_excluded must be TRUE or FALSE.", call. = FALSE)
  }

  required <- c(
    "complex.id", "Q", "n_clashes", "clash_score", "burial",
    "n_contacts_tp", "n_pep_contacted"
  )
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "dat is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  duplicate_ids <- dat |>
    dplyr::count(.data$complex.id, name = "n") |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_ids) > 0L) {
    stop(
      "dat contains duplicated complex.id values: ",
      paste(utils::head(duplicate_ids$complex.id, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(dat)

  if (!"F_good" %in% names(out)) {
    if (all(c("F_tcr_pep", "F_tcr_mhc") %in% names(out))) {
      out <- out |>
        dplyr::mutate(F_good = -(.data$F_tcr_pep + .data$F_tcr_mhc))
    } else {
      out$F_good <- NA_real_
    }
  }

  optional_numeric <- c(
    "T", "P_native", "p_real", "p_real_bn", "F_good", "chain_balance",
    "n_hbond", "fp_b0_frac_r7"
  )
  for (column in setdiff(optional_numeric, names(out))) {
    out[[column]] <- NA_real_
  }

  percentile_rank <- function(x) {
    x <- as.numeric(x)
    answer <- rep(NA_real_, length(x))
    ok <- is.finite(x)
    n_ok <- sum(ok)

    if (n_ok == 1L) {
      answer[ok] <- 0.5
    } else if (n_ok > 1L) {
      answer[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n_ok - 1)
    }
    answer
  }

  evidence_columns <- intersect(names(rank_weights), names(out))
  evidence_columns <- evidence_columns[
    vapply(
      evidence_columns,
      function(column) any(is.finite(out[[column]])),
      logical(1)
    )
  ]
  if (length(evidence_columns) < 1L) {
    stop(
      "None of the rank_weights columns contains finite values.",
      call. = FALSE
    )
  }

  percentile_names <- paste0(evidence_columns, "_percentile")
  for (i in seq_along(evidence_columns)) {
    out[[percentile_names[[i]]]] <- percentile_rank(
      out[[evidence_columns[[i]]]]
    )
  }

  score_matrix <- as.matrix(out[, percentile_names, drop = FALSE])
  storage.mode(score_matrix) <- "double"
  weights <- unname(rank_weights[evidence_columns])
  weight_matrix <- matrix(
    weights,
    nrow = nrow(score_matrix),
    ncol = length(weights),
    byrow = TRUE
  )
  observed <- is.finite(score_matrix)
  weighted_sum <- rowSums(
    ifelse(observed, score_matrix * weight_matrix, 0),
    na.rm = TRUE
  )
  observed_weight <- rowSums(ifelse(observed, weight_matrix, 0))
  out$promise_score <- ifelse(
    observed_weight > 0,
    weighted_sum / observed_weight,
    NA_real_
  )

  critical_matrix <- as.data.frame(out[, required[-1], drop = FALSE])
  out$flag_missing_critical <- apply(
    critical_matrix,
    1L,
    function(x) any(!is.finite(as.numeric(x)))
  )

  out$flag_invalid_values <-
    (!is.na(out$n_clashes) & out$n_clashes < 0) |
    (!is.na(out$clash_score) & out$clash_score < 0) |
    (!is.na(out$burial) & out$burial < 0) |
    (!is.na(out$n_contacts_tp) & out$n_contacts_tp < 0) |
    (!is.na(out$n_pep_contacted) & out$n_pep_contacted < 0) |
    (is.finite(out$chain_balance) &
       (out$chain_balance < 0 | out$chain_balance > 0.5)) |
    (is.finite(out$T) & (out$T < 0 | out$T > 1)) |
    (is.finite(out$P_native) &
       (out$P_native < 0 | out$P_native > 1)) |
    (is.finite(out$p_real) & (out$p_real < 0 | out$p_real > 1))

  out$flag_severe_clashes <-
    (is.finite(out$n_clashes) &
       out$n_clashes >= cut$severe_n_clashes) |
    (is.finite(out$clash_score) &
       out$clash_score >= cut$severe_clash_score)

  out$flag_sparse_peptide_interface <-
    is.finite(out$n_contacts_tp) &
    is.finite(out$n_pep_contacted) &
    out$n_contacts_tp < cut$min_n_contacts_tp &
    out$n_pep_contacted < cut$min_n_pep_contacted

  out$flag_very_low_burial <-
    is.finite(out$burial) & out$burial < cut$min_burial

  out$caution_mild_clashes <-
    is.finite(out$n_clashes) & out$n_clashes > 0 &
    !out$flag_severe_clashes

  out$caution_low_contact_count <-
    !out$flag_sparse_peptide_interface &
    ((is.finite(out$n_contacts_tp) &
        out$n_contacts_tp < cut$min_n_contacts_tp) |
       (is.finite(out$n_pep_contacted) &
          out$n_pep_contacted < cut$min_n_pep_contacted))

  out$caution_low_Q <- is.finite(out$Q) & out$Q < cut$low_Q
  out$caution_low_T <- is.finite(out$T) & out$T < cut$low_T
  out$caution_low_P_native <-
    is.finite(out$P_native) & out$P_native < cut$low_P_native
  out$caution_low_p_real <-
    is.finite(out$p_real) & out$p_real < cut$low_p_real
  out$caution_chain_imbalance <-
    is.finite(out$chain_balance) &
    out$chain_balance < cut$low_chain_balance
  out$caution_no_hydrogen_bonds <-
    is.finite(out$n_hbond) & out$n_hbond == 0
  out$caution_high_burial <-
    is.finite(out$burial) & out$burial > cut$max_burial
  out$caution_fragmented_footprint <-
    is.finite(out$fp_b0_frac_r7) &
    out$fp_b0_frac_r7 > cut$fragmented_fp_b0_frac_r7

  q_pct <- percentile_rank(out$Q)
  f_pct <- percentile_rank(out$F_good)
  out$caution_energy_geometry_discordance <-
    is.finite(q_pct) & is.finite(f_pct) & q_pct <= 0.25 & f_pct >= 0.75

  red_flag_columns <- c(
    "flag_missing_critical",
    "flag_invalid_values",
    "flag_severe_clashes",
    "flag_sparse_peptide_interface",
    "flag_very_low_burial"
  )
  caution_columns <- c(
    "caution_mild_clashes",
    "caution_low_contact_count",
    "caution_low_Q",
    "caution_low_T",
    "caution_low_P_native",
    "caution_low_p_real",
    "caution_chain_imbalance",
    "caution_no_hydrogen_bonds",
    "caution_high_burial",
    "caution_fragmented_footprint",
    "caution_energy_geometry_discordance"
  )

  out$red_flag_count <- rowSums(out[, red_flag_columns, drop = FALSE])
  out$caution_count <- rowSums(out[, caution_columns, drop = FALSE])
  out$keep_candidate <- out$red_flag_count == 0L
  out$adjusted_promise_score <-
    out$promise_score - caution_penalty * out$caution_count

  fmt <- function(x, digits = 2L) {
    if (!is.finite(x)) {
      return("NA")
    }
    formatC(x, digits = digits, format = "f")
  }

  row_red_flags <- function(i) {
    messages <- character()

    if (out$flag_missing_critical[[i]]) {
      missing_here <- required[-1][
        !vapply(
          out[i, required[-1], drop = FALSE],
          function(x) is.finite(as.numeric(x[[1]])),
          logical(1)
        )
      ]
      messages <- c(
        messages,
        paste0(
          "missing critical metric(s): ",
          paste(missing_here, collapse = ", ")
        )
      )
    }
    if (out$flag_invalid_values[[i]]) {
      messages <- c(messages, "invalid or out-of-range metric value")
    }
    if (out$flag_severe_clashes[[i]]) {
      messages <- c(
        messages,
        paste0(
          "severe steric clashes (n_clashes=", fmt(out$n_clashes[[i]], 0L),
          ", clash_score=", fmt(out$clash_score[[i]]), ")"
        )
      )
    }
    if (out$flag_sparse_peptide_interface[[i]]) {
      messages <- c(
        messages,
        paste0(
          "very sparse TCR-peptide interface (", fmt(out$n_contacts_tp[[i]], 0L),
          " contact pairs across ", fmt(out$n_pep_contacted[[i]], 0L),
          " peptide residues)"
        )
      )
    }
    if (out$flag_very_low_burial[[i]]) {
      messages <- c(
        messages,
        paste0("very low buried area (", fmt(out$burial[[i]], 0L), " Å²)")
      )
    }

    if (length(messages) == 0L) "None" else paste(messages, collapse = "; ")
  }

  row_cautions <- function(i) {
    messages <- character()

    if (out$caution_mild_clashes[[i]]) {
      messages <- c(
        messages,
        paste0(
          "non-severe clashes (n=", fmt(out$n_clashes[[i]], 0L),
          ", score=", fmt(out$clash_score[[i]]), ")"
        )
      )
    }
    if (out$caution_low_contact_count[[i]]) {
      messages <- c(
        messages,
        paste0(
          "limited peptide engagement (", fmt(out$n_contacts_tp[[i]], 0L),
          " pairs; ", fmt(out$n_pep_contacted[[i]], 0L), " residues)"
        )
      )
    }
    if (out$caution_low_Q[[i]]) {
      messages <- c(messages, paste0("low Q=", fmt(out$Q[[i]])))
    }
    if (out$caution_low_T[[i]]) {
      messages <- c(messages, paste0("low T=", fmt(out$T[[i]])))
    }
    if (out$caution_low_P_native[[i]]) {
      messages <- c(
        messages,
        paste0("low P_native=", fmt(out$P_native[[i]]))
      )
    }
    if (out$caution_low_p_real[[i]]) {
      messages <- c(messages, paste0("low p_real=", fmt(out$p_real[[i]])))
    }
    if (out$caution_chain_imbalance[[i]]) {
      messages <- c(
        messages,
        paste0("strong alpha/beta imbalance (", fmt(out$chain_balance[[i]]), ")")
      )
    }
    if (out$caution_no_hydrogen_bonds[[i]]) {
      messages <- c(messages, "no detected TCR-peptide hydrogen bonds")
    }
    if (out$caution_high_burial[[i]]) {
      messages <- c(
        messages,
        paste0("unusually high burial (", fmt(out$burial[[i]], 0L), " Å²)")
      )
    }
    if (out$caution_fragmented_footprint[[i]]) {
      messages <- c(
        messages,
        paste0(
          "fragmented footprint (fp_b0_frac_r7=",
          fmt(out$fp_b0_frac_r7[[i]]), ")"
        )
      )
    }
    if (out$caution_energy_geometry_discordance[[i]]) {
      messages <- c(
        messages,
        "favorable contact energy but bottom-quartile geometry"
      )
    }

    if (length(messages) == 0L) "None" else paste(messages, collapse = "; ")
  }

  out$red_flags <- vapply(seq_len(nrow(out)), row_red_flags, character(1))
  out$cautions <- vapply(seq_len(nrow(out)), row_cautions, character(1))

  out$candidate_rank <- NA_integer_
  retained <- which(out$keep_candidate & is.finite(out$adjusted_promise_score))
  if (length(retained) > 0L) {
    retained_order <- retained[
      order(
        -out$adjusted_promise_score[retained],
        out$caution_count[retained],
        -out$Q[retained],
        out$clash_score[retained],
        na.last = TRUE
      )
    ]
    out$candidate_rank[retained_order] <- seq_along(retained_order)
  }

  n_retained <- sum(out$keep_candidate)
  row_description <- function(i) {
    metrics <- paste0(
      "Q=", fmt(out$Q[[i]]),
      ", T=", fmt(out$T[[i]]),
      ", P_native=", fmt(out$P_native[[i]]),
      ", p_real=", fmt(out$p_real[[i]]),
      ", F_good=", fmt(out$F_good[[i]]),
      ", peptide contacts=", fmt(out$n_contacts_tp[[i]], 0L),
      "/", fmt(out$n_pep_contacted[[i]], 0L),
      ", clashes=", fmt(out$n_clashes[[i]], 0L),
      " (score ", fmt(out$clash_score[[i]]), ")"
    )

    if (!out$keep_candidate[[i]]) {
      return(paste0(
        "Exclude from the shortlist because: ", out$red_flags[[i]],
        ". Favorable aggregate scores do not override this structural/QC failure. ",
        "Observed metrics: ", metrics, "."
      ))
    }

    strengths <- character()
    if (is.finite(out$Q[[i]]) && out$Q[[i]] >= 0) {
      strengths <- c(strengths, "geometry at or above the native-reference mean")
    }
    if (is.finite(out$T[[i]]) && out$T[[i]] >= 0.5) {
      strengths <- c(strengths, "supportive footprint topology")
    }
    if (is.finite(out$P_native[[i]]) && out$P_native[[i]] >= 0.5) {
      strengths <- c(strengths, "high cohort-relative P_native")
    }
    if (is.finite(out$p_real[[i]]) && out$p_real[[i]] >= 0.5) {
      strengths <- c(strengths, "real-like frozen-recognizer score")
    }
    if (is.finite(f_pct[[i]]) && f_pct[[i]] >= 0.75) {
      strengths <- c(strengths, "top-quartile contact energetics")
    }
    if (is.finite(out$n_clashes[[i]]) && out$n_clashes[[i]] == 0) {
      strengths <- c(strengths, "no detected clashes")
    }
    if (length(strengths) == 0L) {
      strengths <- "passes all conservative exclusion filters"
    }

    caution_text <- if (out$cautions[[i]] == "None") {
      "No configured cautions."
    } else {
      paste0("Cautions: ", out$cautions[[i]], ".")
    }

    paste0(
      if (out$candidate_rank[[i]] == 1L) "Favor" else "Retain",
      " as rank ", out$candidate_rank[[i]], " of ", n_retained,
      " retained candidates (adjusted promise score=",
      fmt(out$adjusted_promise_score[[i]], 3L),
      "; unpenalized score=", fmt(out$promise_score[[i]], 3L), "). ",
      "Strengths: ", paste(strengths, collapse = "; "), ". ",
      caution_text, " Observed metrics: ", metrics, "."
    )
  }

  out$selection_rationale <- vapply(
    seq_len(nrow(out)),
    row_description,
    character(1)
  )

  out <- out |>
    dplyr::arrange(
      dplyr::desc(.data$keep_candidate),
      .data$candidate_rank,
      dplyr::desc(.data$adjusted_promise_score),
      .data$complex.id
    )

  if (isTRUE(drop_excluded)) {
    out <- out |>
      dplyr::filter(.data$keep_candidate)
  }

  out
}
