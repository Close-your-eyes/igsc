#' Assign Memorable Names to T-Cell Receptor Clonotypes
#'
#' @description
#' Converts a long-format T-cell receptor (TCR) annotation table to paired
#' alpha/beta-chain representations, identifies clonotypes from shared CDR3
#' sequences, and assigns a random human-readable name to each clonotype. The
#' result contains both long- and wide-format tables.
#'
#' When `single_chains_only = TRUE`, the function keeps the annotation with the
#' highest UMI count and then the highest read count for each barcode and chain.
#' Remaining ties are resolved deterministically by CDR3 sequence and returned
#' separately for inspection. When `strict_TCR_biology = TRUE`, cells with more
#' than one beta chain or more than two alpha chains are excluded from naming
#' and later restored with `NA` as their clonotype name. Cells without any
#' non-missing CDR3 sequence also receive `NA`.
#'
#' @param cl_long A data frame in long format with one row per TCR-chain
#'   annotation. It must contain the columns identified by `sample_col`,
#'   `barcode_col`, `chain_col`, `cdr3_col`, and `clonotype_col`. The chain
#'   column must contain only `"TRA"` and `"TRB"`. If
#'   `single_chains_only = TRUE`, columns `umis` and `reads` are also required.
#' @param single_chains_only A single logical value. If `TRUE`, retain at most
#'   one annotation per barcode and chain, prioritizing UMI count and then read
#'   count.
#' @param strict_TCR_biology A single logical value. If `TRUE`, do not assign
#'   clonotype names to cells with more than one TRB chain or more than two TRA
#'   chains.
#' @param same_cl_name_across_donor_only_when_full_annotation A single logical
#'   value. If `TRUE`, a clonotype name may be shared across donors only when
#'   both its TRA and TRB CDR3 sequences are annotated.
#' @param shared_cdr3_mode How shared CDR3 sequences are interpreted within a
#'   sample. `"compatible"` (the default) joins cells that share a TRA and/or
#'   TRB CDR3 only when their other observed chain does not conflict. A
#'   single-chain observation that is compatible with multiple conflicting
#'   paired receptors is left as its own clonotype. `"any"` uses ordinary
#'   graph connectivity and joins cells whenever either chain is shared.
#' @param names_to_avoid A character vector of names that must not be assigned,
#'   or `NULL` to impose no additional exclusions.
#' @param pick_randomNames_args A named list of additional arguments passed to
#'   [igsc::pick_randomNames()] through [Gmisc::fastDoCall()]. The function
#'   supplies `n` and `names_to_avoid`; values supplied here should not conflict
#'   with those arguments.
#' @param sample_fields A character vector naming the fields encoded in the
#'   sample identifier, in order. It must include `"donor"`, and its length
#'   must match the number of fields in every sample identifier.
#' @param sample_field_sep A single character string separating fields in the
#'   sample identifier.
#' @param sample_col A single character string naming the sample column.
#' @param barcode_col A single character string naming the cell-barcode column.
#' @param chain_col A single character string naming the TCR-chain column.
#' @param cdr3_col A single character string naming the amino-acid CDR3 column.
#' @param clonotype_col A single character string naming the source clonotype-ID
#'   column.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{`cl_long`}{The input annotations with the assigned `cl_name`
#'     column joined in.}
#'   \item{`cl_wide`}{One row per sample and barcode, with chain-specific
#'     annotation columns suffixed by `_TRA` or `_TRB`, the parsed
#'     `sample_fields`, and `cl_name`.}
#'   \item{`multi_anno_TCR_equal_umi_read`}{Rows involved in unresolved
#'     multi-annotation ties after filtering by maximum UMI and read counts, or
#'     `NULL` when no such ties are present or `single_chains_only = FALSE`.}
#' }
#'
#' @details
#' Matching is performed once with connected components rather than by a fixed
#' number of alternating TRA/TRB collapse passes. Within samples, compatible
#' shared TRA or TRB sequences provide links. Exact receptor signatures are
#' then propagated across samples from the same donor. Across donors, exact
#' signatures are propagated only for paired TRA+TRB receptors when
#' `same_cl_name_across_donor_only_when_full_annotation = TRUE`.
#'
#' In `"compatible"` mode, every merged component must be explainable by one
#' maximal observed sequence set for each chain. This permits dropout and a
#' directly observed dual-alpha receptor while preventing an unpaired TRA- or
#' TRB-only record from bridging two contradictory paired receptors. Random
#' names are generated with [igsc::pick_randomNames()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' named_clonotypes <- add_clonotype_names(tcr_annotations)
#'
#' named_clonotypes$cl_long
#' named_clonotypes$cl_wide
#' }
add_clonotype_names_chatty <- function(cl_long,
                                single_chains_only = T,
                                strict_TCR_biology = F,
                                same_cl_name_across_donor_only_when_full_annotation = T,
                                shared_cdr3_mode = c("compatible", "any"),
                                names_to_avoid = NULL,
                                pick_randomNames_args = list(max_iter = 10000,
                                                             randomNames_args = list(which.names = "first")),
                                sample_fields = c("virus", "donor", "genotype", "treatment"),
                                sample_field_sep = "_",
                                sample_col = "sample",
                                barcode_col = "barcode",
                                chain_col = "chain",
                                cdr3_col = "cdr3",
                                clonotype_col = "clonotype_id") {

  shared_cdr3_mode <- match.arg(shared_cdr3_mode)

  logical_args <- list(
    single_chains_only = single_chains_only,
    strict_TCR_biology = strict_TCR_biology,
    same_cl_name_across_donor_only_when_full_annotation =
      same_cl_name_across_donor_only_when_full_annotation
  )
  invalid_logical <- names(logical_args)[!vapply(
    logical_args,
    function(x) is.logical(x) && length(x) == 1L && !is.na(x),
    logical(1)
  )]
  if (length(invalid_logical)) {
    stop(
      "These arguments must each be one non-missing logical value: ",
      paste(invalid_logical, collapse = ", "),
      call. = FALSE
    )
  }

  required_cols <- c(
    sample_col,
    barcode_col,
    chain_col,
    cdr3_col,
    clonotype_col
  )

  missing_cols <- setdiff(required_cols, names(cl_long))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(anyNA(cl_long[[sample_col]]),
          anyNA(cl_long[[barcode_col]]),
          anyNA(cl_long[[chain_col]]))) {
    stop("NA found in sample, barcode or chain.")
  }

  invalid_chains <- setdiff(unique(as.character(cl_long[[chain_col]])), c("TRA", "TRB"))
  if (length(invalid_chains)) {
    stop("chain_col must contain TRA and/or TRB only.")
  }

  if (!"donor" %in% sample_fields) {
    stop("sample_fields must have a 'donor' field. if not, then make a dummy.")
  }

  cl_long |>
    dplyr::distinct(!!rlang::sym(sample_col)) |>
    tidyr::separate_wider_delim(
      cols = dplyr::all_of(sample_col),
      delim = sample_field_sep,
      names = sample_fields,
      too_few = "error",
      too_many = "error")

  cdr3tra_col <- paste0(cdr3_col, "_TRA")
  cdr3trb_col <- paste0(cdr3_col, "_TRB")
  clonotypetra_col <- paste0(clonotype_col, "_TRA")
  clonotypetrb_col <- paste0(clonotype_col, "_TRB")

  multi_anno_TCR_equal_umi_read <- NULL
  if (single_chains_only) {
    if (any(!c("umis", "reads") %in% names(cl_long))) {
      stop("columns 'umis' and 'reads' needed in cl_long when single_chains_only.")
    }
    ranked_annotations <- cl_long |>
      dplyr::group_by(
        !!rlang::sym(sample_col),
        !!rlang::sym(barcode_col),
        !!rlang::sym(chain_col)
      ) |>
      dplyr::arrange(
        dplyr::desc(umis),
        dplyr::desc(reads),
        !!rlang::sym(cdr3_col),
        !!rlang::sym(clonotype_col),
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        .top_umis = dplyr::first(umis),
        .top_reads = dplyr::first(reads)
      )

    multi_anno_TCR_equal_umi_read <- ranked_annotations |>
      dplyr::filter(
        (is.na(umis) & is.na(.top_umis)) |
          (!is.na(umis) & !is.na(.top_umis) & umis == .top_umis),
        (is.na(reads) & is.na(.top_reads)) |
          (!is.na(reads) & !is.na(.top_reads) & reads == .top_reads)
      ) |>
      dplyr::filter(dplyr::n() > 1L) |>
      dplyr::ungroup() |>
      dplyr::select(-.top_umis, -.top_reads)

    if (!nrow(multi_anno_TCR_equal_umi_read)) {
      multi_anno_TCR_equal_umi_read <- NULL
    } else {
      message("Some sample/barcode/chain groups have annotations tied on UMI and read count.")
    }

    cl_long <- ranked_annotations |>
      dplyr::slice_head(n = 1L) |>
      dplyr::ungroup() |>
      dplyr::select(-.top_umis, -.top_reads)
  }

  cl_wide <- cl_long |>
    dplyr::select(
      !!rlang::sym(sample_col),
      !!rlang::sym(barcode_col),
      !!rlang::sym(chain_col),
      !!rlang::sym(cdr3_col),
      !!rlang::sym(clonotype_col)
    ) |>
    tidyr::pivot_wider(
      id_cols = c(dplyr::all_of(sample_col), dplyr::all_of(barcode_col)),
      names_from = !!rlang::sym(chain_col),
      values_from = c(!!rlang::sym(cdr3_col), !!rlang::sym(clonotype_col)),
      names_sep = "_",
      values_fn = collapse_order_fun2
    )

  for (missing_wide_col in c(cdr3tra_col, cdr3trb_col, clonotypetra_col, clonotypetrb_col)) {
    if (!missing_wide_col %in% names(cl_wide)) {
      cl_wide[[missing_wide_col]] <- NA_character_
    }
  }

  cl_wide <- cl_wide |>
    tidyr::separate_wider_delim(
      cols = dplyr::all_of(sample_col),
      delim = sample_field_sep,
      names = sample_fields,
      too_few = "error",
      too_many = "error",
      cols_remove = FALSE
    )

  n_tra <- vapply(cl_wide[[cdr3tra_col]], count_cdr3_values, integer(1))
  n_trb <- vapply(cl_wide[[cdr3trb_col]], count_cdr3_values, integer(1))
  biologically_excluded <- rep(FALSE, nrow(cl_wide))
  if (strict_TCR_biology) {
    biologically_excluded <- n_trb > 1L | n_tra > 2L
  }

  cl_wide$.component_id <- NA_integer_
  assignable <- which(!biologically_excluded & (n_tra + n_trb) > 0L)
  if (length(assignable)) {
    cl_wide$.component_id[assignable] <- assign_shared_cdr3_components(
      cl_wide = cl_wide[assignable, , drop = FALSE],
      sample_col = sample_col,
      donor_col = "donor",
      cdr3tra_col = cdr3tra_col,
      cdr3trb_col = cdr3trb_col,
      shared_cdr3_mode = shared_cdr3_mode,
      same_cl_name_across_donor_only_when_full_annotation =
        same_cl_name_across_donor_only_when_full_annotation
    )
  }

  component_ids <- sort(unique(stats::na.omit(cl_wide$.component_id)))
  message("Picking random names for ", length(component_ids), " CDR3-defined clonotypes.")
  cl_wide$cl_name <- NA_character_
  if (length(component_ids)) {
    pick_randomNames_args1 <- c(
      list(n = length(component_ids), names_to_avoid = names_to_avoid),
      pick_randomNames_args
    )
    component_names <- Gmisc::fastDoCall(
      pick_randomNames,
      pick_randomNames_args1
    )
    cl_wide$cl_name <- component_names[match(cl_wide$.component_id, component_ids)]
  }
  cl_wide9 <- dplyr::select(cl_wide, -.component_id)



  cell_assignments <- cl_wide9 |>
    dplyr::select(
      !!rlang::sym(sample_col),
      !!rlang::sym(barcode_col),
      cl_name
    )

  cl_long9_join <- dplyr::left_join(
    cl_long,
    cell_assignments,
    by = c(sample_col, barcode_col)
  )

  message("Making cl_wide from cl_long.")
  add_cols <- setdiff(names(cl_long9_join), c(sample_col, barcode_col, chain_col, "cl_name"))

  cl_wide9_join <- cl_long9_join |>
    dplyr::select(
      !!rlang::sym(sample_col),
      !!rlang::sym(barcode_col),
      !!rlang::sym(chain_col),
      cl_name,
      dplyr::all_of(add_cols)) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym(chain_col),
      values_from = dplyr::all_of(add_cols),
      values_fn = collapse_fun2) |>
    tidyr::separate(!!rlang::sym(sample_col), into = sample_fields, remove = F, sep = sample_field_sep)

  message("Ordering columns of cl_wide.")
  # use alphabetical order from TRA and TRB to order other columns
  for (i in c(cdr3tra_col, cdr3trb_col)) {
    TRAB_order <- sapply(strsplit(cl_wide9_join[[i]], ","), order, simplify = F)
    for (colname in grep(paste0(gsub(cdr3_col, "", i), "$"), names(cl_wide9_join), value = T)) {
      cl_wide9_join[[colname]] <- purrr::map2_chr(.x = cl_wide9_join[[colname]], .y = TRAB_order, function(x,y) {
        if (is.na(x) || !grepl(",", x, fixed = TRUE)) {
          x
        } else {
          paste(strsplit(x, ",", fixed = TRUE)[[1]][y], collapse = ",")
        }
      })
    }
  }

  for (colname in grep("_TRA$|_TRB$", names(cl_wide9_join), value = T)) {
    cl_wide9_join[which(cl_wide9_join[[colname]] == "NA"),colname] <- NA
  }

  tt <- cl_wide9_join |>
    dplyr::filter(cl_name == "Bigback_Davier") |>
    dplyr::distinct(cdr3_TRA, cdr3_TRB, sample)

  return(list(
    cl_long = cl_long9_join,
    cl_wide = cl_wide9_join,
    multi_anno_TCR_equal_umi_read = multi_anno_TCR_equal_umi_read))
}



split_cdr3_values <- function(x) {
  if (!length(x) || all(is.na(x))) {
    return(character())
  }

  values <- unlist(strsplit(as.character(x[!is.na(x)]), ",", fixed = TRUE))
  values <- trimws(values)
  sort(unique(values[nzchar(values) & values != "NA"]))
}


count_cdr3_values <- function(x) {
  length(split_cdr3_values(x))
}


assign_shared_cdr3_components <- function(
    cl_wide,
    sample_col = "sample",
    donor_col = "donor",
    cdr3tra_col = "cdr3_TRA",
    cdr3trb_col = "cdr3_TRB",
    shared_cdr3_mode = c("compatible", "any"),
    same_cl_name_across_donor_only_when_full_annotation = TRUE) {

  shared_cdr3_mode <- match.arg(shared_cdr3_mode)
  required_cols <- c(sample_col, donor_col, cdr3tra_col, cdr3trb_col)
  missing_cols <- setdiff(required_cols, names(cl_wide))
  if (length(missing_cols)) {
    stop(
      "Missing columns for CDR3 component assignment: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(cl_wide[c(sample_col, donor_col)])) {
    stop("Sample and donor cannot be missing during CDR3 component assignment.")
  }

  n <- nrow(cl_wide)
  if (!n) {
    return(integer())
  }

  tra_sets <- lapply(cl_wide[[cdr3tra_col]], split_cdr3_values)
  trb_sets <- lapply(cl_wide[[cdr3trb_col]], split_cdr3_values)
  tra_signature <- vapply(tra_sets, paste, character(1), collapse = "\035")
  trb_signature <- vapply(trb_sets, paste, character(1), collapse = "\035")
  full_annotation <- lengths(tra_sets) > 0L & lengths(trb_sets) > 0L
  has_annotation <- lengths(tra_sets) > 0L | lengths(trb_sets) > 0L

  parent <- seq_len(n)
  members <- lapply(seq_len(n), function(i) i)

  root_of <- function(i) {
    while (parent[[i]] != i) {
      i <- parent[[i]]
    }
    i
  }

  has_covering_set <- function(sequence_sets) {
    sequence_sets <- sequence_sets[lengths(sequence_sets) > 0L]
    if (length(sequence_sets) < 2L) {
      return(TRUE)
    }

    set_keys <- vapply(sequence_sets, paste, character(1), collapse = "\035")
    sequence_sets <- sequence_sets[!duplicated(set_keys)]
    any(vapply(seq_along(sequence_sets), function(i) {
      candidate <- sequence_sets[[i]]
      all(vapply(sequence_sets, function(observed) {
        all(observed %in% candidate)
      }, logical(1)))
    }, logical(1)))
  }

  component_is_consistent <- function(indices) {
    has_covering_set(tra_sets[indices]) && has_covering_set(trb_sets[indices])
  }

  merge_roots <- function(a, b, guard = shared_cdr3_mode == "compatible") {
    root_a <- root_of(a)
    root_b <- root_of(b)
    if (root_a == root_b) {
      return(root_a)
    }

    combined_members <- c(members[[root_a]], members[[root_b]])
    if (guard && !component_is_consistent(combined_members)) {
      return(NA_integer_)
    }

    keep_root <- min(root_a, root_b)
    drop_root <- max(root_a, root_b)
    parent[[drop_root]] <<- keep_root
    members[[keep_root]] <<- sort(combined_members)
    members[[drop_root]] <<- integer()
    keep_root
  }

  merge_node_group <- function(node_ids, guard = shared_cdr3_mode == "compatible") {
    group_roots <- unique(vapply(node_ids, root_of, integer(1)))
    if (length(group_roots) < 2L) {
      return(TRUE)
    }

    combined_members <- unlist(members[group_roots], use.names = FALSE)
    if (guard && !component_is_consistent(combined_members)) {
      return(FALSE)
    }

    current_root <- group_roots[[1]]
    for (next_root in group_roots[-1]) {
      current_root <- merge_roots(current_root, next_root, guard = FALSE)
    }
    TRUE
  }

  make_key <- function(...) {
    do.call(paste, c(list(...), sep = "\034"))
  }

  merge_by_key <- function(keys, eligible = rep(TRUE, n), guard = TRUE) {
    eligible_nodes <- which(eligible)
    if (length(eligible_nodes) < 2L) {
      return(invisible(NULL))
    }

    key_groups <- split(eligible_nodes, keys[eligible_nodes], drop = TRUE)
    for (node_ids in key_groups[lengths(key_groups) > 1L]) {
      merge_node_group(node_ids, guard = guard)
    }
    invisible(NULL)
  }

  # Exact cell signatures are indisputable within a sample and greatly reduce
  # the graph before shared-sequence links are considered.
  sample_signature_key <- make_key(
    cl_wide[[sample_col]],
    tra_signature,
    trb_signature
  )
  merge_by_key(
    sample_signature_key,
    eligible = has_annotation,
    guard = FALSE
  )

  # Build an inverted index from sample/chain/CDR3 to exact-signature
  # components. This avoids comparing every cell with every other cell.
  n_memberships <- sum(lengths(tra_sets)) + sum(lengths(trb_sets))
  shared_rows <- vector("list", n)
  if (n_memberships) {
    for (i in seq_len(n)) {
      shared_rows[[i]] <- data.frame(
        node = rep(i, length(tra_sets[[i]]) + length(trb_sets[[i]])),
        sample = rep(
          as.character(cl_wide[[sample_col]][[i]]),
          length(tra_sets[[i]]) + length(trb_sets[[i]])
        ),
        chain = c(rep("TRA", length(tra_sets[[i]])), rep("TRB", length(trb_sets[[i]]))),
        cdr3 = c(tra_sets[[i]], trb_sets[[i]]),
        stringsAsFactors = FALSE
      )
    }
    shared_rows <- do.call(rbind, shared_rows)
  } else {
    shared_rows <- data.frame(
      node = integer(),
      sample = character(),
      chain = character(),
      cdr3 = character()
    )
  }

  edge_environment <- new.env(hash = TRUE, parent = emptyenv())
  if (nrow(shared_rows)) {
    shared_key <- make_key(shared_rows$sample, shared_rows$chain, shared_rows$cdr3)
    shared_groups <- split(seq_len(nrow(shared_rows)), shared_key, drop = TRUE)

    for (row_ids in shared_groups[lengths(shared_groups) > 1L]) {
      group_roots <- unique(vapply(shared_rows$node[row_ids], root_of, integer(1)))
      if (length(group_roots) < 2L) {
        next
      }

      root_pairs <- utils::combn(sort(group_roots), 2L)
      for (j in seq_len(ncol(root_pairs))) {
        root_a <- root_pairs[1L, j]
        root_b <- root_pairs[2L, j]
        pair_key <- sprintf("%012d:%012d", root_a, root_b)

        combined_members <- c(members[[root_a]], members[[root_b]])
        if (shared_cdr3_mode == "compatible" &&
            !component_is_consistent(combined_members)) {
          next
        }

        shared_chain <- shared_rows$chain[[row_ids[[1]]]]
        if (base::exists(pair_key, envir = edge_environment, inherits = FALSE)) {
          edge <- base::get(pair_key, envir = edge_environment, inherits = FALSE)
          edge$chains <- unique(c(edge$chains, shared_chain))
        } else {
          edge <- list(
            key = pair_key,
            a = root_a,
            b = root_b,
            chains = shared_chain
          )
        }
        base::assign(pair_key, edge, envir = edge_environment)
      }
    }
  }

  edge_names <- sort(base::ls(edge_environment, all.names = TRUE))
  edges <- if (length(edge_names)) {
    base::mget(edge_names, envir = edge_environment, inherits = FALSE)
  } else {
    list()
  }

  # A TRA-only or TRB-only component must not choose arbitrarily among several
  # incompatible paired targets. All such links are withheld unless the target
  # signatures themselves have a directly observed common covering set.
  blocked_edges <- character()
  if (shared_cdr3_mode == "compatible" && length(edges)) {
    stage_roots <- unique(vapply(seq_len(n), root_of, integer(1)))
    for (root in stage_roots) {
      root_members <- members[[root]]
      missing_tra <- all(lengths(tra_sets[root_members]) == 0L)
      missing_trb <- all(lengths(trb_sets[root_members]) == 0L)
      if (missing_tra == missing_trb) {
        next
      }

      adjacent <- Filter(function(edge) edge$a == root || edge$b == root, edges)
      if (!length(adjacent)) {
        next
      }

      candidate_edges <- Filter(function(edge) {
        target <- if (edge$a == root) edge$b else edge$a
        target_members <- members[[target]]
        if (missing_tra) {
          any(lengths(tra_sets[target_members]) > 0L)
        } else {
          any(lengths(trb_sets[target_members]) > 0L)
        }
      }, adjacent)

      candidate_roots <- unique(vapply(candidate_edges, function(edge) {
        if (edge$a == root) edge$b else edge$a
      }, integer(1)))

      if (length(candidate_roots) > 1L) {
        candidate_members <- unlist(members[candidate_roots], use.names = FALSE)
        if (!component_is_consistent(candidate_members)) {
          blocked_edges <- c(
            blocked_edges,
            vapply(candidate_edges, `[[`, character(1), "key")
          )
        }
      }
    }
    blocked_edges <- unique(blocked_edges)
  }

  if (length(edges)) {
    edge_strength <- vapply(edges, function(edge) length(edge$chains), integer(1))
    edge_order <- order(-edge_strength, names(edges))
    for (edge in edges[edge_order]) {
      if (edge$key %in% blocked_edges) {
        next
      }
      merge_roots(
        edge$a,
        edge$b,
        guard = shared_cdr3_mode == "compatible"
      )
    }
  }

  # Stronger, fully paired exact evidence is propagated before partial exact
  # evidence. Batch consistency checks keep a shared partial signature from
  # selecting one of several contradictory paired components by iteration order.
  donor_signature_key <- make_key(
    cl_wide[[donor_col]],
    tra_signature,
    trb_signature
  )
  global_signature_key <- make_key(tra_signature, trb_signature)
  merge_by_key(
    donor_signature_key,
    eligible = full_annotation,
    guard = shared_cdr3_mode == "compatible"
  )
  merge_by_key(
    global_signature_key,
    eligible = full_annotation,
    guard = shared_cdr3_mode == "compatible"
  )

  if (same_cl_name_across_donor_only_when_full_annotation) {
    merge_by_key(
      donor_signature_key,
      eligible = !full_annotation & has_annotation,
      guard = shared_cdr3_mode == "compatible"
    )
  } else {
    merge_by_key(
      global_signature_key,
      eligible = !full_annotation & has_annotation,
      guard = shared_cdr3_mode == "compatible"
    )
  }

  roots <- vapply(seq_len(n), root_of, integer(1))
  match(roots, unique(roots))
}



collapse_fun <- function(x) paste(x, collapse = ",")

collapse_fun2 <- function(x) {
  if (length(unique(x)) == 1) {
    as.character(unique(x))
  } else {
    paste(x, collapse = ",")
  }
}
collapse_order_fun <- function(x) paste(sort(x), collapse = ",")
collapse_order_fun2 <- function(x) {
  values <- sort(unique(as.character(x[!is.na(x)])))
  values <- trimws(values)
  values <- values[nzchar(values) & values != "NA"]
  if (!length(values)) NA_character_ else paste(values, collapse = ",")
}

check.clonotype.id.levels <- function(cl_wide,
                                      id_cols = c("clonotype_id_TRA", "clonotype_id_TRB"),
                                      group_cols = c("sample")) {

  # maybe this function is useless

  if (length(id_cols) != 2) {
    stop("id_cols has to have length 2.")
  }

  ## test if any clonotype_id between TRA and TRB is different
  cl_wide_sub <-
    cl_wide %>%
    tidyr::drop_na(dplyr::all_of(id_cols)) %>%
    dplyr::mutate(!!id_cols[1] := strsplit(!!rlang::sym(id_cols[1]), ","),
                  !!id_cols[2] := strsplit(!!rlang::sym(id_cols[2]), ",")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!id_cols[1] := paste(unique(!!rlang::sym(id_cols[1])), collapse = ","),
                  !!id_cols[2] := paste(unique(!!rlang::sym(id_cols[2])), collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!rlang::syms(c(id_cols[1], group_cols))) %>%
    dplyr::summarise(n = dplyr::n_distinct(!!rlang::sym(id_cols[1])), .groups = "drop")

  return(cl_wide_sub)
}

collapse.clonotypes <- function(cl_wide,
                                clonotype_col = "cl_name",
                                cdr3_cols = c("cdr3_TRB", "cdr3_TRA"),
                                split_cols = NULL,
                                barcode_col = "barcode") {
  # e.g. c("clonotype_id_TRB", "cdr3_TRB", "sample")
  # e.g. c("clonotype_id_TRA", "cdr3_TRA", "sample")
  # e.g. c("cdr3_TRA", "cdr3_TRB", "patient")
  #cdr3_col_to_unnest = NULL) {

  if (is.null(split_cols)) {
    stop("Provide a value for split_cols.")
  }
  '  if (!is.null(cdr3_col_to_unnest)) {
    cdr3_col_to_unnest <- match.arg(cdr3_col_to_unnest, c("cdr3_TRB", "cdr3_TRA"))
    if (!cdr3_col_to_unnest %in% split_cols) {
      stop("Shouldnt cdr3_col_to_unnest be in split_cols as well?")
    }
  }'

  if (length(cdr3_cols) != 2) {
    stop("cdr3_cols should be of length 2.")
  }

  #split_cols <- match.arg(split_cols, c("clonotype_id_TRB", "clonotype_id_TRA"))

  ## this collapses cl_name to one based on same values in split_cols
  ## so within groups of split_cols, the same cl_name is assigned
  ## in principle this uses the annotation from 10X (clonotype_id) to group same clonotypes together

  ## this can only work within sample, as clonotype_id may be equal for completely different clonotypes across samples
  ## so when cl_wide contains multiple samples, these have to be split by split_cols
  ## or this could be done outside of this function

  ## run this function separately on cdr3_TRB and cdr3_TRA
  ## this probably has benefits as NAs in one of these columns cannot be used for grouping

  ## maybe instead of the above; unnesting cdr3_TRB and cdr3_TRA is needed to --> extra fun
  ## alternative: just join all additional columns from cl_long in the very end
  ## catch all collapsible clonotypes (think of multi annotations)

  ## some clonotypes collapse to one clonotype_id, already by 10X, when there is one multi-annotated clone
  ## which then links two somewhat different groups of clonotypes with different CDR3s

  n_barcode_start <- length(unique(cl_wide[[barcode_col]]))
  barcode_count_start <- sort(table(cl_wide[[barcode_col]]))
  n_unique_cl_start <- length(unique(cl_wide[["cl_name"]]))

  '  if (length(cdr3_col) > 1 && unnest_cdr3_col) {
    message("cdr3_col > 1; setting unnest_cdr3_col to FALSE.")
    unnest_cdr3_col <- F
  }

  if (!is.null(cdr3_col_to_unnest)) {
    cl_wide <-
      cl_wide %>%
      dplyr::mutate(!!cdr3_col_to_unnest := strsplit(!!rlang::sym(cdr3_col_to_unnest), ",")) %>%
      tidyr::unnest(!!rlang::sym(cdr3_col_to_unnest))
  }'

  split_val <- purrr::reduce(lapply(split_cols, function(x) cl_wide[,x,drop=T]), paste, sep = "_")
  cl_wide_split <- split(cl_wide, split_val)
  cl_wide_split_NA <- cl_wide_split[which(grepl("NA_NA", names(cl_wide_split)))]
  cl_wide_split_nonNA <- cl_wide_split[which(!grepl("NA_NA", names(cl_wide_split)))] ##only helper
  # one clonotype annotation from 10X per cl_name
  cl_wide_split_1 <- cl_wide_split_nonNA[sapply(cl_wide_split_nonNA, function(x) length(unique(x[,clonotype_col,drop=T])) == 1)]
  # multiple clonotype annotation from 10X per cl_name
  cl_wide_split_multi <- cl_wide_split_nonNA[sapply(cl_wide_split_nonNA, function(x) length(unique(x[,clonotype_col,drop=T])) > 1)]

  #tt <- cl_wide_split_multi[which(grepl("CAENEGSGGGADGLTF", names(cl_wide_split_multi)))][[1]]

  # check how many factor levels there are in the other cdr3
  # max should be max of number differnet cdr3 annotated per cl (this is max 2 TRA when TCR biology is strict)

  # check this example
  #x<-cl_wide_split_multi[[which(grepl("CASSSLGQPIEKLFF", names(cl_wide_split_multi)))]]
  # a third clonotype that has two TRA should not cause collapsing two other clonotypes to one which have one of its chains
  # only as the first with two TRA annotated is present in the same sample 'by coincidence' causes the collapsing
  # this should not happen

  #which(split_val == "check2")
  #x <- cl_wide_split_multi[[488]]

  other_cdr3_col <- setdiff(cdr3_cols, split_cols)
  if (length(cl_wide_split_multi) > 0) {
    if (length(other_cdr3_col) > 0) {
      split_val <- sapply(cl_wide_split_multi, function (x) {
        #!!!
        # this is a critical step which decides which clonotypes are collapsed
        # how to decide that algorithmically?
        #!!!
        # the order of if clauses below matters!

        other_chains <- unique(x[,other_cdr3_col,drop=T])
        if (length(other_chains) == 1) {
          return(T)
        }
        if (length(which(is.na(other_chains))) == 1 && length(other_chains) == 2) {
          # one NA and other only one level
          # can be collapsed
          return(T)
        }
        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && !any(grepl(",", other_chains))) {
          # more than one other chain but no multi annotation
          return(F)
        }
        if (length(which(is.na(other_chains))) == 0 && length(other_chains) >= 2) {
          # one NA and other only one level
          # can be collapsed
          return(F)
        }
        if (length(other_chains) == 2 && length(which(grepl(",", other_chains))) == 1) {
          return(T)
        }
        other_chains_split <- unlist(strsplit(other_chains, ","))
        if (length(which(is.na(other_chains))) == 1 && length(other_chains[which(!is.na(other_chains))]) == 2 && length(which(grepl(",", other_chains))) == 1 && length(unique(other_chains_split[which(!is.na(other_chains_split))])) == 2) { #
          return(T)
        }
        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && any(grepl(",", other_chains)) && length(unique(other_chains_split[which(!is.na(other_chains_split))])) > 2) { #
          # multiple levels with multi annotation
          return(F)
        }

        other_chains_split2 <- strsplit(other_chains[which(!is.na(other_chains))], ",")
        if (length(which(is.na(other_chains))) == 1 && length(other_chains[which(!is.na(other_chains))]) == 3 && length(which(grepl(",", other_chains))) == 1 && length(unique(other_chains_split[which(!is.na(other_chains_split))])) == 2 && all(unlist(lapply(other_chains_split2[which(lengths(other_chains_split2) == 1)], function(y) y %in% other_chains_split2[[which(lengths(other_chains_split2) == 2)]])))) {
          return(T)
        }

        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && !any(grepl(",", other_chains))) {
          return(F)
        }


        stop("New case: check conditions on collapsing clonotypes.")

        '      max <- max(lengths(strsplit(x[,other_cdr3_col,drop=T], ",")))
      unique_other <- unique(unlist(strsplit(x[,other_cdr3_col,drop=T], ",")))
      unique_other <- unique_other[which(!is.na(unique_other))]
      length(unique_other) <= max'
      })
    } else {
      ## when cdr3_cols is cdr3_TRA and cdr3_TRB
      split_val <- rep(T, length(cl_wide_split_multi))
    }

    cl_wide_split_multi_1 <- cl_wide_split_multi[which(split_val)]
    cl_wide_split_multi_2 <- cl_wide_split_multi[which(!split_val)]
  } else {
    cl_wide_split_multi_1 <- NULL
    cl_wide_split_multi_2 <- NULL
  }


  if (length(cl_wide_split_multi_1) > 0) {
    igsc:::print_table(table(sapply(cl_wide_split_multi_1, function(x) length(unique(x[["cl_name"]])))),
                       c("n unique cl_name", "n groups"),
                       paste0(length(cl_wide_split_multi_1), " groups with collapsible clonotypes."))



    #message("Number of different cdr3_cols combinations per group, including NA:")
    #print(table(sapply(cl_wide_split_multi_1, function(x) nrow(unique(x[,cdr3_cols])))))
    #message(sum(sapply(cl_wide_split_multi_1, nrow)), " total clones to be collapsed.")

    cl_wide_split_multi_1 <- lapply(cl_wide_split_multi_1, function(x) {

      x[,clonotype_col] <- sort(unique(x[,clonotype_col,drop=T]))[1] # alphabetic first
      '
      # make this an outside option?
      if (!is.null(cdr3_col_to_unnest)) {
        x[,clonotype_col] <- names(sort(table(x[,clonotype_col,drop=T]), decreasing = T))[1] # most frequent first
      } else {
        x[,clonotype_col] <- sort(unique(x[,clonotype_col,drop=T]))[1] # alphabetic first
      }'
      return(x)
    })
  }


  cl_wide <- Reduce(dplyr::bind_rows, list(dplyr::bind_rows(cl_wide_split_NA),
                                           dplyr::bind_rows(cl_wide_split_1),
                                           dplyr::bind_rows(cl_wide_split_multi_1),
                                           dplyr::bind_rows(cl_wide_split_multi_2)))
  '
  if (!is.null(cdr3_col_to_unnest)) {
    cl_wide <-
      cl_wide %>%
      dplyr::group_by(!!!rlang::syms(names(cl_wide)[which(!names(cl_wide) %in% cdr3_col_to_unnest)])) %>%
      dplyr::summarise(!!cdr3_col_to_unnest := collapse_order_fun(!!rlang::sym(cdr3_col_to_unnest)), .groups = "drop")

    # in rare cases clonotype with different multi-annotated chains were dragged into different split groups above, only when unnest_cdr3_col = T
    # because of that, they got different cl_name; now they cannot be collapse and then their barcodes are duplicated (two rows for one cell which used to be one row only)

    barcode_count_end <- sort(table(cl_wide$barcode))
    barcode_count_df <-
      stack(barcode_count_start) %>%
      dplyr::rename("values1" = values) %>%
      dplyr::mutate(ind = as.character(ind)) %>%
      dplyr::full_join(stack(barcode_count_end) %>%
                         dplyr::rename("values2" = values) %>%
                         dplyr::mutate(ind = as.character(ind)), by = "ind") %>%
      dplyr::mutate(diff = values2 - values1)
    problematic_barcodes <- barcode_count_df %>% dplyr::filter(diff > 0) %>% dplyr::pull(ind)

    if (any(barcode_count_df$diff < 0)) {
      warning("Check barcode_count_df. One or more values are below 0.")
    }

    cl_wide <- rbind(cl_wide_sub1 <-
                       cl_wide %>%
                       dplyr::filter(!barcode %in% problematic_barcodes),
                     cl_wide_sub2 <-
                       cl_wide %>%
                       dplyr::filter(barcode %in% problematic_barcodes) %>%
                       dplyr::mutate(cl_name = cl_names_read) %>%
                       dplyr::group_by(!!!rlang::syms(names(cl_wide)[which(!names(cl_wide) %in% cdr3_col_to_unnest)])) %>%
                       dplyr::summarise(!!cdr3_col_to_unnest := collapse_order_fun(!!rlang::sym(cdr3_col_to_unnest)), .groups = "drop"))
  }'

  n_unique_cl_end <- length(unique(cl_wide[["cl_name"]]))
  message("Number of unique clonotypes at before and after: ", n_unique_cl_start, ", ", n_unique_cl_end, ". (", n_unique_cl_end-n_unique_cl_start, ", ", round(((n_unique_cl_end-n_unique_cl_start)/n_unique_cl_start)*100, 2), " %)")

  n_barcode_end <- length(unique(cl_wide[[barcode_col]]))
  barcode_count_end <- sort(table(cl_wide[[barcode_col]]))

  # run check functions and notify how many different clonotypes per xx still exist
  check1 <- get.cdr3.levels.per.clonotype(cl_wide = cl_wide)
  igsc:::print_table(table(check1[which(check1$total > 1), "n_cdr3_TRA_cdr3_TRB",drop=T]),
                     c("n comb of cdr3_TRA+cdr3_TRB", "n observation"),
                     caption = "CDR3 levels per cl_name")


  check2 <- get.clname.levels.per.clonotypeid(cl_wide)
  igsc:::print_table(table(dplyr::distinct(check2, cl_names, n)$n),
                     c("n unique cl_name per clonotype id", "n observation"),
                     caption = "cl_name levels per clonotype_id")

  print("\n")
  print("\n")
  print("\n")

  if (n_barcode_end != n_barcode_start) {
    warning("Number of unique barcodes has changed.")
  }
  if (!identical(barcode_count_start, barcode_count_end)) {
    warning("Count of barcodes has changed.")
  }

  return(cl_wide)
}


collapse.clonotypes2 <- function(cl_wide,
                                 clonotype_col = "cl_name",
                                 cdr3_cols = c("cdr3_TRA", "cdr3_TRB"),
                                 grouping_cols = c("patient")) {

  ## this function has become obsolete
  ## it is now achieved by collapse.clonotypes

  # collapse.clonotypes relied on clonotype_id within samples
  # so equal clonotypes across samples could not be identified
  # here grouping relies on equal pairings of cdr3_cols within grouping_cols
  # so across sample collapsing is possible

  # though same collapses are made which I thought would have happened in
  # the first function, collapse.clonotypes
  # but anyway ...

  cl_wide <-
    cl_wide %>%
    dplyr::group_by(!!!rlang::syms(c(grouping_cols, cdr3_cols))) %>%
    dplyr::mutate(cl_names = list(sort(unique(!!rlang::sym(clonotype_col))))) %>%
    dplyr::mutate(cl_names_read = paste(sort(unique(!!rlang::sym(clonotype_col))), collapse = ",")) %>%
    dplyr::ungroup()

  # get groups of cl_name which have same combination of cdr3_cols with grouping_cols
  # but why has this not been catched by collapse.clonotypes?
  grouped_cls <-
    cl_wide %>%
    dplyr::filter(lengths(cl_names) > 1) %>%
    dplyr::distinct(cl_names_read, cl_names)

  if (nrow(grouped_cls) > 1) {
    # this gets intersecting/overlapping groups of cl_names within grouped_cls
    # i.e. Laila is in one pair with Maria but also in another pair with Madalyn
    # they will be joined to a triplet
    # get intersecting vectors 1
    grouped_cls_intersect1 <- purrr::flatten(lapply(grouped_cls$cl_names, function(x) {
      purrr::discard(lapply(grouped_cls$cl_names, function(y) {
        if (length(intersect(x,y)) > 0 && !identical(x,y)) {
          sort(unique(c(x,y)))
        }
      }), is.null)
    }))
    grouped_cls_intersect1 <- grouped_cls_intersect1[!duplicated(grouped_cls_intersect1)]

    # this is run a second time
    # but can we be sure that 2 times are enough to catch all intersections?
    # get intersecting vectors 2
    grouped_cls_intersect2 <- purrr::flatten(lapply(grouped_cls_intersect1, function(x) {
      purrr::discard(lapply(grouped_cls_intersect1, function(y) {
        if (length(intersect(x,y)) > 0 && !identical(x,y)) {
          sort(unique(c(x,y)))
        }
      }), is.null)
    }))
    grouped_cls_intersect2 <- grouped_cls_intersect2[!duplicated(grouped_cls_intersect2)]

    # reduce redundancy
    # i.e. those that were grouped into triplets or quadruples or so: remove their underlying pairs
    grouped_cls_intersect1 <- purrr::discard(lapply(grouped_cls_intersect1, function(x) {
      if (all(sapply(grouped_cls_intersect2, function(y) length(intersect(x,y)) == 0))) {
        return(x)
      } else {
        NULL
      }
    }), is.null)

    grouped_cls_intersect0 <- purrr::discard(lapply(grouped_cls$cl_names, function(x) {
      if (all(sapply(grouped_cls_intersect1, function(y) length(intersect(x,y)) == 0))) {
        return(x)
      } else {
        NULL
      }
    }), is.null)

    grouped_cls_intersect0 <- purrr::discard(lapply(grouped_cls_intersect0, function(x) {
      if (all(sapply(grouped_cls_intersect2, function(y) length(intersect(x,y)) == 0))) {
        return(x)
      } else {
        NULL
      }
    }), is.null)

    # join all non-intersecting groupings
    # these groups of cl_name have overlapping assignment of equal cdr3_cols
    grouped_cls_intersect_final <- c(grouped_cls_intersect0, grouped_cls_intersect1, grouped_cls_intersect2)

    for (i in 1:length(grouped_cls_intersect_final)) {
      temp <-
        cl_wide %>%
        dplyr::filter(cl_name %in% grouped_cls_intersect_final[[i]]) %>%
        dplyr::group_by(!!!rlang::syms(cdr3_cols)) %>%
        dplyr::count()
      ## clonotype_col is only collapsed when per group of cdr3_cols there are max. two different combination: nrow(temp) == 2
      ## and one of the combination has to contain an NA, meaning that one chain has not been annotated
      ## in principle clonotypes are collapsed when a cdr3 has only appeared one other clonotype (only one level/sequence of other chain)
      ## but there are additional ones, where the other chain is missing
      ## then it is somewhat like that the chain is missing due to technical reason and that it confers to that one other level that only exists in the data set

      ## this can be run once across samples
      ## and then again within samples to catch as many as possible
      if (nrow(temp) == 2 && (anyNA(temp[,cdr3_cols[1],drop=T]) || anyNA(temp[,cdr3_cols[2],drop=T]))) {
        cl_wide[which(cl_wide[,clonotype_col,drop=T] %in% grouped_cls_intersect_final[[i]]), clonotype_col] <- sort(grouped_cls_intersect_final[[i]])[1]
      }
    }
  }

  return(cl_wide[,-which(names(cl_wide) %in% c("cl_names"))]) #"cl_names_read"
}

check.clonotype.changes <- function(cl_wide_before,
                                    cl_wide_after,
                                    clonotype_col = "cl_name") {

  combined_cl <-
    stack(table(cl_wide_before[,clonotype_col,drop=T])) %>%
    dplyr::rename("prev" = values) %>%
    dplyr::left_join(stack(table(cl_wide_after[,clonotype_col,drop=T])) %>%
                       dplyr::rename("new" = values), by = "ind") %>%
    dplyr::mutate(new = ifelse(is.na(new), 0, new)) %>% # now missing cl_name (assigned completely to another cl_name) get a 0
    dplyr::mutate(diff = new-prev)

  return(combined_cl)
}

compare.cl.wide.df <- function(cl_wide1,
                               cl_wide2,
                               ref_cols = c("clonotype_id_TRB", "clonotype_id_TRA", "patient"),
                               clonotype_col = "cl_name") {

  cl_wide1_summ <- get.clonotype.levels.per.ref(cl_wide = cl_wide1,
                                                ref_cols = ref_cols,
                                                clonotype_col = clonotype_col)
  names(cl_wide1_summ)[which(names(cl_wide1_summ) %in% c("cl_names_str", "cl_names"))] <- paste0(names(cl_wide1_summ)[which(names(cl_wide1_summ) %in% c("cl_names_str", "cl_names"))], "1")

  cl_wide2_summ <- get.clonotype.levels.per.ref(cl_wide = cl_wide2,
                                                ref_cols = ref_cols,
                                                clonotype_col = clonotype_col)
  names(cl_wide2_summ)[which(names(cl_wide2_summ) %in% c("cl_names_str", "cl_names"))] <- paste0(names(cl_wide2_summ)[which(names(cl_wide2_summ) %in% c("cl_names_str", "cl_names"))], "2")

  cl_wide_summ <-
    cl_wide1_summ %>%
    dplyr::left_join(cl_wide2_summ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cl1_len = length(cl_names1), cl2_len = length(cl_names2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cl_len_diff = cl1_len - cl2_len)

  return(cl_wide_summ)
}

get.clonotype.levels.per.ref <- function(cl_wide,
                                         ref_cols = c("clonotype_id_TRB", "clonotype_id_TRA", "patient"),
                                         clonotype_col = "cl_name") {
  #cl_names_count = stack(table(cl_name)))
  cl_wide %>%
    dplyr::group_by(!!!rlang::syms(ref_cols)) %>%
    dplyr::summarise(cl_names_str = paste(unique(cl_name), collapse = ","),
                     cl_names = list(unique(cl_name)), .groups = "drop")
}

get.cdr3.levels.per.clonotype <- function(cl_wide,
                                          clonotype_col = "cl_name",
                                          cdr3_cols = c("cdr3_TRA", "cdr3_TRB")) {

  # cdr3_cols must be length 2 currently
  ## 3 different comb of cdr3_TRA and cdr3_TRB is possible (NA in TRA and TRB)
  ## 4 is also possible in case of multi annotations

  cl_wide %>%
    tidyr::drop_na(!!rlang::sym(clonotype_col)) %>%
    dplyr::group_by(!!!rlang::syms(clonotype_col)) %>%
    dplyr::summarise(!!paste0("n_", paste(cdr3_cols, collapse = "_")) := dplyr::n_distinct(!!!rlang::syms(cdr3_cols)), # , na.rm = T # triple bang needed to work in n_distinct; like in group_by
                     !!cdr3_cols[1] := dplyr::n_distinct(!!!rlang::syms(cdr3_cols[1]), na.rm = T),
                     !!cdr3_cols[2] := dplyr::n_distinct(!!!rlang::syms(cdr3_cols[2]), na.rm = T),
                     total = dplyr::n())
}

get.clname.levels.per.clonotypeid <- function(cl_wide,
                                              clonotype_col = "cl_name",
                                              clonotype_id_col = c("clonotype_id_TRB", "clonotype_id_TRA"),
                                              grouping_cols = c("sample")) {

  cl_wide %>%
    tidyr::drop_na(!!rlang::sym(clonotype_col)) %>%
    dplyr::group_by(!!!rlang::syms(c(clonotype_id_col, grouping_cols))) %>%
    dplyr::summarise(cl_names = paste(unique(cl_name), collapse = ","), n = dplyr::n_distinct(!!!rlang::syms(clonotype_col)), .groups = "drop")

}


compare.cl.wide.by.barcode <- function(cl_wide1, cl_wide2) {
  cl_wide1 %>%
    dplyr::select(barcode, sample, cl_name) %>%
    dplyr::rename("cl_name1" = cl_name) %>%
    dplyr::left_join(cl_wide2 %>%
                       dplyr::select(barcode, sample, cl_name) %>%
                       dplyr::rename("cl_name2" = cl_name))
}



get.number.of.cell.and.unique.cdr3.per.cl <- function(cl_long) {

  cl_long %>%
    dplyr::group_by(patient, cl_name, chain, cdr3) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(patient, cl_name, chain) %>%
    dplyr::mutate(n_unique_cdr3 = dplyr::n_distinct(cdr3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count_and_ncdr3_greate_1 = count > 1 & n_unique_cdr3 > 1)
}
