#' Assign Memorable Names to T-Cell Receptor Clonotypes
#'
#' @description
#' Converts a long-format T-cell receptor (TCR) annotation table to paired
#' alpha/beta-chain representations, assigns a random human-readable name to
#' each clonotype, and repeatedly merges compatible clonotypes within samples
#' and donors. The result contains both long- and wide-format tables.
#'
#' When `single_chains_only = TRUE`, the function keeps the annotation with the
#' highest UMI count and then the highest read count for each barcode and chain.
#' Remaining ties are resolved deterministically by CDR3 sequence and returned
#' separately for inspection. When `strict_TCR_biology = TRUE`, cells with more
#' than one beta chain or more than two alpha chains are excluded from naming
#' and later restored with `NA` as their clonotype name.
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
#'   chains, and split named clonotypes that still contain multiple TRB CDR3
#'   sequences.
#' @param same_cl_name_across_donor_only_when_full_annotation A single logical
#'   value. If `TRUE`, a clonotype name may be shared across donors only when
#'   both its TRA and TRB CDR3 sequences are annotated.
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
#' This function relies on package-internal helpers `collapse_order_fun`,
#' `collapse_fun2`, `collapse.clonotypes`, and
#' `get.cdr3.levels.per.clonotype`. Random names are generated with
#' [igsc::pick_randomNames()]. Messages report progress while annotations are
#' filtered, named, and collapsed.
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
add_clonotype_names <- function(cl_long,
                                single_chains_only = T,
                                strict_TCR_biology = F,
                                same_cl_name_across_donor_only_when_full_annotation = T,
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

  if (!identical(sort(unique(cl_long[[chain_col]])), c("TRA", "TRB"))) {
    stop("chain_col must contain TRA and TRB only.")
  }

  if (!"donor" %in% sample_fields) {
    stop("sample_fields must have a 'donor' field. if not, then make a dummy.")
  }

  # verify sample_fields are consistent
  cl_long |>
    dplyr::distinct(!!rlang::sym(sample_col)) |>
    #dplyr::mutate(!!rlang::sym(sample_col) := "A_B_C_D_E") |>
    tidyr::separate_wider_delim(
      cols = dplyr::all_of(sample_col),
      delim = sample_field_sep,
      names = sample_fields,
      too_few = "error",
      too_many = "error")

  ident_cols <- c(sample_col, barcode_col, chain_col, cdr3_col)
  cdr3tra_col <- paste0(cdr3_col, "_TRA")
  cdr3trb_col <- paste0(cdr3_col, "_TRB")
  clonotypetra_col <- paste0(clonotype_col, "_TRA")
  clonotypetrb_col <- paste0(clonotype_col, "_TRB")

  multi_anno_TCR_equal_umi_read <- NULL
  if (single_chains_only) {
    if (any(!c("umis", "reads") %in% names(cl_long))) {
      stop("columns 'umis' and 'reads' needed in cl_long when single_chains_only.")
    }
    cl_long <- cl_long |>
      dplyr::filter(umis == max(umis), .by = c(!!rlang::sym(barcode_col), !!rlang::sym(chain_col))) |>
      dplyr::filter(reads == max(reads), .by = c(!!rlang::sym(barcode_col), !!rlang::sym(chain_col))) |>
      dplyr::add_count(!!rlang::sym(barcode_col), name = "countxx")
    if (any(cl_long$countxx>2)) {
      message("some barcodes with multi annotations with equal umis and reads.")
      multi_anno_TCR_equal_umi_read <- dplyr::filter(cl_long, countxx>2)
      cl_long <- dplyr::filter(cl_long, !!rlang::sym(cdr3_col) == sort(!!rlang::sym(cdr3_col))[1],
                               .by = c(!!rlang::sym(barcode_col), !!rlang::sym(chain_col)))
    }
    cl_long <- dplyr::select(cl_long, -countxx)
  }
browser()
  cl_wide_names <-
    cl_long |>
    ## use sample here, not patient, to avoid grouping by same barcode from different samples!!!
    # or make a column "ID" of sample and barcode -- this could be checked for uniqueness within the function first
    # use collapse_order_fun to have the same order for every transcriptome in case of multi annotations
    dplyr::select(dplyr::all_of(ident_cols)) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(chain_col), values_from = !!rlang::sym(cdr3_col), values_fn = collapse_order_fun) |>
    tidyr::nest(data = c(!!rlang::sym(sample_col), !!rlang::sym(barcode_col)))

  if (strict_TCR_biology) {
    cl_wide__TCR_biology_excluded <-
      cl_wide_names |>
      dplyr::filter(stringr::str_count(TRB, ",") > 0 | stringr::str_count(TRA, ",") > 1) |>
      tidyr::unnest(data) |>
      dplyr::mutate(TRA = strsplit(TRA, ","), TRB = strsplit(TRB, ",")) |>
      # keep_empty = T is not strictly needed when there are NAs and not empty values; below drop_na is used to remove NAs
      tidyr::unnest(TRA, keep_empty = T) |>
      tidyr::unnest(TRB, keep_empty = T) |>
      tidyr::pivot_longer(cols = c(TRA, TRB), names_to = chain_col, values_to = cdr3_col) |>
      tidyr::drop_na() |>
      dplyr::distinct() |>
      dplyr::mutate(cl_name = NA)
    cl_wide_names <- cl_wide_names |>
      dplyr::filter(stringr::str_count(TRB, ",") == 0 | is.na(TRB)) |>
      dplyr::filter(stringr::str_count(TRA, ",") <= 1 | is.na(TRA))
  }
  message("Picking random names.")
  pick_randomNames_args1 <- c(list(n = nrow(cl_wide_names), names_to_avoid = names_to_avoid), pick_randomNames_args)
  cl_wide_names <- cl_wide_names |>
    dplyr::mutate(cl_name = Gmisc::fastDoCall(igsc::pick_randomNames, pick_randomNames_args1)) |>
    tidyr::unnest(data)
  #any(cl_wide_names[["cl_name"]] %in% names_to_avoid)

  cl_long_join <- cl_wide_names |>
    dplyr::mutate(TRA = strsplit(TRA, ","), TRB = strsplit(TRB, ",")) |>
    # keep_empty = T is probably strictly not needed when there are NAs and not empty values; below drop_na is used to remove NAs
    tidyr::unnest(TRA, keep_empty = T) |>
    tidyr::unnest(TRB, keep_empty = T) |>
    tidyr::pivot_longer(cols = c(TRA, TRB), names_to = chain_col, values_to = cdr3_col) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::left_join(cl_long, by = ident_cols) # to get clonotype_id

  cl_wide <-
    cl_long_join |>
    dplyr::select(
      cl_name,
      !!rlang::sym(chain_col),
      !!rlang::sym(barcode_col),
      !!rlang::sym(sample_col),
      !!rlang::sym(cdr3_col),
      !!rlang::sym(clonotype_col)
    ) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym(chain_col),
      values_from = c(!!rlang::sym(cdr3_col), !!rlang::sym(clonotype_col)),
      names_sep = "_",
      values_fn = collapse_fun2
    ) |>
    tidyr::separate(!!rlang::sym(sample_col), into = sample_fields, remove = F, sep = sample_field_sep)

  # length(unique(cl_wide$cl_name))
  # sum(cl_wide2$cl_name == "Kirby_Parker")
  #
  # tt <- cl_wide2 |> dplyr::filter(cl_name == "Kirby_Parker") |> dplyr::distinct()
  # tt2 <- cl_wide3 |> dplyr::filter(cl_name == "Kirby_Parker") |> dplyr::distinct()

browser()
anyNA(cl_wide$clonotype_id_TRA)
any(cl_wide$clonotype_id_TRA != cl_wide$clonotype_id_TRB)
cl_wide <- cl_wide |> dplyr::mutate(diff = clonotype_id_TRA == clonotype_id_TRB)
any(!cl_wide$diff, na.rm = T)

tt <- cl_wide |>
  tidyr::drop_na(clonotype_id_TRB, clonotype_id_TRA) |>
  dplyr::summarise(n_TRB_id = dplyr::n_distinct(clonotype_id_TRB), .by = cl_name)

  ## this is not beautiful but it worked somehow;
  ## toggling strict_TCR_biology between T and F gave comparable results for clonotypes which are not actually affected by 'forbidden' multi-annotations
  message("Collapsing clonotypes.")
  cl_wide_start <- cl_wide

  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide2 <- collapse.clonotypes(cl_wide, barcode_col = barcode_col, split_cols = c(cdr3trb_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide3 <- collapse.clonotypes(cl_wide2, barcode_col = barcode_col, split_cols = c(cdr3tra_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide4 <- collapse.clonotypes(cl_wide3, barcode_col = barcode_col, split_cols = c(cdr3trb_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide5 <- collapse.clonotypes(cl_wide4, barcode_col = barcode_col, split_cols = c(cdr3tra_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRA+cdr3_TRB and donor ----------------------------")
  cl_wide6 <- collapse.clonotypes(cl_wide5, barcode_col = barcode_col, split_cols = c(cdr3tra_col, cdr3trb_col, "donor"), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide7 <- collapse.clonotypes(cl_wide6, barcode_col = barcode_col, split_cols = c(cdr3trb_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide8 <- collapse.clonotypes(cl_wide7, barcode_col = barcode_col, split_cols = c(cdr3tra_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  message("---------------------------- cdr3_TRA+cdr3_TRB and donor ----------------------------")
  cl_wide9 <- collapse.clonotypes(cl_wide8, barcode_col = barcode_col, split_cols = c(cdr3tra_col, cdr3trb_col, "donor"), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  # run this at very last, once
  message("---------------------------- cdr3_TRA+cdr3_TRB and sample ----------------------------")
  cl_wide9 <- collapse.clonotypes(cl_wide9, barcode_col = barcode_col, split_cols = c(cdr3tra_col, cdr3trb_col, sample_col), cdr3_cols = c(cdr3tra_col, cdr3trb_col))
  #any(cl_wide9[["cl_name"]] %in% names_to_avoid)

  # only allow same clonotype name across patients, when TRA and TRB are annotated
  if (same_cl_name_across_donor_only_when_full_annotation) {
    shared_cl_name <- cl_wide9 |>
      dplyr::summarise(n_donor = dplyr::n_distinct(donor), .by = cl_name) |>
      dplyr::filter(n_donor > 1) |>
      dplyr::pull(cl_name)

    if (length(shared_cl_name)) {
      cl_wide9_split <- split(cl_wide9, cl_wide9[["cl_name"]] %in% shared_cl_name)
      if ("TRUE" %in% names(cl_wide9_split)) {
        cl_wide9_split_TRUE <- split(cl_wide9_split[["TRUE"]], cl_wide9_split[["TRUE"]][["cl_name"]])
        for (i in seq_along(cl_wide9_split_TRUE)) {
          if (anyNA(cl_wide9_split_TRUE[[i]][[cdr3tra_col]]) || anyNA(cl_wide9_split_TRUE[[i]][[cdr3trb_col]])) {
            cl_wide9_split_TRUE[[i]] <- split(cl_wide9_split_TRUE[[i]], cl_wide9_split_TRUE[[i]][["donor"]])
          }
        }

        cl_wide9_split_TRUE <- purrr::list_flatten(cl_wide9_split_TRUE)
        pick_randomNames_args2 <- c(list(n = length(cl_wide9_split_TRUE),
                                         names_to_avoid = c(cl_wide9[["cl_name"]], names_to_avoid)),
                                    pick_randomNames_args)
        new_names <- Gmisc::fastDoCall(igsc::pick_randomNames, pick_randomNames_args2)

        for (i in seq_along(cl_wide9_split_TRUE)) {
          cl_wide9_split_TRUE[[i]][["cl_name"]] <- new_names[i]
        }
        cl_wide9_split_TRUE <- dplyr::bind_rows(cl_wide9_split_TRUE)
        cl_wide9 <- dplyr::bind_rows(cl_wide9_split_TRUE, cl_wide9_split[["FALSE"]])
      }
    }
  }



  # run and notify of result
  # get.clonotype.levels.per.ref
  # get.cdr3.levels.per.clonotype
  # get.clname.levels.per.clonotypeid
  # prep cl_wide for joining

  cl_long9 <-
    cl_wide9 |>
    dplyr::select(-c(!!rlang::sym(clonotypetrb_col), !!rlang::sym(clonotypetra_col))) |>
    dplyr::mutate(!!rlang::sym(cdr3tra_col) := strsplit(!!rlang::sym(cdr3tra_col), ","),
                  !!rlang::sym(cdr3trb_col) := strsplit(!!rlang::sym(cdr3trb_col), ",")) |>
    ## very important to set keep_empty = T; otherwise empty rows (e.g. no !!rlang::sym(cdr3tra_col)) are dropped
    tidyr::unnest(!!rlang::sym(cdr3tra_col), keep_empty = T) |>
    tidyr::unnest(!!rlang::sym(cdr3trb_col), keep_empty = T) |>
    tidyr::pivot_longer(cols = c(!!rlang::sym(cdr3tra_col), !!rlang::sym(cdr3trb_col)),
                        names_to = chain_col,
                        values_to = cdr3_col) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::select(-dplyr::all_of(sample_fields)) |>
    dplyr::mutate(!!rlang::sym(chain_col) := gsub(paste0(cdr3_col, "_"), "", !!rlang::sym(chain_col)))

  if (strict_TCR_biology) {
    # run get.cdr3.levels.per.clonotype(cl_wide9)
    # then split those clonotypes which have multiple TRB
    multi_TRB <-
      get.cdr3.levels.per.clonotype(cl_wide9) |>
      dplyr::filter(!!rlang::sym(cdr3trb_col) > 1) |>
      dplyr::pull(cl_name)
    if (length(multi_TRB) > 0) {
      cl_wide9_1 <- dplyr::filter(cl_wide9, !cl_name %in% multi_TRB) #| is.na(cl_name)
      cl_wide9_2 <- dplyr::filter(cl_wide9, cl_name %in% multi_TRB)
      split_val <- purrr::reduce(lapply(unique(c("cl_name", cdr3trb_col)), function(x) cl_wide9_2[[x]]), paste, sep = "_")
      cl_wide9_2_split <- split(cl_wide9_2, split_val)
      new_names <- igsc::pick_randomNames(n = length(cl_wide9_2_split), names_to_avoid = unique(c(cl_wide9_1[["cl_name"]], names_to_avoid)), max_iter = 10000)
      for (i in seq_along(cl_wide9_2_split)) {
        cl_wide9_2_split[[i]]$cl_name <- new_names[i]
      }
      cl_wide9_2 <- dplyr::bind_rows(cl_wide9_2_split)
      cl_wide9 <- dplyr::bind_rows(cl_wide9_1, cl_wide9_2)
    }
    cl_wide9 <- dplyr::bind_rows(cl_wide9, cl_wide__TCR_biology_excluded)
  }

  cl_long9_join <- dplyr::left_join(cl_long, cl_long9, by = ident_cols)

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

  if (strict_TCR_biology) {
    if (length(intersect(cl_wide9_join |> dplyr::filter(is.na(cl_name)) |> dplyr::distinct(!!rlang::sym(barcode_col)) |> dplyr::pull(!!rlang::sym(barcode_col)),
                         cl_long9_join |> dplyr::filter(is.na(cl_name)) |> dplyr::distinct(!!rlang::sym(barcode_col)) |> dplyr::pull(!!rlang::sym(barcode_col)))) !=
        length(cl_wide9_join |> dplyr::filter(is.na(cl_name)) |> dplyr::distinct(!!rlang::sym(barcode_col)) |> dplyr::pull(!!rlang::sym(barcode_col)))) {
      warning("barcodes associated with is.na(cl_name) are not equal between cl_wide and cl_long.")
    }
  }

  #cl_nest and cl_wide_join should have same nrow

  # cl_nest = cl_long9_join |>
  #   tidyr::nest(data = c(chain, cdr3, cdr3_nt, contig_id, v_gene, d_gene, j_gene, c_gene,
  #                        reads, umis, clonotype_id, consensus_id, contig_seq, consensus_seq, V_imgt, J_imgt, length))

  return(list(# cl_wide_start = cl_wide_start,
    cl_long = cl_long9_join,
    cl_wide = cl_wide9_join,
    multi_anno_TCR_equal_umi_read = multi_anno_TCR_equal_umi_read))
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
  if (length(unique(x)) == 1) {
    as.character(unique(x))
  } else {
    paste(sort(x), collapse = ",")
  }
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

  browser()
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

