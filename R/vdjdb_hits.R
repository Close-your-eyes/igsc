#' Compare CDR3 of TCR from sequencing to known CDR3 in the vdjdb
#'
#' Matching CDR3 with references known antigen specificity may help to
#' narrow potential antigen specificity. Here this is done simply by comparing
#' aa sequences with levensthein distance (lv). Export a tsv of all or selected
#' entries here https://vdjdb.cdr3.net/search. Provide that tsv as data frame
#' (vdjdb). From vdjdb and tcrs distinct rows are compared by lv and if selected
#' by amino acid similarity (BLOSUM62_similarity). Returned data frame may be
#' joined to original columns from vdjdb and tcrs followed by plotting with
#' ggplot2 (see example).
#'
#' @param tcrs data frame with CDR3 seq data; e.g. cl_long from clonotype prep
#' or see system.file("extdata", "example_TCR.tsv.gz", package = "igsc")
#' @param vdjdb data frame of TCR data from vdjdb;
#' if not provided system.file("extdata", "vdjdb.tsv.gz", package = "igsc")
#' is used (downloaded 23/09/2025)
#' @param vdj_tr_col column with TCR chain in the vdjdb reference data frame
#' (should contain TRA and/or TRB only);
#' if NULL: guessed
#' @param tcr_tr_col column with TCR chain in the tcrs data frame
#' (should contain TRA and/or TRB only)
#' if NULL: guessed
#' @param vdj_cdr3_col column with CDR3s in the vdjdb reference data frame
#' if NULL: guessed
#' @param tcr_cdr3_col column with CDR3s in the tcrs data frame
#' if NULL: guessed
#' @param TRAxTRB FALSE: compare TRA vs TRA and TRB vs TRB only,
#' TRUE: TRA vs TRB only;
#' provide T, F or c(T,F); the latter does both
#' @param max_lvdist maximum levensthein distance between CDR3s for return
#' (also filtering before BLOSUM62 similarity calculation)
#' @param BLOSUM62_similarity calculate similarity (pairwise alignment) between
#' CDR3s based on BLOSUM62 substitution matrix
#' @param ... arguments passed to mapply_fun and lapply_fun,
#' most relevant: mc.cores (e.g parallel::detectCores())
#' @param lapply_fun lapply function for levensthein distance calculation,
#' suggested are lapply, pblapply::pblapply, parallel::mclapply;
#' however stringdist package uses multiple threads: getOption("sd_num_thread")
#'
#' @return data frame of matched CDR3 from tcrs and vdjdb
#' @export
#'
#' @examples
#' \dontrun{
#' # tcrs data frame with sequencing data from TCRs has to be created (cl_long is appropriate)
#' vdjdb_complete <- as.data.frame(vroom::vroom(system.file("extdata", "vdjdb.tsv.gz", package = "igsc"), .name_repair = make.names))
#' tcr_data <- as.data.frame(vroom::vroom(system.file("extdata", "example_TCR.tsv.gz", package = "igsc")))
#' vdjdb_hs <- vdjdb_complete |>
#'   dplyr::filter(Species == "HomoSapiens") |>
#'   dplyr::rename("chain" = Gene)
#' # dont differ MHC (CD4 vs CD8 origin) here, but you may
#' uni_res <- igsc::vdjdb_hits(vdjdb = vdjdb_hs,
#'                             tcrs = tcr_data,
#'                             TRAxTRB = c(F, T), # TRA vs TRA, TRB vs TRB and TRA vs TRB
#'                             max_lvdist = 8) # max distance (max n of unequal aa in CDR3)
#' res <- uni_res |>
#'   dplyr::left_join(vdjdb_hs, by = c("chain_vdjdb" = "chain", "CDR3_vdjdb" = "CDR3", "MHC.class" = "MHC.class"))
#' # filter for most reliable hits
#' res2 <- res |>
#'   #dplyr::filter(Score>1) |> # only most confident entries from vdjdb, see: https://github.com/antigenomics/vdjdb-db
#'   #dplyr::filter(lv == 0) |> # exact matches only
#'   dplyr::group_by(CDR3_vdjdb, CDR3_tcr, lv, chain_tcr, chain_vdjdb, MHC.class, Epitope, Epitope.gene, Epitope.species) |>
#'   dplyr::summarise(Score = paste(unique(Score), collapse = ","), count = dplyr::n(), .groups = "drop") |>
#'   dplyr::add_count(CDR3_tcr) # some CDR3 may have multiple Epitope hits
#'
#' plot <- ggplot(res2, aes(x = Epitope.gene)) +
#'   geom_bar(aes(fill = Epitope.species), color = "black") +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#'   colrr::scale_fill_custom() +
#'   facet_wrap(vars(MHC.class), scales = "free")
#' }
vdjdb_hits <- function(tcrs,
                       vdjdb = NULL,
                       vdj_tr_col = NULL,
                       tcr_tr_col = NULL,
                       vdj_cdr3_col = NULL,
                       tcr_cdr3_col = NULL,
                       TRAxTRB = F,
                       max_lvdist = 3,
                       BLOSUM62_similarity = F,
                       lapply_fun = lapply,
                       ...) {


  if (missing(tcrs)) {
    stop("Please provide a tcrs data frame.")
  }

  if (is.null(vdjdb)) {
    vdjdb <- as.data.frame(vroom::vroom(system.file("extdata", "vdjdb.tsv.gz", package = "igsc"),
                                        show_col_types = F,
                                        progress = F,
                                        .name_repair = make.names))
  }


  # detect column names
  if (is.null(vdj_tr_col)) {
    vdj_tr_col <- names(which.max(prop_TRA_TRB(vdjdb)))
  }
  if (is.null(vdj_cdr3_col)) {
    vdj_cdr3_col <- names(which.max(prop_CDR3_start(vdjdb)))
  }
  if (is.null(tcr_tr_col)) {
    tcr_tr_col <- names(which.max(prop_TRA_TRB(tcrs)))
  }
  if (is.null(tcr_cdr3_col)) {
    tcr_cdr3_col <- names(which.max(prop_CDR3_start(tcrs)))
  }

  # check if columns exist
  if (any(!c(vdj_tr_col, vdj_cdr3_col) %in% names(vdjdb))) {
    stop("vdj_tr_col or vdj_cdr3_col not found in vdjdb.")
  }
  if (any(!c(tcr_tr_col, tcr_cdr3_col) %in% names(tcrs))) {
    stop("tcr_tr_col or tcr_cdr3_col not found in tcrs.")
  }

  if (vdj_tr_col == tcr_tr_col) {
    message("Making tcr_tr_col and vdj_tr_col unique.")
    names(vdjdb)[which(names(vdjdb) == vdj_tr_col)] <- paste0(vdj_tr_col, "_vdjdb")
    vdj_tr_col <- paste0(vdj_tr_col, "_vdjdb")
    names(tcrs)[which(names(tcrs) == tcr_tr_col)] <- paste0(tcr_tr_col, "_tcr")
    tcr_tr_col <- paste0(tcr_tr_col, "_tcr")
  }

  if (vdj_cdr3_col == tcr_cdr3_col) {
    message("Making tcr_cdr3_col and vdj_cdr3_col unique.")
    names(vdjdb)[which(names(vdjdb) == vdj_cdr3_col)] <- paste0(vdj_cdr3_col, "_vdjdb")
    vdj_cdr3_col <- paste0(vdj_cdr3_col, "_vdjdb")
    names(tcrs)[which(names(tcrs) == tcr_cdr3_col)] <- paste0(tcr_cdr3_col, "_tcr")
    tcr_cdr3_col <- paste0(tcr_cdr3_col, "_tcr")
  }

  message("tcr_tr_col: ", tcr_tr_col)
  message("vdj_tr_col: ", vdj_tr_col)
  message("tcr_cdr3_col: ", tcr_cdr3_col)
  message("vdj_cdr3_col: ", vdj_cdr3_col)


  vdjdb[[vdj_tr_col]] <- toupper(vdjdb[[vdj_tr_col]])
  vdjdb[[vdj_cdr3_col]] <- toupper(vdjdb[[vdj_cdr3_col]])
  tcrs[[tcr_tr_col]] <- toupper(tcrs[[tcr_tr_col]])
  tcrs[[tcr_cdr3_col]] <- toupper(tcrs[[tcr_cdr3_col]])


  if (!all(unique(vdjdb[[vdj_tr_col]]) %in% c("TRA", "TRB"))) {
    stop(paste0(vdj_tr_col, " column should only contain TRA and/or TRB."))
  }
  if (!all(unique(tcrs[[tcr_tr_col]]) %in% c("TRA", "TRB"))) {
    stop(paste0(tcr_tr_col, " column should only contain TRA and/or TRB."))
  }

  lapply_fun <- match.fun(lapply_fun)

  vdjdb <- dplyr::distinct(vdjdb, !!rlang::sym(vdj_cdr3_col), !!rlang::sym(vdj_tr_col))
  tcrs <- dplyr::distinct(tcrs, !!rlang::sym(tcr_tr_col), !!rlang::sym(tcr_cdr3_col))

  matches <- do.call(rbind, lapply(TRAxTRB,
                                   vdjdb_tcrs_match_fun,
                                   vdjdb = vdjdb,
                                   tcrs = tcrs,
                                   vdj_tr_col = vdj_tr_col,
                                   tcr_tr_col = tcr_tr_col,
                                   vdj_cdr3_col = vdj_cdr3_col,
                                   tcr_cdr3_col = tcr_cdr3_col,
                                   max_lvdist = max_lvdist,
                                   lapply_fun = lapply_fun,
                                   ...))


  matches <- dplyr::bind_cols(matches,
                              get_alignment_scores(x = matches[[vdj_cdr3_col]],
                                                   y = matches[[tcr_cdr3_col]],
                                                   substitution_matrix = "BLOSUM62",
                                                   ...)
  )



  return(matches)

}


#' Calculate Pairwise Alignment Scores
#'
#' Aligns corresponding elements of `x` and `y` and optionally normalizes each
#' alignment score by the mean of the two sequences' self-alignment scores. If
#' `y` is `NULL`, all ordered pairs of sequences in `x` are aligned.
#'
#' @param x A character vector or sequence-set object containing the first set
#'   of biological sequences. It must not be empty.
#' @param y A character vector or sequence-set object containing the second set
#'   of biological sequences, or `NULL`. When supplied, it must have the same
#'   length as `x`. When `NULL`, all ordered pairs in `x` are compared.
#' @param substitution_matrix A substitution matrix or a single character
#'   string naming a matrix accepted by [pwalign::pairwiseAlignment()]. Defaults
#'   to `"BLOSUM62"`.
#' @param score_col A single, non-missing character string used as the name of
#'   the raw-score result. When `NULL`, the substitution-matrix name is used if
#'   available; otherwise `"score"` is used.
#' @param normalized_col A single, non-missing character string used as the
#'   name of the normalized-score result. When `NULL`, `"_rel"` is appended to
#'   `score_col`.
#' @param return Output format: `"mat"` returns a named list of matrices and
#'   `"df"` returns a data frame. In matrix output, rows represent `x` sequences
#'   and columns represent `y` sequences; unobserved combinations are `NA`.
#' @param return_score One or both of `"score"` and `"normalized"`, specifying
#'   which score types to return. Defaults to both.
#' @param ... Additional arguments passed to [pwalign::pairwiseAlignment()],
#'   such as `type`, `gapOpening`, and `gapExtension`. `scoreOnly` and
#'   `substitutionMatrix` are set internally and must not be supplied.
#'
#' @return If `return = "df"`, a data frame containing `x`, `y`, and the
#'   requested score columns. If `return = "mat"`, a named list containing one
#'   matrix per requested score type. A normalized score is the pairwise score
#'   divided by the mean of the two corresponding self-alignment scores. It is
#'   `NA` when that mean is zero.
#'
#' @export
#'
#' @examples
#' get_alignment_scores(
#'   x = c("CASSLGQETQYF", "CASSPGQGDNEQFF"),
#'   y = c("CASSIRSSYEQYF", "CASSPGQGTNEQFF"),
#'   return = "df"
#' )
#'
#' get_alignment_scores(
#'   x = c("CASSLGQETQYF", "CASSIRSSYEQYF"),
#'   return = "mat",
#'   return_score = "normalized",
#'   gapOpening = 10,
#'   gapExtension = 0.5
#' )
get_alignment_scores <- function(
    x,
    y = NULL,
    substitution_matrix = "BLOSUM62",
    score_col = substitution_matrix,
    normalized_col = paste0(score_col, "_rel"),
    return = c("mat", "df"),
    return_score = c("score", "normalized"),
    make_unique = F,
    ...) {

  return <- rlang::arg_match(return)
  return_score <- rlang::arg_match(return_score, multiple = T)

  x <- as.character(x)
  if (!length(x) || anyNA(x)) {
    rlang::abort("`x` must contain at least one non-missing sequence.")
  }

  if (!is.null(y)) {
    y <- as.character(y)
    if (!length(y) || anyNA(y)) {
      rlang::abort("`y` must contain at least one non-missing sequence.")
    }
  }

  if (!is.null(y) && anyDuplicated(data.frame(x,y)) && make_unique) {
    message("making x,y pairwise unique.")
    z <- unique(data.frame(x,y))
    x <- z$x
    y <- z$y
  }

  # if (anyDuplicated(x)) {
  #   message("making x unique.")
  #   x <- unique(x)
  # }
  # if (!is.null(y) && anyDuplicated(y)) {
  #   message("making y unique.")
  #   y <- unique(y)
  # }

  for (arg in c("score_col", "normalized_col")) {
    value <- get(arg)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      rlang::abort(paste0("`", arg, "` must be a single, non-empty string."))
    }
  }
  if (identical(score_col, normalized_col)) {
    rlang::abort("`score_col` and `normalized_col` must be different.")
  }

  if (is.null(y)) {
    grid <- expand.grid(x = x, y = x, stringsAsFactors = FALSE)
    x <- grid$x
    y <- grid$y
  } else {
    y <- as.character(y)
    if (length(x) != length(y)) {
      rlang::abort("`x` and `y` must have the same length.")
    }
    if (anyNA(y)) {
      rlang::abort("`y` must not contain missing sequences.")
    }
  }

  score <- pwalign::pairwiseAlignment(
    x,
    y,
    substitutionMatrix = substitution_matrix,
    scoreOnly = TRUE,
    ...
  )

  # ref1 <- pwalign::pairwiseAlignment(
  #   x,
  #   x,
  #   substitutionMatrix = substitution_matrix,
  #   scoreOnly = TRUE,
  #   ...)
  #
  # ref2 <- pwalign::pairwiseAlignment(
  #   y,
  #   y,
  #   substitutionMatrix = substitution_matrix,
  #   scoreOnly = TRUE,
  #   ...)

  if ("normalized" %in% return_score) {
    reference_score <- (pwalign::pairwiseAlignment(
      x,
      x,
      substitutionMatrix = substitution_matrix,
      scoreOnly = TRUE,
      ...) + pwalign::pairwiseAlignment(
        y,
        y,
        substitutionMatrix = substitution_matrix,
        scoreOnly = TRUE,
        ...)) / 2
  }

  if (return == "df") {
    if (length(return_score) == 2) {
      df <- data.frame(score, score / reference_score, x, y)
      names(df) <- c(score_col, normalized_col, "x", "y")
    } else if (return_score == "score") {
      df <- data.frame(score, x, y)
      names(df) <- c(score_col, "x", "y")
    } else if (return_score == "normalized") {
      df <- data.frame(score / reference_score, x, y)
      names(df) <- c(normalized_col, "x", "y")
    }

    return(df)
  }

  if (return == "mat") {
    if (length(return_score) == 2) {
      mats <- list(brathering::df_long_to_mat(data.frame(y,x,score)),
                   brathering::df_long_to_mat(data.frame(y,x,score / reference_score)))
    } else if (return_score == "score") {
      mats <- list(brathering::df_long_to_mat(data.frame(y,x,score)),
                   NULL)
    } else if (return_score == "normalized") {
      mats <- list(NULL,
                   brathering::df_long_to_mat(data.frame(y,x,score / reference_score)))
    }

    names(mats) <- c(score_col, normalized_col)
    return(mats)
  }

}

vdjdb_tcrs_match_fun <- function(TRAxTRB,
                                 vdjdb,
                                 tcrs,
                                 vdj_tr_col,
                                 tcr_tr_col,
                                 vdj_cdr3_col,
                                 tcr_cdr3_col,
                                 max_lvdist,
                                 lapply_fun,
                                 ...) {

  dots <- list(...)
  if (is.null(dots[["mc.cores"]])) {
    dots[["mc.cores"]] <- 1
  }

  # split into chunks of 100 rows
  tcrs <- split_min100(tcrs, n = dots[["mc.cores"]])
  tcrs <- purrr::map(tcrs, ~split(.x[[tcr_cdr3_col]], .x[[tcr_tr_col]])[c("TRA", "TRB")])
  tcrs <- purrr::map(tcrs, ~purrr::discard(.x, lengths(.x) == 0))
  tcrs <- purrr::list_flatten(tcrs)

  vdjdb <- split(vdjdb[[vdj_cdr3_col]], vdjdb[[vdj_tr_col]])[c("TRA", "TRB")]
  vdjdb <- purrr::discard(vdjdb, lengths(vdjdb) == 0)

  ## comparing TRA with TRA and TRB with TRB
  matches <- purrr::map_dfr(purrr::set_names(c("TRA", "TRB")), function(TRX) {
    tcrs <- tcrs[which(grepl(TRX, names(tcrs)))]
    # pick TRX or its counterpart
    vdjdb <- vdjdb[[ifelse(TRAxTRB, setdiff(c("TRA", "TRB"), TRX), TRX)]]

    matches <- lapply_fun(tcrs, function(x) {
      matches <- stringdist::stringdistmatrix(
        a = vdjdb,
        b = x,
        method = "lv",
        useNames = T,
        nthread = getOption("sd_num_thread"))
      userows <- which(apply(matches, 1, function(x) any(x<=max_lvdist)))
      usecols <- which(apply(matches, 2, function(x) any(x<=max_lvdist)))

      brathering::mat_to_df_long(
        matches[userows, usecols, drop = F],
        rownames_to = vdj_cdr3_col,
        colnames_to = tcr_cdr3_col,
        values_to = "lv") |>
        dplyr::filter(lv <= max_lvdist)
    }, ...) |>
      dplyr::bind_rows()
    matches[[tcr_tr_col]] <- TRX
    matches[[vdj_tr_col]] <- ifelse(TRAxTRB, setdiff(c("TRA", "TRB"), TRX), TRX)
    return(matches)
  })
  rownames(matches) <- NULL

  # matches <- lapply_fun(tcrs, function(x) {
  #   #print(which(mapply(identical, tcrs, list(x))))
  #   matches <- mapply(stringdist::stringdistmatrix,
  #                     a = vdjdb,
  #                     b = x,
  #                     method = "lv",
  #                     useNames = "strings",
  #                     SIMPLIFY = F,
  #                     nthread = 1)
  #
  #   matches <- sapply(matches, function(y) y[which(apply(y, 1, min) <= max_lvdist), which(apply(y, 2, min) <= max_lvdist), drop = F], simplify = F)
  #   matches <- sapply(matches, reshape2::melt, simplify = F, c(vdj_cdr3_col, tcr_cdr3_col), value.name = "lv")
  #   names(matches) <- paste(sort(unique(vdjdb[,vdj_tr_col,drop=T]), decreasing = sort_desc), sort(unique(x[,tcr_tr_col,drop=T])), sep = "_")
  #   matches <- do.call(rbind, matches)
  #   matches <- matches[which(matches$lv <= max_lvdist),]
  #   if (nrow(matches) == 0) {
  #     return(NULL)
  #   }
  #   matches[,vdj_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 1)
  #   matches[,tcr_tr_col] <- sapply(strsplit(sapply(strsplit(rownames(matches), "\\."), "[", 1), "_"), "[", 2)
  #   rownames(matches) <- NULL
  #   return(matches)
  # }, ...)
  #
  # matches <- do.call(rbind, matches)
  # rownames(matches) <- NULL

  return(matches)
}

prop_TRA_TRB <- function(df, minrow = 100) {
  rows <- min(minrow, nrow(df))
  sapply(df, function(col) mean(toupper(col[1:rows]) %in% c("TRA", "TRB"), na.rm = TRUE))
}


prop_CDR3_start <- function(df, minrow = 100) {
  CDR3_start <- readRDS(system.file("extdata", "cdr3starts.rds", package = "igsc"))
  rows <- min(minrow, nrow(df))

  charcols <- which(unlist(lapply(df, is.character)))
  candidatecols <- which(unlist(lapply(df[,charcols], function(x) mean(nchar(x))))>9)
  candidatecols2 <- intersect(candidatecols, charcols)
  if (!length(candidatecols2)) {
    stop("CDR3 column could not be guessed.")
  }
  sapply(df[,candidatecols2,drop = F],
         function(col) mean(grepl(paste(paste0("^", names(CDR3_start)), collapse = "|"), col, ignore.case = T), na.rm = T))
}

split_min100 <- function(df, n) {

  # tcrs_split <- split(tcrs, ceiling(1:nrow(tcrs)/100))
  total <- nrow(df)

  if (total / n < 100) {
    n <- floor(total / 100)
  }

  groups <- ceiling(seq_along(1:total) / (total / n))

  return(split(df, groups))
}


