#' Align sequences within clusters and combine their consensuses
#'
#' @description
#' Optionally align the input sequences, cluster them using DECIPHER sequence
#' distances, and independently align each cluster. Construct one consensus per
#' cluster, align these consensuses, and construct their overall consensus.
#' Return a distance matrix and heatmap of the aligned cluster consensuses.
#' Optionally create alignment plots using caller-provided helper functions.
#'
#' @param seqs A nonempty character vector or an object inheriting from
#'   \code{XStringSet}. Sequences must not be missing, empty, or contain
#'   whitespace. Names must either be absent or be nonempty and unique.
#'   Missing names are generated as \code{seq_1}, \code{seq_2}, and so on.
#' @param treeline_args A plain list with unique, nonempty argument names,
#'   forwarded to \code{DECIPHER::Treeline()} (\code{TreeLine()} in older
#'   releases). Defaults are \code{method = "complete"} and
#'   \code{cutoff = 0.01}; omitted defaults are retained. The cutoff must be
#'   one finite, nonnegative number. The arguments \code{myXStringSet},
#'   \code{myDistMatrix}, and \code{type} are controlled by this function
#'   and cannot be supplied here.
#' @param alignseqs_args A plain list with unique, nonempty argument names,
#'   forwarded to \code{DECIPHER::AlignSeqs()} for every alignment. Names
#'   must match formal arguments of \code{AlignSeqs()} or
#'   \code{DECIPHER::AlignProfiles()}, excluding \code{pattern} and
#'   \code{subject}. The argument \code{myXStringSet} cannot be overridden.
#'   The arguments \code{guideTree} and \code{structures} must be absent or
#'   \code{NULL}, because the sequence inputs differ between alignments.
#'   Do not pass \code{p.weight} or \code{s.weight}: \code{AlignSeqs()}
#'   supplies them internally, so passing them causes a duplicate-argument
#'   error. Other profile arguments are checked by DECIPHER when an
#'   alignment is run.
#' @param init_with_alignment A single nonmissing logical value. If
#'   \code{TRUE}, align the input before calculating distances. If
#'   \code{FALSE}, treat the input as an existing alignment; all input
#'   sequences must have equal widths. Cluster members are independently
#'   realigned in either case.
#' @param consensus_thresh A single finite numeric value from zero up to,
#'   but excluding, one. Passed as \code{threshold} to
#'   \code{DECIPHER::ConsensusSequence()} for each cluster. This controls
#'   the fraction of sequence information that may be discarded at a
#'   position; it is not a minimum agreement proportion.
#' @param consensus_consensus_thresh A threshold with the same requirements
#'   and interpretation as \code{consensus_thresh}, used to construct the
#'   overall consensus from the aligned cluster consensuses.
#' @param seq_type \code{NULL}, \code{"DNA"}, \code{"RNA"}, or \code{"AA"},
#'   ignoring case. With \code{NULL}, typed DNA, RNA, and amino-acid sets
#'   retain their biological type; character vectors and other
#'   \code{XStringSet} subclasses default to DNA. For other biological
#'   types, specify this argument explicitly. An explicit type must agree
#'   with an input \code{DNAStringSet}, \code{RNAStringSet}, or
#'   \code{AAStringSet}.
#' @param make_plots A single nonmissing logical value. If \code{TRUE},
#'   \code{compare_seq_df_wide()} and \code{aln_plot()} must be available
#'   in the caller's environment or an enclosing environment. If
#'   \code{FALSE}, plotting helpers are unnecessary and both plot return
#'   components are \code{NULL}. The consensus distance heatmap is still
#'   computed when this is \code{FALSE}.
#'
#' @details
#' Inputs are converted to uppercase and validated against the selected
#' biological alphabet. Generic sets, including \code{BStringSet}, are
#' accepted through conversion to a DNA, RNA, or amino-acid set; arbitrary
#' text alphabets are not supported by the DECIPHER alignment functions.
#'
#' Gap characters \code{"-"} and \code{"."} are removed before each
#' alignment. Sequences containing only these characters are rejected.
#' Existing gaps are retained for the distance calculation when
#' \code{init_with_alignment = FALSE}. Unknown argument names and reserved
#' argument overrides are rejected; errors from downstream operations are
#' reported with the operation's context.
#'
#' Cluster results are ordered by decreasing member count, and members
#' within each cluster are ordered by decreasing ungapped sequence length.
#' Cluster consensuses are ordered by decreasing width. Their names identify
#' their clusters as \code{cluster_<id>}; cluster assignments retain input
#' order. Consensus positions without sufficient information use
#' \code{"N"} for DNA/RNA and \code{"X"} for amino acids.
#'
#' Single-sequence alignments bypass \code{AlignSeqs()}. A single input
#' sequence is assigned to cluster 1 and has no dendrogram. Each cluster
#' contributes one sequence to the final consensus alignment, so cluster
#' size does not directly weight the overall consensus. DECIPHER checks
#' some forwarded argument values only when the corresponding operation
#' runs; those values may go unchecked for a single input sequence.
#'
#' Alignment plots convert each \code{XStringSet} using
#' \code{xstringset_to_df()}, then pivot it to wide form. The first sequence
#' name is the reference passed to \code{compare_seq_df_wide()}, along with
#' \code{change_nonref = TRUE}, \code{nonref_mismatch_as = "mismatch_symbol"},
#' and \code{return_as_long = TRUE}. The result is passed to
#' \code{aln_plot()}; each cluster plot receives a title with its cluster ID
#' and member count. Plotting also uses \code{tidyr} and \code{ggplot2}.
#'
#' The consensus distance heatmap is always calculated using
#' \code{DECIPHER::DistanceMatrix()}, \code{brathering::mat_to_df_long()},
#' and \code{fcexpr::heatmap_long_df()}. These packages must be available
#' even when \code{make_plots = FALSE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{clusters}{A named vector of cluster assignments in input order.}
#'   \item{dendrogram}{The DECIPHER dendrogram, or \code{NULL} for a single
#'     input sequence.}
#'   \item{clustalns}{A named list of aligned biological sequence sets, with
#'     cluster IDs as list names and original sequence names retained.}
#'   \item{clustaln_plots}{A list of helper-generated plot objects with the
#'     same names and order as \code{clustalns}, or \code{NULL}.}
#'   \item{clustconsensuses}{A biological sequence set containing one
#'     consensus per cluster, ordered by decreasing width. These sequences
#'     have not yet been aligned to one another and may contain gaps.}
#'   \item{consensusaln}{The aligned cluster consensuses.}
#'   \item{consensusaln_plot}{The helper-generated consensus alignment plot,
#'     or \code{NULL}.}
#'   \item{consensusconsensus}{A biological sequence set of length one
#'     containing the overall consensus.}
#'   \item{consensus_distmat}{A distance matrix calculated from the aligned
#'     cluster consensuses.}
#'   \item{consensus_distmat_htmp}{A heatmap of \code{consensus_distmat},
#'     with cluster consensus names on both axes.}
#' }
#' All returned biological sequence sets use the selected DNA, RNA, or
#' amino-acid type.
#'
#' @examples
#' if (requireNamespace("Biostrings", quietly = TRUE) &&
#'     requireNamespace("DECIPHER", quietly = TRUE) &&
#'     requireNamespace("brathering", quietly = TRUE) &&
#'     requireNamespace("fcexpr", quietly = TRUE)) {
#'   dna <- c(first = "ACGTACGTACGT", second = "ACGTACGTACGA")
#'   result <- align_sequences_clusterwise(dna, make_plots = FALSE)
#'   result$clusters
#'   result$consensusconsensus
#'   result$consensus_distmat
#'
#'   protein <- Biostrings::BStringSet(c(
#'     first = "MKTLLAVAVAAA", second = "MKTLLAVAVAAT"
#'   ))
#'   protein_result <- align_sequences_clusterwise(
#'     protein, seq_type = "AA", make_plots = FALSE
#'   )
#' }
#'
#' @export
align_sequences_clusterwise <- function(
    seqs,
    treeline_args = list(method = "complete", cutoff = 0.01),
    alignseqs_args = list(),
    init_with_alignment = TRUE,
    consensus_thresh = 0.3,
    consensus_consensus_thresh = 0.7,
    seq_type = NULL,
    make_plots = TRUE) {

  fail <- function(...) stop(..., call. = FALSE)
  run <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      fail(label, ": ", conditionMessage(e))
    })
  }
  check_flag <- function(x, label) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      fail(label, " must be a single TRUE or FALSE.")
    }
  }
  check_threshold <- function(x, label) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || x >= 1) {
      fail(label, " must be a finite number from zero up to, but excluding, one.")
    }
  }
  check_args <- function(x, label, allowed, reserved) {
    if (!is.list(x) || is.object(x)) {
      fail(label, " must be a plain named list.")
    }
    if (!length(x)) return(invisible(NULL))
    nm <- names(x)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
      fail(label, " must have nonempty, unique argument names.")
    }
    bad <- intersect(nm, reserved)
    if (length(bad)) {
      fail(label, " cannot override: ", paste(bad, collapse = ", "), ".")
    }
    bad <- setdiff(nm, allowed)
    if (length(bad)) {
      fail("Unknown argument(s) in ", label, ": ", paste(bad, collapse = ", "), ".")
    }
  }

  check_flag(init_with_alignment, "init_with_alignment")
  check_flag(make_plots, "make_plots")
  check_threshold(consensus_thresh, "consensus_thresh")
  check_threshold(consensus_consensus_thresh, "consensus_consensus_thresh")

  igsc:::.ensure_packages(c("Biostrings", "DECIPHER", "fcexpr", "brathering"))
  # for (pkg in c("Biostrings", "DECIPHER")) {
  #   if (!requireNamespace(pkg, quietly = TRUE)) fail("Install package '", pkg, "'.")
  # }

  if (missing(seqs)) fail("seqs is required.")
  if (!(is.character(seqs) && is.null(dim(seqs))) &&
      !methods::is(seqs, "XStringSet")) {
    fail("seqs must be a character vector or an XStringSet.")
  }
  if (!length(seqs)) fail("seqs must contain at least one sequence.")

  if (!is.null(seq_type)) {
    if (!is.character(seq_type) || length(seq_type) != 1L || is.na(seq_type)) {
      fail("seq_type must be NULL, 'DNA', 'RNA', or 'AA'.")
    }
    seq_type <- toupper(seq_type)
    if (!seq_type %in% c("DNA", "RNA", "AA")) {
      fail("seq_type must be NULL, 'DNA', 'RNA', or 'AA'.")
    }
  }
  input_type <- c("DNA", "RNA", "AA")[
    vapply(c("DNAStringSet", "RNAStringSet", "AAStringSet"),
           function(cls) methods::is(seqs, cls), logical(1))
  ]
  if (length(input_type)) {
    if (!is.null(seq_type) && seq_type != input_type) {
      fail("seq_type conflicts with the input's ", input_type, "StringSet class.")
    }
    seq_type <- input_type
  } else if (is.null(seq_type)) {
    seq_type <- "DNA"
  }
  constructor <- getExportedValue("Biostrings", paste0(seq_type, "StringSet"))
  values <- as.character(seqs)
  if (anyNA(values) || any(!nzchar(values))) {
    fail("seqs cannot contain missing or empty sequences.")
  }
  if (any(grepl("[[:space:]]", values))) {
    fail("seqs cannot contain whitespace.")
  }
  seq_names <- names(seqs)
  seqs <- run(paste0("Invalid ", seq_type, " sequence alphabet"),
              constructor(toupper(values)))
  if (is.null(seq_names)) seq_names <- paste0("seq_", seq_along(seqs))
  if (anyNA(seq_names) || any(!nzchar(seq_names)) || anyDuplicated(seq_names)) {
    fail("Sequence names must be nonempty and unique, or entirely absent.")
  }
  names(seqs) <- seq_names

  # Treeline was named TreeLine in older DECIPHER releases.
  exports <- getNamespaceExports("DECIPHER")
  tree_name <- intersect(c("Treeline", "TreeLine"), exports)
  if (!length(tree_name)) fail("This DECIPHER version lacks Treeline/TreeLine.")
  tree_fun <- getExportedValue("DECIPHER", tree_name[1L])
  check_args(treeline_args, "treeline_args", names(formals(tree_fun)),
             c("myXStringSet", "myDistMatrix", "type"))
  check_args(alignseqs_args, "alignseqs_args",
             setdiff(union(names(formals(DECIPHER::AlignSeqs)),
                           names(formals(DECIPHER::AlignProfiles))),
                     c("...", "pattern", "subject")),
             "myXStringSet")
  # These sequence-specific objects cannot be reused for every cluster.
  for (arg in c("guideTree", "structures")) {
    if (!is.null(alignseqs_args[[arg]])) {
      fail("alignseqs_args$", arg, " must be NULL because alignment inputs change.")
    }
  }
  tree_args <- utils::modifyList(list(method = "complete", cutoff = 0.01),
                                 treeline_args, keep.null = TRUE)
  cutoff <- tree_args$cutoff
  if (!is.numeric(cutoff) || is.complex(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0) {
    fail("treeline_args$cutoff must be one finite, nonnegative number.")
  }

  if (make_plots) {
    caller <- parent.frame()
    helpers <- c("compare_seq_df_wide", "aln_plot")
    plot_helpers <- lapply(helpers, function(nm) {
      get0(nm, envir = caller, mode = "function", inherits = TRUE)
    })
    names(plot_helpers) <- helpers
    missing_helpers <- helpers[vapply(plot_helpers, is.null, logical(1))]
    if (length(missing_helpers)) {
      fail("Missing plotting helper(s): ", paste(missing_helpers, collapse = ", "),
           ". Define them or set make_plots = FALSE.")
    }
  }


  plot_alignment <- function(x, label, id = NULL) {
    # Local column IDs avoid assumptions about the caller's sequence names.
    # wide <- as.data.frame(as.list(unname(as.character(x))),
    #                       stringsAsFactors = FALSE)
    # names(wide) <- paste0("seq_", seq_along(x))
    # names(x) <- paste0("seq_", seq_along(x))

    run(label, {
      long <- plot_helpers$compare_seq_df_wide(
        tidyr::pivot_wider(xstringset_to_df(x), names_from = seq.name, values_from = seq),
        change_nonref = T,
        ref = names(x)[1],
        nonref_mismatch_as = "mismatch_symbol",
        return_as_long = T
      )
      pp <- plot_helpers$aln_plot(long)
      if (!is.null(id)) {
        pp <- pp + ggplot2::labs(title = paste0("cluster ", id, ", n = ", length(x)))
      }
      return(pp)
    })
  }

  ungap <- function(x, label) {
    value <- gsub("[-.]", "", as.character(x))
    if (any(!nzchar(value))) fail(label, " contains a sequence with only gaps.")
    result <- constructor(value)
    names(result) <- names(x)
    result
  }
  align <- function(x, label) {
    x <- ungap(x, label)
    if (length(x) == 1L) return(x)
    run(label, do.call(DECIPHER::AlignSeqs,
                       c(list(myXStringSet = x), alignseqs_args)))
  }
  consensus <- function(x, threshold, label) {
    run(label, DECIPHER::ConsensusSequence(
      x, threshold = threshold,
      noConsensusChar = if (seq_type == "AA") "X" else "N"
    ))
  }

  # Validate biological content even when the initial alignment is skipped.
  raw_seqs <- ungap(seqs, "seqs")
  if (!init_with_alignment && length(unique(Biostrings::width(seqs))) != 1L) {
    fail("With init_with_alignment = FALSE, all sequences must have equal widths.")
  }
  seqsaln <- if (init_with_alignment) align(raw_seqs, "Initial alignment") else seqs
  if (length(seqs) == 1L) {
    clusters <- stats::setNames(1L, seq_names)
    dend <- NULL
  } else {
    distmat <- run("Distance calculation", DECIPHER::DistanceMatrix(seqsaln))
    if (any(!is.finite(distmat))) {
      fail("Distance matrix contains nonfinite values; check sequence overlap.")
    }
    cluster_table <- run("Clustering", do.call(
      tree_fun, c(list(myDistMatrix = distmat, type = "clusters"), tree_args)
    ))
    if (NROW(cluster_table) != length(seqs) || NCOL(cluster_table) != 1L) {
      fail("Treeline returned an unexpected cluster table.")
    }
    clusters <- cluster_table[, 1L]
    if (anyNA(clusters)) fail("Treeline returned missing cluster assignments.")
    names(clusters) <- seq_names
    dend <- run("Dendrogram construction", do.call(
      tree_fun, c(list(myDistMatrix = distmat, type = "dendrogram"), tree_args)
    ))
  }

  # Keep assignments in input order; sort only the per-cluster results.
  groups <- split(seq_along(raw_seqs), clusters)
  groups <- groups[order(lengths(groups), decreasing = TRUE)]
  message(length(groups), " clusters.")
  clustalns <- lapply(names(groups), function(id) {
    x <- raw_seqs[groups[[id]]]
    x <- x[order(Biostrings::width(x), decreasing = TRUE)]
    align(x, paste0("Alignment of cluster ", id))
  })
  names(clustalns) <- names(groups)

  clustaln_plots <- if (make_plots) {
    lapply(names(clustalns), function(id) {
      plot_alignment(clustalns[[id]], paste0("Plot of cluster ", id), id)
    })
  } else NULL
  if (make_plots) names(clustaln_plots) <- names(clustalns)



  cons <- lapply(names(clustalns), function(id) {
    consensus(clustalns[[id]], consensus_thresh, paste0("Consensus of cluster ", id))
  })
  consensuses <- constructor(vapply(cons, function(x) as.character(x)[1L], character(1)))
  names(consensuses) <- paste0("cluster_", names(clustalns))
  consensuses <- consensuses[order(Biostrings::width(consensuses), decreasing = TRUE)]
  consensaln <- align(consensuses, "Alignment of cluster consensuses")

  consensus_distmat <- NULL
  consensus_distmat_htmp <- NULL
  if (length(consensaln)>1) {
    consensus_distmat <- DECIPHER::DistanceMatrix(consensaln)
    consensus_distmat_htmp <- fcexpr::heatmap_long_df(df = brathering::mat_to_df_long(consensus_distmat,
                                                                                      rownames_to = "cl1",
                                                                                      colnames_to = "cl2",
                                                                                      values_to = "distance"),
                                                      groups = "cl1",
                                                      features = "cl2",
                                                      values = "distance",
                                                      colorsteps = 7)
  }

  consensaln_plot <- if (make_plots) {
    plot_alignment(consensaln, "Plot of cluster consensuses")
  } else NULL
  consensconsens <- consensus(consensaln, consensus_consensus_thresh,
                              "Consensus of cluster consensuses")

  list(
    clusters = clusters,
    dendrogram = dend,
    clustalns = clustalns,
    clustaln_plots = clustaln_plots,
    clustconsensuses = consensuses,
    consensusaln = consensaln,
    consensusaln_plot = consensaln_plot,
    consensusconsensus = consensconsens,
    consensus_distmat = consensus_distmat,
    consensus_distmat_htmp = consensus_distmat_htmp)
}


# align_sequences_clusterwise <- function(seqs,
#                                         treeline_args = list(method = "complete",
#                                                              cutoff = 0.01),
#                                         alignseqs_args = list(),
#                                         init_with_alignment = T,
#                                         consensus_thresh = 0.3,
#                                         consensus_consensus_thresh = 0.7) {
#
#   # seqs: character or XStringSet
#
#   # init with alignment
#   if (init_with_alignment) {
#     seqsaln <- do.call(DECIPHER::AlignSeqs, args = c(list(myXStringSet = seqs), alignseqs_args))
#   } else {
#     seqsaln <- seqs
#   }
#
#   distmat <- DECIPHER::DistanceMatrix(seqsaln)
#
#   clusters <- do.call(DECIPHER::Treeline, args = c(list(myDistMatrix = distmat,
#                                                         type = "clusters"),
#                                                    treeline_args))[,1]
#
#   dend <- do.call(DECIPHER::Treeline, args = c(list(myDistMatrix = distmat),
#                                                treeline_args))
#   # plot(dend)
#
#   clust <- table(clusters)
#   message(length(clust), " clusters.")
#
#   seqs2 <- split(seqs, clusters)
#   # order sequences by length
#   seqs2 <- purrr::map(seqs2, ~.x[order(purrr::map_dbl(as.character(.x), nchar), decreasing = T)])
#   # order groups by number of members
#   seqs2 <- seqs2[order(purrr::map_dbl(seqs2, ~length(.x)), decreasing = T)]
#   clustalns <- purrr::map(seqs2, function(x) {
#     if (length(x) == 1) {
#       return(x)
#     }
#     do.call(DECIPHER::AlignSeqs, args = c(list(myXStringSet = x), alignseqs_args))
#   })
#
#   clustaln_dfs <- purrr::map(clustalns, ~compare_seq_df_wide(
#     tidyr::pivot_wider(xstringset_to_df(.x), names_from = seq.name, values_from = seq),
#     change_nonref = T,
#     ref = "seq_1",
#     nonref_mismatch_as = "mismatch_symbol",
#     return_as_long = T
#   )
#   )
#
#   clustaln_plots <- purrr::map(clustaln_dfs, ~aln_plot(.x))
#
#   consensuses <- purrr::map_chr(clustalns, ~as.character(DECIPHER::ConsensusSequence(.x, threshold = consensus_thresh)))
#   consensuses <- consensuses[order(purrr::map_dbl(consensuses, nchar), decreasing = T)]
#   consensuses <- Biostrings::DNAStringSet(consensuses)
#
#   consensaln <- do.call(DECIPHER::AlignSeqs, args = c(list(myXStringSet = consensuses), alignseqs_args))
#   consensaln_df <- compare_seq_df_wide(
#     tidyr::pivot_wider(xstringset_to_df(consensaln), names_from = seq.name, values_from = seq),
#     change_nonref = T,
#     ref = names(consensaln)[1],
#     nonref_mismatch_as = "mismatch_symbol",
#     return_as_long = T
#   )
#
#   consensaln_plot <- aln_plot(consensaln_df)
#
#   # DECIPHER::BrowseSeqs(consensaln)
#   consensconsens <- DECIPHER::ConsensusSequence(consensaln, threshold = consensus_consensus_thresh)
#
#   return(list(clusters = clusters,
#               dendrogram = dend,
#               clustalns = clustalns,
#               clustaln_plots = clustaln_plots,
#               clustconsensuses = consensuses,
#               consensusaln = consensaln,
#               consensusaln_plot = consensaln_plot,
#               consensusconsensus = consensconsens))
# }
