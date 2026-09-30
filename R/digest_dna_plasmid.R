#' Digest plasmid DNA while checking the sequence boundary
#'
#' Simulate restriction digestion with \code{DECIPHER::DigestDNA}, using
#' rotation to check for cuts missed at the artificial boundary of a circular
#' sequence. Retain the original orientation unless rotation reveals
#' additional cuts on the requested strand or strands.
#'
#' @param sites A nonempty character vector of restriction recognition
#'   sequences and cut locations in the notation accepted by
#'   \code{DECIPHER::DigestDNA}, for example \code{"G/AATTC"} or
#'   \code{"GGTCTC(1/5)"}. Missing values are not allowed.
#' @param seqs A character vector or \code{Biostrings::DNAStringSet} containing
#'   one or more nonempty plasmid sequences in the top strand's 5' to 3'
#'   orientation. Each element is treated as a separate circular molecule.
#'   Sequence names are preserved.
#' @param type Output type: \code{"fragments"} (the default) or
#'   \code{"positions"}.
#' @param strand Strand or strands to digest: \code{"both"} (the default),
#'   \code{"top"}, or \code{"bottom"}. The bottom strand is the reverse
#'   complement of the top strand.
#' @param join_fragments A single logical value. If \code{TRUE}, join the
#'   last and first fragments on each strand to reconstruct the fragment
#'   spanning the artificial sequence boundary. Defaults to \code{FALSE},
#'   which keeps the terminal fragments separate.
#'   Ignored when \code{type = "positions"}.
#' @param reverse_bottom A single logical value, defaulting to \code{TRUE}.
#'   Reverse the nucleotide order within each bottom-strand fragment so it
#'   is displayed in the 3' to 5' direction alongside the top strand.
#'   This reverses the sequence without complementing it. If \code{FALSE},
#'   bottom fragments retain their 5' to 3' sequence orientation.
#'   Bottom-fragment ordering is adjusted in either case. Ignored for
#'   position output and when \code{strand = "top"}.
#'
#' @details
#' Cut positions from rotated sequences are mapped back to the original
#' coordinates for comparison. Cuts are compared by location on each
#' requested strand, rather than by their total count.
#'
#' If rotation reveals no additional cuts, the sequence is passed to
#' \code{DECIPHER::DigestDNA} in its original orientation. Otherwise, the
#' function chooses a rotation that retains every detected cut on the
#' requested strands. Sequences are processed independently, so only those
#' requiring rotation are changed.
#'
#' The boundary check normally uses the original sequence and a rotation
#' by half its length. Short sequences require more rotations. If additional
#' cuts are found but the checked rotations omit other cuts, further
#' rotations are searched until a suitable orientation is found.
#'
#' When \code{type = "fragments"} and \code{join_fragments = TRUE}, the
#' last and first linear fragments from
#' \code{DECIPHER::DigestDNA} are joined, in that order, to reconstruct the
#' fragment spanning the artificial sequence boundary. Joining is performed
#' separately for each sequence and strand, only when that strand has more
#' than one fragment. This applies to both original and rotated sequences.
#' With \code{join_fragments = FALSE}, terminal fragments remain separate.
#'
#' Bottom fragments are ordered along the top strand: their order from
#' \code{DECIPHER::DigestDNA} is reversed. If terminal fragments have been
#' joined, the joined boundary fragment remains first on both strands and
#' only the internal bottom fragments are reordered. The first bottom
#' fragment therefore corresponds to the first top fragment, and so on
#' for digests with corresponding cuts on both strands. Sticky ends can
#' give corresponding fragments different lengths; no padding is added.
#' The \code{reverse_bottom} argument controls nucleotide order within
#' bottom fragments independently of this fragment ordering.
#'
#' Position output is returned directly from \code{DECIPHER::DigestDNA}.
#' Positions refer to the orientation used for the final digest; they are
#' not mapped back when a sequence has been rotated.
#'
#' @return
#' For \code{type = "fragments"}, an ordinary R list with one
#' \code{Biostrings::DNAStringSet} per input sequence, preserving input
#' order and sequence names. Within each set, all \code{"top"} fragments
#' precede all \code{"bottom"} fragments, according to \code{strand}.
#' The bottom fragments follow the corresponding top-fragment order.
#' When joining is
#' enabled, the joined boundary fragment is placed first within each strand,
#' followed by any internal fragments in their original order. A strand cut
#' once then yields one full-length linear fragment beginning immediately
#' after the cut on that strand, before any bottom-sequence reversal. An
#' uncut strand is returned once as its full sequence, with bottom-strand
#' nucleotide order controlled by \code{reverse_bottom}.
#'
#' For \code{type = "positions"}, a list with one element per input sequence.
#' Each element is a named list containing numeric cut-position vectors
#' for \code{top} and/or \code{bottom}. Positions identify the nucleotide
#' immediately after the cut, counted from the respective strand's 5' end.
#' A strand with no cuts has an empty vector.
#'
#' @section Messages and errors:
#' A message identifies each rotated sequence and its leftward shift in
#' bases. For position output, the message also states that coordinates
#' refer to the rotated sequence.
#'
#' When fragment joining is enabled and a strand has exactly one cut, a
#' message identifies the sequence and affected strand or strands. The
#' joined fragment has the full sequence length even though a cut occurred,
#' so the result may otherwise appear unchanged. No such message is emitted
#' for uncut strands, when joining is disabled, or for position output.
#' Messages can be silenced with \code{suppressMessages()}.
#'
#' Each plasmid must be longer than the region spanning every supplied
#' recognition sequence and its associated cuts. The function stops with
#' an error if this requirement is not met, or if no single rotation retains
#' every detected cut on the requested strands.
#'
#' @seealso \code{\link[DECIPHER]{DigestDNA}},
#'   \code{\link[Biostrings]{DNAStringSet}}
#'
#' @examples
#' # An internal EcoRI site: retain the original orientation.
#' tt <- digest_dna_plasmid("G/AATTC", c(plasmid = "AAAAAGAATTCAAAAAAAAA"))
#'
#' # print the cut sites
#' dna_double_strand_print(tt[1], tt[3])
#' dna_double_strand_print(tt[2], tt[4], align = "right")
#'
#' # Join the ends: one cut yields a full-length fragment and a message.
#' digest_dna_plasmid(
#'   "G/AATTC", c(plasmid = "AAAAAGAATTCAAAAAAAAA"), join_fragments = TRUE
#' )
#'
#' # An EcoRI site spanning the boundary: rotate and report the shift.
#' seqs <- c(plasmid = "AATTCAAAAAAAAAAAAAAG")
#' digest_dna_plasmid("G/AATTC", seqs, type = "fragments")
#'
#' # get EcoRI site from DECIPHER package
#' library(DECIPHER)
#' data(RESTRICTION_ENZYMES)
#' ecori <- RESTRICTION_ENZYMES["EcoRI"]
#' names(RESTRICTION_ENZYMES)
#'
#' # DNAStringSet input and top-strand fragments.
#' digest_dna_plasmid(
#'   "G/AATTC", Biostrings::DNAStringSet(seqs), strand = "top"
#' )
#'
#' # Keep bottom sequences in 5' to 3' orientation, but align fragment order.
#' fragments <- digest_dna_plasmid("G/AATTC", seqs, reverse_bottom = FALSE)
#' fragments[[1L]]
#'
#' @export
digest_dna_plasmid <- function(sites,
                               seqs,
                               type = c("fragments", "positions"),
                               strand = c("both", "top", "bottom"),
                               join_fragments = FALSE,
                               reverse_bottom = TRUE) {

  type <- rlang::arg_match(type)
  strand <- rlang::arg_match(strand)
  if (!is.logical(join_fragments) || length(join_fragments) != 1L ||
      is.na(join_fragments)) {
    stop("join_fragments must be a single TRUE or FALSE.")
  }
  if (!is.logical(reverse_bottom) || length(reverse_bottom) != 1L ||
      is.na(reverse_bottom)) {
    stop("reverse_bottom must be a single TRUE or FALSE.")
  }

  if (is.character(seqs)) seqs <- Biostrings::DNAStringSet(seqs)
  if (!methods::is(seqs, "DNAStringSet")) {
    stop("seqs must be a character vector or DNAStringSet.")
  }
  if (!length(seqs) || any(BiocGenerics::width(seqs) == 0L)) {
    stop("seqs must contain nonempty sequences.")
  }
  if (!is.character(sites) || !length(sites) || anyNA(sites)) {
    stop("sites must be a nonempty character vector without NA.")
  }

  # This also lets DECIPHER validate the restriction-site notation.
  cuts <- DECIPHER::DigestDNA(sites, seqs, "positions", strand)

  # Bound the region needed to recognize a site and make its cuts.
  spans <- vapply(sites, function(site) {
    if (grepl("(", site, fixed = TRUE)) {
      motif <- sub("\\(.*$", "", site)
      offsets <- as.numeric(strsplit(
        sub(".*\\((.*)\\)$", "\\1", site), "/", fixed = TRUE
      )[[1L]])
      boundaries <- c(0, nchar(motif), nchar(motif) + offsets)
    } else {
      boundaries <- c(0, nchar(site) - 1L)
    }
    diff(range(boundaries))
  }, numeric(1))
  span <- max(spans)
  if (any(BiocGenerics::width(seqs) <= span)) {
    stop("Each plasmid must be longer than every recognition-and-cut span.")
  }

  rotate <- function(x, shift) {
    if (shift == 0L) return(x)
    paste0(substring(x, shift + 1L), substr(x, 1L, shift))
  }

  adjusted <- seqs
  notices <- character()
  for (i in seq_along(seqs)) {
    x <- as.character(seqs[[i]])
    n <- nchar(x)
    half <- n %/% 2L
    original <- cuts[[i]]
    all_cuts <- original
    mapped_cuts <- function(shift) {
      rotated <- DECIPHER::DigestDNA(
        sites, rotate(x, shift), "positions", strand
      )[[1L]]
      for (s in names(rotated)) {
        # The reverse complement rotates in the opposite direction.
        offset <- if (s == "top") shift else -shift
        rotated[[s]] <- sort(unique(
          (rotated[[s]] - 1L + offset) %% n + 1L
        ))
      }
      rotated
    }
    same_cuts <- function(a, b) {
      all(vapply(names(a), function(s) setequal(a[[s]], b[[s]]), logical(1)))
    }

    # Two orientations suffice for normal plasmids. Very short circles
    # need more rotations to avoid hiding a site or cut in both views.
    shifts <- if (half > span) half else seq_len(n - 1L)
    checked <- lapply(shifts, mapped_cuts)
    for (candidate in checked) {
      for (s in names(candidate)) {
        all_cuts[[s]] <- union(all_cuts[[s]], candidate[[s]])
      }
    }
    if (same_cuts(original, all_cuts)) next

    # Use an orientation that retains every known cut, including cuts
    # visible only across the original sequence boundary.
    complete <- vapply(checked, same_cuts, logical(1), b = all_cuts)
    chosen <- if (any(complete)) shifts[which(complete)[1L]] else NA_integer_
    if (is.na(chosen)) {
      for (shift in setdiff(seq_len(n - 1L), shifts)) {
        if (same_cuts(mapped_cuts(shift), all_cuts)) {
          chosen <- shift
          break
        }
      }
    }
    if (is.na(chosen)) {
      stop("Sequence ", i, ": no rotation retains every cut. ",
           "DigestDNA cannot represent this circular digest in one linear orientation.")
    }

    adjusted[i] <- Biostrings::DNAStringSet(rotate(x, chosen))
    label <- if (is.null(names(seqs))) as.character(i) else
      paste0(i, " ('", names(seqs)[i], "')")
    notices <- c(notices, paste0(
      "Sequence ", label, " was rotated left by ", chosen,
      " bases to include cuts spanning the original boundary.",
      if (type == "positions") " Positions refer to the rotated sequence." else ""
    ))
  }
  names(adjusted) <- names(seqs)
  result <- DECIPHER::DigestDNA(sites, adjusted, type, strand)
  if (type == "fragments" && join_fragments) {
    joined <- vector("list", length(result))
    names(joined) <- names(result)
    for (i in seq_along(result)) {
      fragments <- result[[i]]
      strand_names <- unique(names(fragments))
      single_cut <- strand_names[vapply(strand_names, function(s) {
        sum(names(fragments) == s) == 2L
      }, logical(1))]
      if (length(single_cut)) {
        label <- if (is.null(names(seqs))) as.character(i) else
          paste0(i, " ('", names(seqs)[i], "')")
        notices <- c(notices, paste0(
          "Sequence ", label, ": exactly one cut on each listed strand (",
          paste(single_cut, collapse = ", "), "). Joining the end fragments ",
          "returns one full-length fragment per listed strand. The strands ",
          "were cut, even though their lengths are unchanged."
        ))
      }
      strands <- lapply(strand_names, function(s) {
        pieces <- fragments[names(fragments) == s]
        if (length(pieces) > 1L) {
          last <- length(pieces)
          boundary <- Biostrings::DNAStringSet(paste0(
            as.character(pieces[c(last, 1L)]), collapse = ""
          ))
          pieces <- c(boundary, pieces[-c(1L, last)])
          names(pieces) <- rep(s, length(pieces))
        }
        pieces
      })
      joined[[i]] <- do.call(c, strands)
    }
    result <- joined
  }
  if (type == "fragments") {
    result <- lapply(result, function(fragments) {
      top <- fragments[names(fragments) == "top"]
      bottom <- fragments[names(fragments) == "bottom"]
      if (join_fragments && length(bottom) > 1L) {
        # Both joined boundary fragments must remain first in their groups.
        bottom <- c(bottom[1L], rev(bottom[-1L]))
      } else {
        bottom <- rev(bottom)
      }
      if (reverse_bottom) bottom <- Biostrings::reverse(bottom)
      names(top) <- make.unique(names(top))
      names(bottom) <- make.unique(names(bottom))
      c(top, bottom)
    })
  }
  for (notice in notices) message(notice)
  result
}

# An internal EcoRI site: retain the original orientation.
# tt <- digest_dna_plasmid("G/AATTC", c(plasmid = "AAAAAGAATTCAAAAAAAAA"), reverse_bottom = T)[[1]]
# tt
# dna_double_strand_print(tt[1],tt[3])
# dna_double_strand_print(tt[2],tt[4], align = "right")
#
# tt <- digest_dna_plasmid("G/AATTC", c(plasmid = "AAAAAGAATTCAAAAAAAAA"), reverse_bottom = F)[[1]]
# tt
# dna_double_strand_print(tt[1],tt[3])


