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
#'   \code{"positions"}. Unambiguous abbreviations are accepted.
#' @param strand Strand or strands to digest: \code{"both"} (the default),
#'   \code{"top"}, or \code{"bottom"}. The bottom strand is the reverse
#'   complement of the top strand. Unambiguous abbreviations are accepted.
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
#' The final result is returned directly from \code{DECIPHER::DigestDNA}.
#' Its linear digestion behavior is preserved: terminal fragments are not
#' joined across the artificial boundary. Returned positions refer to the
#' orientation used for the final digest; they are not mapped back when a
#' sequence has been rotated.
#'
#' @return
#' For \code{type = "fragments"}, a \code{DNAStringSetList} with one element
#' per input sequence. Each element contains fragments named \code{"top"}
#' and/or \code{"bottom"}, according to \code{strand}. An uncut strand is
#' returned as its full sequence in the orientation used for digestion.
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
#' refer to the rotated sequence. Messages can be silenced with
#' \code{suppressMessages()}.
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
#' digest_dna_plasmid("G/AATTC", c(plasmid = "AAAAAGAATTCAAAAAAAAA"))
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
#' @export
digest_dna_plasmid <- function(sites,
                               seqs,
                               type = c("fragments", "positions"),
                               strand = c("both", "top", "bottom")) {

  type <- rlang::arg_match(type)
  strand <- rlang::arg_match(strand)

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
  for (notice in notices) message(notice)
  result
}
