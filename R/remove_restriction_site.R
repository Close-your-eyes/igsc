#' Remove a restriction site using synonymous codon substitutions
#'
#' Modify a DNA sequence to remove occurrences of a restriction recognition
#' sequence while preserving translation under the standard genetic code.
#' An optional protected interval permits sites fully contained within it and
#' prevents changes to any of its bases.
#'
#' @param dna A single, nonempty character string containing the DNA sequence.
#'   Only A, C, G and T are supported. Lowercase letters are converted to
#'   uppercase and whitespace is removed before coordinates are interpreted.
#' @param frame_start A single integer giving the 1-based position of the first
#'   complete codon in the normalized sequence. Must be within \code{dna}.
#'   Defaults to \code{1L}.
#' @param site A single, nonempty character string containing the restriction
#'   recognition sequence, not an enzyme name. The same normalization as for
#'   \code{dna} is applied. Ambiguous bases and cleavage markers are unsupported.
#' @param protected \code{NULL}, or a numeric vector \code{c(start, end)} giving
#'   one protected interval. Coordinates must be integers, 1-based, inclusive,
#'   ordered, and within the normalized sequence. All bases in this interval
#'   remain unchanged. Only sites fully contained in the interval are permitted;
#'   sites crossing either boundary must still be removed.
#' @param both_strands A single logical value. If \code{TRUE} (the default),
#'   search for both \code{site} and its reverse complement in \code{dna}.
#'   If \code{FALSE}, search only for \code{site} as supplied.
#' @param max_states A single positive integer limiting the number of sequence
#'   states with forbidden sites that the search expands. Defaults to
#'   \code{10000L}. Increasing this limit may require more time and memory.
#' @param verbose A single logical value. If \code{TRUE} (the default), report
#'   each final codon replacement using \code{message()}, including its 1-based
#'   start position in the normalized DNA, original and replacement codons, and
#'   amino acid code and name. Report when no replacements are needed.
#'   Intermediate search edits are not reported. Use \code{FALSE} to silence
#'   these messages. The returned DNA string is unaffected.
#'
#' @details
#' The sequence is treated as linear. All complete codons from
#' \code{frame_start} to the end are eligible for synonymous substitutions.
#' Bases before \code{frame_start} and any trailing incomplete codon remain
#' unchanged. Site detection covers the entire sequence, including those bases.
#'
#' Translation is preserved in the specified forward reading frame only.
#' The search does not stop at stop codons; stop codons may be replaced with
#' other stop codons. Alternative genetic codes and special initiation-codon
#' semantics are not supported.
#'
#' Backtracking checks overlapping occurrences and sites newly created by edits.
#' The first valid sequence found is returned; the result is not guaranteed to
#' minimize nucleotide changes or optimize codon usage.
#'
#' An error is raised for invalid input, when no synonymous solution exists
#' under the constraints, or when the search limit is reached. Reaching the
#' limit does not establish that a solution is impossible. No partially edited
#' sequence is returned on error.
#'
#' @return A single uppercase DNA string of the same length as the normalized
#'   input, with no forbidden occurrences of the requested site. If no forbidden
#'   sites are present initially, the normalized input is returned unchanged.
#'
#' @examples
#' # Remove every EcoRI recognition site.
#' remove_restriction_site("ATGGAATTCGAATTC", site = "GAATTC")
#'
#' # Keep the first site and remove the second.
#' remove_restriction_site(
#'   "ATGGAATTCGAATTC", site = "GAATTC", protected = c(4, 9)
#' )
#' # Returns "ATGGAATTCGAATTT".
#'
#' # The first complete codon starts at position 2.
#' remove_restriction_site("CGAATTCGG", frame_start = 2, site = "GAATTC")
#'
#' # Return the edited sequence without replacement messages.
#' remove_restriction_site("ATGGAATTC", site = "GAATTC", verbose = FALSE)
#'
#' @export
remove_restriction_site <- function(dna, frame_start = 1L, site,
                                    protected = NULL,
                                    both_strands = TRUE,
                                    max_states = 10000L,
                                    verbose = TRUE) {
  clean <- function(x, label) {
    if (!is.character(x) || length(x) != 1L || is.na(x))
      stop(label, " must be a single DNA string.")
    x <- toupper(gsub("[[:space:]]", "", x))
    if (!grepl("^[ACGT]+$", x))
      stop(label, " must contain only A, C, G and T.")
    x
  }
  integer_value <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
  }
  dna <- clean(dna, "dna")
  site <- clean(site, "site")
  n <- nchar(dna)
  width <- nchar(site)
  if (!integer_value(frame_start) || frame_start < 1L || frame_start > n)
    stop("frame_start must be a 1-based position within dna.")
  if (!integer_value(max_states) || max_states < 1L)
    stop("max_states must be a positive integer.")
  if (!is.logical(both_strands) || length(both_strands) != 1L ||
      is.na(both_strands)) stop("both_strands must be TRUE or FALSE.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.")
  if (!is.null(protected) &&
      (!is.numeric(protected) || length(protected) != 2L ||
       any(!is.finite(protected)) || any(protected != floor(protected)) ||
       protected[1L] < 1L || protected[2L] > n ||
       protected[1L] > protected[2L]))
    stop("protected must be NULL or c(start, end), within dna.")

  # Codons in T, C, A, G order, with the third base changing fastest.
  bases <- c("T", "C", "A", "G")
  codons <- unlist(lapply(bases, function(a)
    unlist(lapply(bases, function(b) paste0(a, b, bases)))))
  amino_acids <- strsplit(paste0(
    "FFLLSSSSYY**CC*W", "LLLLPPPPHHQQRRRR",
    "IIIMTTTTNNKKSSRR", "VVVVAAAADDEEGGGG"), "")[[1L]]
  code <- setNames(amino_acids, codons)
  synonyms <- split(codons, amino_acids)

  reverse_complement <- function(x)
    paste(rev(strsplit(chartr("ACGT", "TGCA", x), "")[[1L]]), collapse = "")
  motifs <- if (both_strands) unique(c(site, reverse_complement(site))) else site
  windows <- if (width <= n) seq_len(n - width + 1L) else integer()
  allowed <- rep(FALSE, length(windows))
  if (!is.null(protected))
    allowed <- windows >= protected[1L] & windows + width - 1L <= protected[2L]
  forbidden_hits <- function(x) {
    if (!length(windows)) return(integer())
    windows[!allowed & substring(x, windows, windows + width - 1L) %in% motifs]
  }
  starts <- if (frame_start <= n - 2L)
    seq.int(frame_start, n - 2L, by = 3L) else integer()

  # Backtracking explores alternative synonymous edits and checks for new sites.
  pending <- list(dna)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  assign(dna, TRUE, envir = seen)
  visited <- 0L
  while (length(pending)) {
    x <- pending[[length(pending)]]
    pending[[length(pending)]] <- NULL
    hits <- forbidden_hits(x)
    if (!length(hits)) {
      if (verbose) {
        aa_names <- c(
          A = "alanine", R = "arginine", N = "asparagine",
          D = "aspartate", C = "cysteine", Q = "glutamine",
          E = "glutamate", G = "glycine", H = "histidine",
          I = "isoleucine", L = "leucine", K = "lysine",
          M = "methionine", F = "phenylalanine", P = "proline",
          S = "serine", T = "threonine", W = "tryptophan",
          Y = "tyrosine", V = "valine", "*" = "stop"
        )
        changed <- starts[vapply(starts, function(p)
          substr(dna, p, p + 2L) != substr(x, p, p + 2L), logical(1L))]
        if (!length(changed)) message("No codon replacements needed.")
        for (p in changed) {
          original <- substr(dna, p, p + 2L)
          replacement <- substr(x, p, p + 2L)
          aa <- code[[original]]
          message(sprintf("Position %d: %s -> %s (%s, %s)",
                          p, original, replacement, aa, aa_names[[aa]]))
        }
      }
      return(x)
    }
    if (visited >= max_states)
      stop("Search limit reached; increase max_states. No sequence returned.")
    visited <- visited + 1L
    hit <- hits[1L]
    overlapping <- starts[starts <= hit + width - 1L & starts + 2L >= hit]
    for (p in overlapping) {
      old <- substr(x, p, p + 2L)
      for (replacement in setdiff(synonyms[[code[[old]]]], old)) {
        candidate <- x
        substr(candidate, p, p + 2L) <- replacement
        if (!is.null(protected) &&
            substr(candidate, protected[1L], protected[2L]) !=
            substr(dna, protected[1L], protected[2L])) next
        if (!exists(candidate, envir = seen, inherits = FALSE)) {
          assign(candidate, TRUE, envir = seen)
          pending[[length(pending) + 1L]] <- candidate
        }
      }
    }
  }
  stop("No synonymous solution exists with these constraints.")
}
