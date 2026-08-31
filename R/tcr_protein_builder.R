# Reconstruct a full human TCR protein chain from V/(D)/J calls and an
# amino-acid CDR3, using an IMGT-style reference data frame.
#
# Load the supplied reference dump with, for example:
#   imgt_ref <- dget("Pasted text(20260827-134039).txt")
#
# Then build a beta chain:
#   tcr <- build_tcr_protein(
#     ref = imgt_ref,
#     v_gene = "TRBV7-9*01",
#     d_gene = "TRBD1*01",
#     j_gene = "TRBJ2-7*01",
#     cdr3 = "CASSLGQAYEQYF"
#   )
#   tcr$sequence
#   cat(as_fasta(tcr))
#
# The D allele is validated and reported, but it is not appended separately:
# its protein contribution is already present in the CDR3 and is not uniquely
# identifiable from an amino-acid CDR3 alone.

.tcr_codon_table <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

.tcr_translate <- function(dna, frame = 0L) {
  dna <- toupper(gsub("\\s+", "", dna))
  frame <- as.integer(frame)
  if (length(frame) != 1L || is.na(frame) || !frame %in% 0:2) {
    stop("frame must be 0, 1, or 2", call. = FALSE)
  }
  if (nchar(dna) - frame < 3L) return("")
  starts <- seq.int(frame + 1L, nchar(dna) - 2L, by = 3L)
  if (!length(starts)) return("")
  codons <- substring(dna, starts, starts + 2L)
  amino_acids <- unname(.tcr_codon_table[codons])
  amino_acids[is.na(amino_acids)] <- "X"
  paste0(amino_acids, collapse = "")
}

.tcr_check_reference <- function(ref) {
  needed <- c("seq.nt", "meta", "Allele", "Gene", "Functionality",
              "seq.aa", "CDR3")
  missing <- setdiff(needed, names(ref))
  if (length(missing)) {
    stop("Reference is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!nrow(ref)) stop("Reference has no rows", call. = FALSE)
  invisible(TRUE)
}

.tcr_clean_gene <- function(x, argument) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(argument, " must be one non-empty gene or allele name",
         call. = FALSE)
  }
  toupper(gsub("\\s+", "", x))
}

.tcr_is_functional <- function(x) {
  x %in% c("F", "[F]", "(F)")
}

.tcr_closest_gene <- function(genes, query) {
  family_match <- regexpr("^(TRAV|TRBV|TRAJ|TRBJ)", query, perl = TRUE)
  if (family_match[[1L]] < 0L) return(NULL)

  family <- regmatches(query, family_match)
  candidates <- sort(unique(genes[
    !is.na(genes) & nzchar(genes) & startsWith(genes, family)
  ]))
  if (!length(candidates)) return(NULL)

  distances <- as.integer(adist(query, candidates))
  closest <- candidates[distances == min(distances)]
  if (length(closest) > 1L) {
    stop("No unique closest ", family, " gene segment for ", query,
         "; equally close matches are: ", paste(closest, collapse = ", "),
         call. = FALSE)
  }
  closest[[1L]]
}

.tcr_resolve_segment <- function(ref, query, segment, strict) {
  query <- .tcr_clean_gene(query, segment)
  alleles <- toupper(as.character(ref$Allele))
  genes <- toupper(as.character(ref$Gene))
  inferred <- FALSE
  used_closest <- FALSE
  requested <- query

  if (grepl("*", query, fixed = TRUE)) {
    index <- which(alleles == query)
  } else {
    candidates <- which(genes == query)
    preferred <- candidates[alleles[candidates] == paste0(query, "*01")]
    if (length(preferred) == 1L) {
      index <- preferred
      inferred <- TRUE
    } else if (length(candidates) == 1L) {
      index <- candidates
      inferred <- TRUE
    } else if (length(candidates) > 1L) {
      stop(segment, " is allele-ambiguous; choose one of: ",
           paste(unique(ref$Allele[candidates]), collapse = ", "),
           call. = FALSE)
    } else {
      index <- integer()
    }
  }

  if (!length(index)) {
    query_gene <- sub("\\*.*$", "", query)
    closest_gene <- .tcr_closest_gene(genes, query_gene)
    if (is.null(closest_gene)) {
      stop(segment, " was not found in the reference: ", query,
           call. = FALSE)
    }

    candidates <- which(genes == closest_gene)
    requested_suffix <- if (grepl("*", query, fixed = TRUE)) {
      sub("^[^*]+", "", query)
    } else {
      ""
    }
    requested_allele <- paste0(closest_gene, requested_suffix)
    same_allele <- if (nzchar(requested_suffix)) {
      candidates[alleles[candidates] == requested_allele]
    } else {
      integer()
    }
    preferred <- candidates[alleles[candidates] == paste0(closest_gene, "*01")]

    if (length(same_allele) == 1L) {
      index <- same_allele
    } else if (length(preferred) == 1L) {
      index <- preferred
    } else if (length(candidates) == 1L) {
      index <- candidates
    } else {
      stop(segment, " matched closest gene segment ", closest_gene,
           " but its allele is ambiguous; choose one of: ",
           paste(unique(ref$Allele[candidates]), collapse = ", "),
           call. = FALSE)
    }

    inferred <- TRUE
    used_closest <- TRUE
    warning(segment, " ", requested, " was not found; using closest match ",
            ref$Allele[index], call. = FALSE)
  }
  if (length(index) != 1L) {
    stop(segment, " resolves to multiple reference rows: ", query,
         call. = FALSE)
  }
  if (strict && !.tcr_is_functional(as.character(ref$Functionality[index]))) {
    stop(segment, " is not annotated as functional: ", ref$Allele[index],
         " (", ref$Functionality[index], ")", call. = FALSE)
  }

  list(index = index, allele = as.character(ref$Allele[index]),
       inferred = inferred, requested = requested,
       used_closest = used_closest)
}

.tcr_chain_from_j <- function(j_allele) {
  if (grepl("^TRAJ", j_allele)) return("alpha")
  if (grepl("^TRBJ", j_allele)) return("beta")
  if (grepl("^TRDJ", j_allele)) return("delta")
  if (grepl("^TRGJ", j_allele)) return("gamma")
  stop("Cannot determine TCR chain from J allele: ", j_allele,
       call. = FALSE)
}

.tcr_validate_v_chain <- function(v_allele, chain) {
  v_gene <- sub("\\*.*$", "", v_allele)
  valid <- switch(
    chain,
    alpha = grepl("^TRAV", v_gene),
    beta = grepl("^TRBV", v_gene),
    gamma = grepl("^TRGV", v_gene),
    delta = grepl("^TRDV", v_gene) || grepl("^TRAV.*DV", v_gene),
    FALSE
  )
  if (!valid) {
    stop("V allele ", v_allele, " is incompatible with a ", chain,
         "-chain J allele", call. = FALSE)
  }
}

.tcr_validate_d_chain <- function(d_allele, chain) {
  if (is.null(d_allele)) return(invisible(TRUE))
  expected <- switch(chain, beta = "^TRBD", delta = "^TRDD", NULL)
  if (is.null(expected)) {
    stop("A ", chain, " chain does not use a D segment", call. = FALSE)
  }
  if (!grepl(expected, d_allele)) {
    stop("D allele ", d_allele, " is incompatible with a ", chain,
         " chain", call. = FALSE)
  }
  invisible(TRUE)
}

.tcr_default_constant <- function(chain, j_allele) {
  switch(
    chain,
    alpha = "TRAC*01",
    beta = if (grepl("^TRBJ1-", j_allele)) "TRBC1*01" else "TRBC2*01",
    delta = "TRDC*01",
    gamma = stop(
      "c_gene is required for gamma chains because the J call does not ",
      "uniquely determine TRGC1 versus TRGC2", call. = FALSE
    )
  )
}

.tcr_resolve_constant <- function(ref, query, chain, strict) {
  query <- .tcr_clean_gene(query, "c_gene")
  alleles <- toupper(as.character(ref$Allele))
  genes <- toupper(as.character(ref$Gene))
  inferred <- FALSE

  if (grepl("*", query, fixed = TRUE)) {
    allele <- query
  } else {
    candidates <- unique(alleles[genes == query])
    preferred <- candidates[candidates == paste0(query, "*01")]
    if (length(preferred) == 1L) {
      allele <- preferred
      inferred <- TRUE
    } else if (length(candidates) == 1L) {
      allele <- candidates
      inferred <- TRUE
    } else if (length(candidates) > 1L) {
      stop("c_gene is allele-ambiguous; choose one of: ",
           paste(candidates, collapse = ", "), call. = FALSE)
    } else {
      allele <- query
    }
  }

  index <- which(alleles == allele)
  if (!length(index)) {
    stop("c_gene was not found in the reference: ", query, call. = FALSE)
  }
  if (strict && any(!.tcr_is_functional(as.character(ref$Functionality[index])))) {
    stop("c_gene is not annotated as functional: ", allele, call. = FALSE)
  }

  expected <- switch(chain, alpha = "^TRAC", beta = "^TRBC",
                     delta = "^TRDC", gamma = "^TRGC")
  if (!grepl(expected, allele)) {
    stop("Constant allele ", allele, " is incompatible with a ", chain,
         " chain", call. = FALSE)
  }

  headers <- as.character(ref$meta[index])
  features <- vapply(strsplit(headers, "|", fixed = TRUE), function(parts) {
    if (length(parts) >= 5L) parts[[5L]] else ""
  }, character(1L))

  # Gamma and delta constants may occupy separate reference rows. Exclude the
  # alternative gamma EX2T/EX2R exons and retain the canonical EX2.
  if (length(index) > 1L) {
    canonical <- features %in% c("EX1", "EX2", "EX3", "EX4", "EX4UTR")
    index <- index[canonical]
    features <- features[canonical]
  }
  if (!length(index)) {
    stop("No canonical coding exons found for ", allele, call. = FALSE)
  }

  dna <- paste0(ref$seq.nt[index], collapse = "")
  protein <- .tcr_translate(dna, frame = 0L)
  stop_position <- regexpr("*", protein, fixed = TRUE)[[1L]]
  if (stop_position > 0L) protein <- substr(protein, 1L, stop_position - 1L)
  protein <- sub("^X+", "", protein)

  if (!nzchar(protein)) {
    stop("Constant-region translation is empty for ", allele,
         call. = FALSE)
  }
  if (strict && grepl("X", protein, fixed = TRUE)) {
    stop("Constant-region translation contains an ambiguous residue for ",
         allele, call. = FALSE)
  }

  list(allele = allele, inferred = inferred, protein = protein,
       features = features)
}

.tcr_v_prefix <- function(ref, index) {
  v_protein <- as.character(ref$seq.aa[index])
  if (is.na(v_protein) || !nzchar(v_protein)) {
    stop("The selected V allele has no translated protein in the reference",
         call. = FALSE)
  }

  germline_tail <- as.character(ref$CDR3[index])
  if (!is.na(germline_tail) && nzchar(germline_tail) &&
      endsWith(v_protein, germline_tail)) {
    prefix <- substr(v_protein, 1L, nchar(v_protein) - nchar(germline_tail))
  } else {
    cysteines <- gregexpr("C", v_protein, fixed = TRUE)[[1L]]
    cysteines <- cysteines[cysteines > nchar(v_protein) - 20L]
    if (!length(cysteines)) {
      stop("Could not locate the conserved CDR3 cysteine in the V allele",
           call. = FALSE)
    }
    prefix <- substr(v_protein, 1L, max(cysteines))
  }

  if (!endsWith(prefix, "C")) {
    stop("V reference prefix does not end at the conserved CDR3 cysteine",
         call. = FALSE)
  }
  prefix
}

.tcr_j_translation <- function(dna) {
  options <- lapply(0:2, function(frame) {
    protein <- .tcr_translate(dna, frame)
    first_stop <- regexpr("*", protein, fixed = TRUE)[[1L]]
    before_stop <- if (first_stop > 0L) {
      substr(protein, 1L, first_stop - 1L)
    } else {
      protein
    }
    score <- 1000L * grepl("[FW]G", before_stop) +
      200L * grepl("[FW]G.G", before_stop) +
      20L * (first_stop < 0L) + nchar(before_stop)
    list(frame = frame, protein = before_stop, score = score)
  })
  scores <- vapply(options, `[[`, numeric(1L), "score")
  best <- options[[which.max(scores)]]
  if (!nzchar(best$protein) || !grepl("[FW]", best$protein)) {
    stop("Could not determine a productive J-segment translation frame",
         call. = FALSE)
  }
  best
}

.tcr_j_overlap <- function(cdr3, j_protein) {
  terminal <- substr(cdr3, nchar(cdr3), nchar(cdr3))
  j_chars <- strsplit(j_protein, "", fixed = TRUE)[[1L]]
  positions <- which(j_chars == terminal)
  if (!length(positions)) {
    stop("CDR3 terminal residue ", terminal,
         " cannot be aligned to the selected J allele", call. = FALSE)
  }

  candidates <- lapply(positions, function(position) {
    maximum <- min(nchar(cdr3), position)
    overlap <- 0L
    for (size in seq_len(maximum)) {
      if (substr(cdr3, nchar(cdr3) - size + 1L, nchar(cdr3)) ==
          substr(j_protein, position - size + 1L, position)) {
        overlap <- size
      }
    }
    list(position = position, overlap = overlap)
  })
  overlaps <- vapply(candidates, `[[`, integer(1L), "overlap")
  selected <- candidates[[which.max(overlaps)]]
  suffix <- if (selected$position < nchar(j_protein)) {
    substr(j_protein, selected$position + 1L, nchar(j_protein))
  } else {
    ""
  }
  list(overlap = selected$overlap, position = selected$position,
       suffix = suffix)
}

#' Build a full T-cell receptor protein chain
#'
#' @param ref An IMGT-style data frame such as the supplied dget() dump.
#' @param v_gene V gene or allele (for example, "TRBV7-9*01"). If the allele
#'   is omitted, *01 is selected when uniquely available. An unavailable TRAV
#'   or TRBV name falls back to the unique closest name in the same family.
#' @param j_gene J gene or allele. An unavailable TRAJ or TRBJ name falls back
#'   to the unique closest name in the same family.
#' @param cdr3 Full amino-acid CDR3 in standard notation, including the
#'   conserved N-terminal C and terminal F or W.
#' @param d_gene Optional D gene or allele for beta/delta chains. It is
#'   validated and reported, not separately inserted into the protein.
#' @param c_gene Optional constant gene or allele. Human alpha, beta, and delta
#'   defaults are inferred; gamma chains require an explicit constant allele.
#' @param strict Reject non-functional alleles, incomplete leaders, ambiguous
#'   translations, and non-canonical CDR3 anchors when TRUE.
#' @return An object of class tcr_protein with sequence, components, segment
#'   calls, inference metadata, and reference translation details.
#' @export
build_tcr_protein <- function(ref, v_gene, j_gene, cdr3, d_gene = NULL,
                              c_gene = NULL, strict = TRUE) {
  .tcr_check_reference(ref)
  strict <- isTRUE(strict)

  cdr3 <- toupper(gsub("\\s+", "", cdr3))
  if (length(cdr3) != 1L || is.na(cdr3) || !nzchar(cdr3)) {
    stop("cdr3 must be one non-empty amino-acid sequence", call. = FALSE)
  }
  if (!grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cdr3)) {
    stop("cdr3 contains a non-standard amino-acid symbol", call. = FALSE)
  }
  if (strict && !startsWith(cdr3, "C")) {
    stop("cdr3 must begin with the conserved cysteine (C)", call. = FALSE)
  }
  terminal <- substr(cdr3, nchar(cdr3), nchar(cdr3))
  if (strict && !terminal %in% c("F", "W")) {
    stop("cdr3 must end with the conserved J residue F or W",
         call. = FALSE)
  }

  v <- .tcr_resolve_segment(ref, v_gene, "v_gene", strict)
  j <- .tcr_resolve_segment(ref, j_gene, "j_gene", strict)
  chain <- .tcr_chain_from_j(j$allele)
  .tcr_validate_v_chain(v$allele, chain)

  d <- NULL
  if (!is.null(d_gene)) {
    d <- .tcr_resolve_segment(ref, d_gene, "d_gene", strict)
    .tcr_validate_d_chain(d$allele, chain)
  }

  if (is.null(c_gene)) c_gene <- .tcr_default_constant(chain, j$allele)
  constant <- .tcr_resolve_constant(ref, c_gene, chain, strict)

  v_prefix <- .tcr_v_prefix(ref, v$index)
  if (strict && !startsWith(v_prefix, "M")) {
    stop("The selected V allele does not contain a complete leader sequence",
         call. = FALSE)
  }
  v_before_cdr3 <- substr(v_prefix, 1L, nchar(v_prefix) - 1L)

  j_translation <- .tcr_j_translation(as.character(ref$seq.nt[j$index]))
  junction <- .tcr_j_overlap(cdr3, j_translation$protein)

  components <- c(
    leader_and_v_before_cdr3 = v_before_cdr3,
    cdr3 = cdr3,
    j_after_cdr3 = junction$suffix,
    constant = constant$protein
  )
  sequence <- paste0(components, collapse = "")

  inferred <- c(
    V = v$inferred,
    D = if (is.null(d)) FALSE else d$inferred,
    J = j$inferred,
    C = constant$inferred
  )
  segments <- c(
    V = v$allele,
    D = if (is.null(d)) NA_character_ else d$allele,
    J = j$allele,
    C = constant$allele
  )
  closest_name_matches <- c(
    V = if (v$used_closest) paste0(v$requested, " -> ", v$allele) else NA_character_,
    D = if (!is.null(d) && d$used_closest) {
      paste0(d$requested, " -> ", d$allele)
    } else {
      NA_character_
    },
    J = if (j$used_closest) paste0(j$requested, " -> ", j$allele) else NA_character_
  )

  result <- list(
    sequence = sequence,
    length_aa = nchar(sequence),
    chain = chain,
    cdr3 = cdr3,
    segments = segments,
    inferred_alleles = inferred,
    closest_name_matches = closest_name_matches,
    components = components,
    j_translation = j_translation$protein,
    j_frame = j_translation$frame,
    j_overlap_aa = junction$overlap,
    constant_features = constant$features,
    note = if (is.null(d)) NULL else
      "D is metadata only; its translated contribution is embedded in CDR3."
  )
  class(result) <- "tcr_protein"
  result
}

print.tcr_protein <- function(x, ...) {
  cat("Full ", x$chain, " TCR protein: ", x$length_aa, " aa\n", sep = "")
  cat("Segments: ", paste(names(x$segments), x$segments,
                          sep = "=", collapse = ", "), "\n", sep = "")
  cat("CDR3: ", x$cdr3, "\n", sep = "")
  cat(x$sequence, "\n", sep = "")
  invisible(x)
}

as_fasta <- function(x, ...) UseMethod("as_fasta")

as_fasta.tcr_protein <- function(x, header = NULL, width = 80L, ...) {
  width <- as.integer(width)
  if (length(width) != 1L || is.na(width) || width < 1L) {
    stop("width must be a positive integer", call. = FALSE)
  }
  if (is.null(header)) {
    segment_text <- paste(na.omit(unname(x$segments)), collapse = "_")
    header <- paste0("TCR_", x$chain, "|", segment_text,
                     "|CDR3=", x$cdr3)
  }
  starts <- seq.int(1L, nchar(x$sequence), by = width)
  lines <- substring(x$sequence, starts,
                     pmin(starts + width - 1L, nchar(x$sequence)))
  paste0(">", header, "\n", paste(lines, collapse = "\n"), "\n")
}
