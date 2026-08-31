#' Prepare an HLA allele data frame from an XML file
#'
#' Reads allele information from an IPD-IMGT/HLA XML file. The input may be
#' uncompressed (`.xml`), gzip-compressed (`.xml.gz`), or ZIP-compressed
#' (`.xml.zip`).
#'
#' ZIP archives must contain a file named `hla.xml`, or exactly one XML file.
#'
#' @param file_path Path to an `.xml`, `.xml.gz`, or `.xml.zip` file.
#' @param lapply_fun Function used to process allele nodes. Suitable choices
#'   include [base::lapply()], `pbapply::pblapply`, and
#'   `parallel::mclapply`.
#' @param replace_none_pg Logical. Replace `"None"` in `p_group` and `g_group`
#'   with `allele_protein` and `allele_coding`, respectively.
#' @param ... Additional arguments passed to `lapply_fun`, such as `mc.cores`.
#'
#' @return A data frame with one row per allele.
#' @export
#'
#' @examples
#' \dontrun{
#' hla <- hla_df_from_xml("hla.xml")
#' hla <- hla_df_from_xml("hla.xml.gz")
#' hla <- hla_df_from_xml("hla.xml.zip")
#'
#' hla <- hla_df_from_xml(
#'   "hla.xml.zip",
#'   lapply_fun = parallel::mclapply,
#'   mc.cores = parallel::detectCores()
#' )
#' }
hla_df_from_xml <- function(file_path,
                            lapply_fun = lapply,
                            replace_none_pg = TRUE,
                            ...) {
  required <- c("xml2", "dplyr")
  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing)) {
    stop(
      "Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(file_path) != 1L ||
      is.na(file_path) ||
      !nzchar(file_path)) {
    stop("`file_path` must be a single, non-empty path.", call. = FALSE)
  }

  file_path <- path.expand(file_path)

  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path, call. = FALSE)
  }

  if (!is.logical(replace_none_pg) ||
      length(replace_none_pg) != 1L ||
      is.na(replace_none_pg)) {
    stop("`replace_none_pg` must be TRUE or FALSE.", call. = FALSE)
  }

  lapply_fun <- match.fun(lapply_fun)
  document <- .read_hla_xml(file_path)

  alleles <- xml2::xml_find_all(
    document,
    ".//*[local-name() = 'allele']"
  )

  if (!length(alleles)) {
    stop("The XML document contains no allele nodes.", call. = FALSE)
  }

  rows <- lapply_fun(alleles, read_child, ...)

  if (!is.list(rows)) {
    stop("`lapply_fun` must return a list.", call. = FALSE)
  }

  rows <- Filter(Negate(is.null), rows)

  if (!length(rows)) {
    stop("No allele records could be parsed.", call. = FALSE)
  }

  df <- dplyr::bind_rows(rows)
  rownames(df) <- NULL

  # Normalize to four allele fields while preserving expression suffixes,
  # e.g. HLA-A*01:01N -> HLA-A*01:01:01:01N.
  df$allele <- .pad_hla_allele(df$allele)

  df$prefix <- substr(df$allele, 1L, 3L)

  allele_name <- sub("^HLA-", "", df$allele)
  mic <- grepl("^MIC[^*]*\\*", allele_name)
  allele_name[mic] <- sub("^MIC", "", allele_name[mic])
  df$gene <- sub("\\*.*$", "", allele_name)

  # Remove the fourth allele field and any expression suffix.
  df$allele_coding <- sub(
    ":[[:digit:]]+[[:alpha:]]*$",
    "",
    df$allele
  )
  df$allele_coding <- sub("^HLA-", "", df$allele_coding)
  df$allele_coding <- sub("^MIC", "", df$allele_coding)

  fields <- strsplit(df$allele_coding, ":", fixed = TRUE)

  df$allele_group <- vapply(
    fields,
    function(x) x[[1L]],
    character(1)
  )

  df$allele_protein <- vapply(
    fields,
    function(x) paste(x[seq_len(min(2L, length(x)))], collapse = ":"),
    character(1)
  )

  df$seq_length <- nchar(df$seq)

  exon2 <- if ("Exon2" %in% names(df)) {
    df$Exon2
  } else {
    rep(NA_character_, nrow(df))
  }

  exon3 <- if ("Exon3" %in% names(df)) {
    df$Exon3
  } else {
    rep(NA_character_, nrow(df))
  }

  exon2_start <- suppressWarnings(
    as.integer(sub("_.*$", "", exon2))
  )
  exon3_end <- suppressWarnings(
    as.integer(sub("^.*_", "", exon3))
  )

  valid <- !is.na(df$seq) &
    !is.na(exon2_start) &
    !is.na(exon3_end) &
    exon2_start >= 1L &
    exon3_end >= exon2_start

  df$seq_Exon2_3 <- NA_character_
  df$seq_Exon2_3[valid] <- substr(
    df$seq[valid],
    exon2_start[valid],
    exon3_end[valid]
  )
  df$seq_Exon2_3_length <- nchar(df$seq_Exon2_3)

  if (replace_none_pg) {
    replace_g <- !is.na(df$g_group) & df$g_group == "None"
    replace_p <- !is.na(df$p_group) & df$p_group == "None"

    df$g_group[replace_g] <- df$allele_coding[replace_g]
    df$p_group[replace_p] <- df$allele_protein[replace_p]
  }

  df
}


# Read an XML document from .xml, .xml.gz, or .xml.zip.
.read_hla_xml <- function(file_path) {
  lower_path <- tolower(file_path)

  if (grepl("\\.xml$", lower_path)) {
    return(xml2::read_xml(file_path))
  }

  if (grepl("\\.xml\\.gz$", lower_path)) {
    connection <- gzfile(file_path, open = "rb")
    on.exit(close(connection), add = TRUE)

    return(xml2::read_xml(connection))
  }

  if (grepl("\\.xml\\.zip$", lower_path)) {
    archive <- utils::unzip(file_path, list = TRUE)

    xml_members <- archive$Name[
      grepl("\\.xml$", archive$Name, ignore.case = TRUE)
    ]

    preferred <- xml_members[
      tolower(basename(xml_members)) == "hla.xml"
    ]

    member <- if (length(preferred) == 1L) {
      preferred
    } else if (length(xml_members) == 1L) {
      xml_members
    } else if (!length(xml_members)) {
      stop("The ZIP archive contains no XML file.", call. = FALSE)
    } else {
      stop(
        paste0(
          "The ZIP archive contains multiple XML files and no unique ",
          "`hla.xml`: ",
          paste(xml_members, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    connection <- unz(file_path, member, open = "rb")
    on.exit(close(connection), add = TRUE)

    return(xml2::read_xml(connection))
  }

  stop(
    "`file_path` must end in .xml, .xml.gz, or .xml.zip.",
    call. = FALSE
  )
}


# Parse a single allele node.
read_child <- function(x) {
  allele <- xml2::xml_attr(x, "name")

  if (is.na(allele) || !nzchar(allele)) {
    stop("Encountered an allele node without a name.", call. = FALSE)
  }

  sequence_node <- xml2::xml_find_first(
    x,
    "./*[local-name() = 'sequence']"
  )

  if (.is_xml_missing(sequence_node)) {
    sequence <- NA_character_
    features <- list()
  } else {
    nucleotide_node <- xml2::xml_find_first(
      sequence_node,
      "./*[local-name() = 'nucsequence']"
    )

    sequence <- if (.is_xml_missing(nucleotide_node)) {
      NA_character_
    } else {
      gsub(
        "[[:space:]]+",
        "",
        xml2::xml_text(nucleotide_node)
      )
    }

    features <- xml2::xml_find_all(
      sequence_node,
      "./*[local-name() = 'feature']"
    )
  }

  result <- data.frame(
    allele = allele,
    seq = sequence,
    g_group = .read_hla_group(x, "hla_g_group"),
    p_group = .read_hla_group(x, "hla_p_group"),
    stringsAsFactors = FALSE
  )

  if (length(features)) {
    feature_names <- xml2::xml_attr(features, "name")

    keep <- !is.na(feature_names) &
      feature_names != "Translation"

    features <- features[keep]
    feature_names <- gsub(
      "[[:space:]]+",
      "",
      feature_names[keep]
    )

    for (i in seq_along(features)) {
      # Keep the first occurrence if a feature name is duplicated.
      if (!feature_names[[i]] %in% names(result)) {
        result[[feature_names[[i]]]] <-
          .read_feature_bounds(features[[i]])
      }
    }
  }

  result
}


.read_hla_group <- function(allele_node, group_name) {
  node <- xml2::xml_find_first(
    allele_node,
    sprintf("./*[local-name() = '%s']", group_name)
  )

  if (.is_xml_missing(node)) {
    return("None")
  }

  attributes <- xml2::xml_attrs(node)

  value <- if ("name" %in% names(attributes)) {
    unname(attributes[["name"]])
  } else if (length(attributes)) {
    unname(attributes[[1L]])
  } else {
    trimws(xml2::xml_text(node))
  }

  if (!length(value) || is.na(value) || !nzchar(value)) {
    "None"
  } else {
    value
  }
}


.read_feature_bounds <- function(feature) {
  coordinates <- xml2::xml_find_first(
    feature,
    "./*[local-name() = 'coordinates']"
  )

  if (.is_xml_missing(coordinates)) {
    return(NA_character_)
  }

  start <- xml2::xml_attr(coordinates, "start")
  end <- xml2::xml_attr(coordinates, "end")

  if (!is.na(start) && !is.na(end)) {
    return(paste(start, end, sep = "_"))
  }

  attributes <- unname(xml2::xml_attrs(coordinates))

  if (length(attributes) >= 2L) {
    paste(attributes[[1L]], attributes[[2L]], sep = "_")
  } else {
    NA_character_
  }
}


.pad_hla_allele <- function(x) {
  result <- x
  valid <- !is.na(x)

  suffix <- ifelse(
    grepl("[NLSCAQ]$", x[valid]),
    sub("^.*([NLSCAQ])$", "\\1", x[valid]),
    ""
  )

  stem <- ifelse(
    nzchar(suffix),
    substr(x[valid], 1L, nchar(x[valid]) - 1L),
    x[valid]
  )

  colon_count <- nchar(gsub("[^:]", "", stem))
  missing_fields <- pmax.int(0L, 3L - colon_count)

  padding <- vapply(
    missing_fields,
    function(n) paste(rep(":01", n), collapse = ""),
    character(1)
  )

  result[valid] <- paste0(stem, padding, suffix)
  result
}


.is_xml_missing <- function(x) {
  length(x) == 0L || inherits(x, "xml_missing")
}


#' Extract and optionally translate concatenated exon sequences
#'
#' Extracts specified exons from genomic sequences, removes intervening
#' introns, and concatenates the exon sequences. The resulting nucleotide
#' sequence can optionally be restricted to a subsequence and translated.
#'
#' @param data A data frame containing a character column named `seq` and exon
#'   coordinate columns named `Exon1`, `Exon2`, etc. Coordinates must use the
#'   one-based inclusive format `"start_end"`.
#' @param exons An optional positive integer vector specifying the exons to
#'   extract. Exons are concatenated in the supplied order. If `NULL`, all
#'   available exons are used in numeric order.
#' @param cds_start An optional positive integer giving the first nucleotide
#'   to retain from the concatenated sequence.
#' @param cds_end An optional positive integer giving the last nucleotide
#'   to retain from the concatenated sequence.
#' @param translate A single logical value indicating whether the extracted
#'   sequence should be translated.
#' @param no.init.codon A logical value passed to [Biostrings::translate()].
#'   When Exon1 is missing and the reading frame is inferred, `TRUE` is used
#'   because the extracted sequence is an internal CDS fragment.
#' @param if.fuzzy.codon Controls how [Biostrings::translate()] handles fuzzy
#'   codons. The default, `"X"`, translates fuzzy codons to `"X"`.
#'
#' @return A character vector containing one nucleotide or amino-acid sequence
#'   per row of `data`.
#'
#' @details
#' A missing Exon1 is allowed. When `translate = TRUE` and Exon1 is missing,
#' all three possible reading frames are tested. Frames containing a stop codon
#' before the final translated position are rejected. Among acceptable frames,
#' a frame ending in a terminal stop codon is preferred; otherwise, the longest
#' translation is selected. Ties are resolved using the smallest nucleotide
#' offset and generate a warning.
#'
#' Missing trailing exon annotations are allowed. However, a missing exon from
#' Exon2 through the highest available exon is treated as an internal gap and
#' produces an error.
#'
#' @examples
#' alleles <- data.frame(
#'   seq = "ATGAAAGGGCCCTAA",
#'   Exon1 = "1_6",
#'   Exon2 = "10_15"
#' )
#'
#' extract_cds(alleles)
#'
#' extract_cds(
#'   alleles,
#'   translate = TRUE,
#'   no.init.codon = FALSE
#' )
#'
#' @export
extract_cds <- function(data,
                        exons = NULL,
                        cds_start = NULL,
                        cds_end = NULL,
                        translate = FALSE,
                        no.init.codon = TRUE,
                        if.fuzzy.codon = "X") {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!"seq" %in% names(data)) {
    stop("`data` must contain a column named `seq`.", call. = FALSE)
  }

  if (!is.character(data$seq)) {
    stop("`data$seq` must be a character vector.", call. = FALSE)
  }

  validate_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  validate_position <- function(x, name) {
    if (is.null(x)) {
      return(invisible(NULL))
    }

    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < 1 ||
        x != floor(x)) {
      stop(
        "`", name, "` must be NULL or one positive integer.",
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  validate_logical(translate, "translate")
  validate_logical(no.init.codon, "no.init.codon")
  validate_position(cds_start, "cds_start")
  validate_position(cds_end, "cds_end")

  exon_cols <- grep(
    "^Exon[0-9]+$",
    names(data),
    value = TRUE
  )

  if (!length(exon_cols)) {
    stop("No exon columns were found.", call. = FALSE)
  }

  exon_numbers <- as.integer(sub("^Exon", "", exon_cols))
  exon_order <- order(exon_numbers)

  exon_cols <- exon_cols[exon_order]
  exon_numbers <- exon_numbers[exon_order]

  if (anyDuplicated(exon_numbers)) {
    stop("Duplicate exon numbers were found.", call. = FALSE)
  }

  explicit_exons <- !is.null(exons)

  if (explicit_exons) {
    if (!is.numeric(exons) ||
        !length(exons) ||
        anyNA(exons) ||
        any(!is.finite(exons)) ||
        any(exons < 1) ||
        any(exons != floor(exons))) {
      stop(
        "`exons` must be NULL or a vector of positive integers.",
        call. = FALSE
      )
    }

    if (anyDuplicated(exons)) {
      stop("`exons` must not contain duplicates.", call. = FALSE)
    }

    exon_cols_req <- paste0("Exon", as.integer(exons))

    missing_columns <- setdiff(exon_cols_req, names(data))
    disallowed_missing <- setdiff(
      missing_columns,
      if (translate) "Exon1" else character()
    )

    if (length(disallowed_missing)) {
      stop(
        "Missing exon columns: ",
        paste(disallowed_missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  } else {
    exon_cols_req <- exon_cols
  }

  if (translate && !requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
      "Package `Biostrings` is required when `translate = TRUE`.",
      call. = FALSE
    )
  }

  translate_inferred_frame <- function(cds, row_number) {
    translate_candidate <- function(offset) {
      available_length <- nchar(cds) - offset
      translated_width <- available_length -
        (available_length %% 3L)

      if (translated_width < 3L) {
        return(NULL)
      }

      candidate_dna <- substr(
        cds,
        offset + 1L,
        offset + translated_width
      )

      candidate_aa <- as.character(
        Biostrings::translate(
          Biostrings::DNAString(toupper(candidate_dna)),
          no.init.codon = TRUE,
          if.fuzzy.codon = if.fuzzy.codon
        )
      )

      aa_letters <- strsplit(
        candidate_aa,
        "",
        fixed = TRUE
      )[[1]]

      stop_positions <- which(aa_letters == "*")
      terminal_position <- length(aa_letters)

      premature_stop <- any(
        stop_positions < terminal_position
      )

      terminal_stop <- (
        length(aa_letters) > 0L &&
          aa_letters[[terminal_position]] == "*"
      )

      list(
        offset = offset,
        aa = candidate_aa,
        length = length(aa_letters),
        premature_stop = premature_stop,
        terminal_stop = terminal_stop
      )
    }

    candidates <- lapply(0:2, translate_candidate)

    acceptable <- which(vapply(
      candidates,
      function(candidate) {
        !is.null(candidate) && !candidate$premature_stop
      },
      logical(1)
    ))

    if (!length(acceptable)) {
      stop(
        "Cannot infer a reading frame in row ", row_number,
        ": all three frames contain a premature stop codon.",
        call. = FALSE
      )
    }

    terminal_candidates <- acceptable[vapply(
      candidates[acceptable],
      function(candidate) candidate$terminal_stop,
      logical(1)
    )]

    candidate_pool <- if (length(terminal_candidates)) {
      terminal_candidates
    } else {
      acceptable
    }

    candidate_lengths <- vapply(
      candidates[candidate_pool],
      function(candidate) candidate$length,
      integer(1)
    )

    longest_candidates <- candidate_pool[
      candidate_lengths == max(candidate_lengths)
    ]

    candidate_offsets <- vapply(
      candidates[longest_candidates],
      function(candidate) candidate$offset,
      integer(1)
    )

    selected <- longest_candidates[
      which.min(candidate_offsets)
    ]

    if (length(candidate_pool) > 1L) {
      warning(
        "Multiple acceptable reading frames were found in row ",
        row_number, "; selected nucleotide offset ",
        candidates[[selected]]$offset, ".",
        call. = FALSE
      )
    }

    candidates[[selected]]$aa
  }

  extract_one <- function(i) {
    sequence <- data$seq[[i]]

    if (is.na(sequence) || !nzchar(sequence)) {
      stop(
        "Missing or empty sequence in row ", i, ".",
        call. = FALSE
      )
    }

    # Create a complete coordinate vector so absent exon columns can also
    # be detected as internal gaps.
    maximum_exon <- max(c(
      exon_numbers,
      if (explicit_exons) as.integer(exons) else integer()
    ))

    all_ranges <- setNames(
      rep(NA_character_, maximum_exon),
      paste0("Exon", seq_len(maximum_exon))
    )

    available_columns <- intersect(
      names(all_ranges),
      exon_cols
    )

    all_ranges[available_columns] <- unlist(
      data[i, available_columns, drop = FALSE],
      use.names = FALSE
    )

    available <- !is.na(all_ranges) & nzchar(all_ranges)
    available_numbers <- which(available)

    if (!length(available_numbers)) {
      return(NA_character_)
    }

    highest_available <- max(available_numbers)

    # Exon1 may be absent, but no gaps are permitted after Exon22 begins.
    required_internal <- if (highest_available >= 2L) {
      seq.int(2L, highest_available)
    } else {
      integer()
    }

    missing_internal <- required_internal[
      !available[required_internal]
    ]

    if (length(missing_internal)) {
      stop(
        "Missing internal exon coordinates in row ", i, ": ",
        paste0(
          "Exon",
          missing_internal,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    exon1_missing <- !available[1L]

    requested_ranges <- all_ranges[exon_cols_req]
    requested_unavailable <- (
      is.na(requested_ranges) |
        !nzchar(requested_ranges)
    )

    if (explicit_exons && any(requested_unavailable)) {
      unavailable_names <- names(requested_ranges)[
        requested_unavailable
      ]

      allowed_missing <- (
        translate &&
          unavailable_names == "Exon1"
      )

      disallowed_names <- unavailable_names[
        !allowed_missing
      ]

      if (length(disallowed_names)) {
        stop(
          "No coordinates available for ",
          paste(disallowed_names, collapse = ", "),
          " in row ", i, ".",
          call. = FALSE
        )
      }
    }

    requested_ranges <- requested_ranges[
      !is.na(requested_ranges) &
        nzchar(requested_ranges)
    ]

    if (!length(requested_ranges)) {
      stop(
        "No extractable exon sequence remains in row ", i, ".",
        call. = FALSE
      )
    }

    exon_sequences <- vapply(
      names(requested_ranges),
      function(exon_name) {
        range <- requested_ranges[[exon_name]]

        if (!grepl("^[0-9]+_[0-9]+$", range)) {
          stop(
            "Invalid coordinate for ", exon_name,
            " in row ", i, ": ", range, ".",
            call. = FALSE
          )
        }

        coordinate <- as.integer(
          strsplit(range, "_", fixed = TRUE)[[1]]
        )

        start <- coordinate[1]
        end <- coordinate[2]
        sequence_length <- nchar(sequence)

        if (start < 1L ||
            end < start ||
            end > sequence_length) {
          stop(
            "Coordinate ", range, " for ", exon_name,
            " is outside the sequence in row ", i,
            " (sequence length: ", sequence_length, ").",
            call. = FALSE
          )
        }

        substr(sequence, start, end)
      },
      character(1),
      USE.NAMES = TRUE
    )

    cds <- paste0(exon_sequences, collapse = "")

    from <- if (is.null(cds_start)) {
      1L
    } else {
      as.integer(cds_start)
    }

    to <- if (is.null(cds_end)) {
      nchar(cds)
    } else {
      as.integer(cds_end)
    }

    if (to < from || to > nchar(cds)) {
      stop(
        "Requested CDS positions ", from, "-", to,
        " are outside the concatenated sequence in row ", i,
        " (length: ", nchar(cds), ").",
        call. = FALSE
      )
    }

    cds <- substr(cds, from, to)

    if (!translate) {
      return(cds)
    }

    if (exon1_missing) {
      return(translate_inferred_frame(cds, i))
    }

    # if (nchar(cds) %% 3L != 0L) {
    #   stop(
    #     "Cannot translate row ", i, ": extracted sequence length ",
    #     nchar(cds), " is not divisible by three.",
    #     call. = FALSE
    #   )
    # }

    as.character(
      Biostrings::translate(
        Biostrings::DNAString(toupper(cds)),
        no.init.codon = no.init.codon,
        if.fuzzy.codon = if.fuzzy.codon
      )
    )
  }

  vapply(
    seq_len(nrow(data)),
    extract_one,
    character(1),
    USE.NAMES = FALSE
  )
}


#' Prepare a data frame from a xml file containing hla allele information
#'
#' The xml file required for this function may be downloaded from the GitHub repository of the IMGT/HLA organization:
#' https://github.com/ANHIG/IMGTHLA/blob/Latest/xml/hla.xml.zip. The function then creates a dataframe from this xml
#' which can be used subsequently to match reads of HLA-loci from bulk or single cell RNA seq. By default the sequence
#' of the 2nd + 3rd exon are returned in an extra column as these will most likely be used to infer the HLA type in case of MHC-I.
#'
#' HLA nomenclature: http://hla.alleles.org/nomenclature/naming.html
#'
#' @param file_path path to the xml file, may be in zipped format
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param replace_none_pg replace "None" in p_group and g_group with allele_protein and allele_coding, respectively
#' @param ... additional argument to the lapply function; e.g. mc.cores may be passed when parallel::mclapply is chosen above
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' hla_df_from_xml(file_path = "path/hla.xml", lapply_fun = parallel::mclapply,
#' mc.cores = parallel::detectCores())
#' }
hla_df_from_xml_legacy <- function(file_path,
                            lapply_fun = lapply,
                            replace_none_pg = T,
                            ...) {

  if (!requireNamespace("xml2", quietly = T)) {
    utils::install.packages("xml2")
  }

  if (grepl("zip$", file_path)) {
    utils::unzip(file_path, exdir = tempdir())
    file_path <- file.path(tempdir(), "hla.xml")
  }

  lapply_fun <- match.fun(lapply_fun)
  children <- xml2::xml_children(xml2::xml_child(xml2::read_xml(file_path)))
  df <- dplyr::bind_rows(lapply_fun(children, read_child, ...))
  rownames(df) <- NULL

  df$allele <- unlist(lapply(df$allele, function(x) {
    while(stringr::str_count(x, ":") < 3) {
      x <- paste0(x, ":01")
    }
    return(x)
  }))

  df$prefix <- substr(df$allele, 1,3)
  df$gene <- stringr::str_sub(stringr::str_extract(df$allele, "-.{1,}\\*"), 2, -2)
  df$gene[which(is.na(df$gene))] <- substr(gsub("-", "", df$allele[which(is.na(df$gene))]), 4,4) # for MIC

  df$allele_coding <- gsub(":[[:digit:]]{2}[[:alpha:]]{0,}$", "", df$allele)
  df$allele_coding <- gsub("MIC", "", gsub("HLA-", "", df$allele_coding))
  df$allele_group <- sapply(strsplit(df$allele_coding, ":"), "[", 1)
  df$allele_protein <- paste0(df$allele_group, ":", sapply(strsplit(df$allele_coding, ":"), "[", 2))
  df$seq_length <- nchar(df$seq)

  df$seq_Exon2_3 <- substr(df$seq, as.numeric(sapply(strsplit(df$Exon2, "_"), "[", 1)), as.numeric(sapply(strsplit(df$Exon3, "_"), "[", 2)))
  df$seq_Exon2_3_length <- nchar(df$seq_Exon2_3)

  if (replace_none_pg) {
    df$g_group <- ifelse(df$g_group == "None", df$allele_coding, df$g_group)
    df$p_group <- ifelse(df$p_group == "None", df$allele_protein, df$p_group)
  }

  return(df)
}

read_child <- function(x) {
  cc <- xml2::xml_children(x)
  attrs <- xml2::xml_attrs(x)
  tryCatch({
    seq <- as.character(xml2::xml_contents(xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]]))[1])
    g_group <- xml2::xml_attrs(cc[[which(xml2::xml_name(cc) == "hla_g_group")]])
    p_group <- xml2::xml_attrs(cc[[which(xml2::xml_name(cc) == "hla_p_group")]])

    # features
    ff <- xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]])[which(xml2::xml_name(xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]])) == "feature")]
    type <- sapply(xml2::xml_attrs(ff), function(x) x[["name"]])
    bounds <- sapply(ff, function(x) xml2::xml_attrs(xml2::xml_child(x, 1)))[which(type != "Translation")]
    bounds <- sapply(bounds, function(x) paste(x, collapse = "_"))
    type <- gsub(" ", "", (type[which(type != "Translation")]))
    df <- data.frame(allele = attrs[["name"]], seq = seq, g_group = g_group, p_group = p_group)
    df <- cbind(df, tidyr::pivot_wider(data.frame(type, bounds), names_from = type, values_from = bounds))
    return(df)
  }, error = function(e) {
    return(NULL)
  })
}


