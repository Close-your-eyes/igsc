#' Read and process a GTF file
#'
#' Read a Gene Transfer Format (GTF) file into memory, optionally restricting
#' the input by sequence name, feature type, or gene name before parsing its
#' attribute column. Files may be plain text or gzip-compressed.
#'
#' Parsing the attribute column is usually the most expensive step. Applying
#' \code{seqnames}, \code{features}, or \code{gene_names} can therefore reduce
#' both run time and memory use.
#'
#' @param file_path A length-one character string giving the path to a GTF file.
#'   Gzip-compressed files with a \code{.gz} extension are supported.
#' @param seqnames An optional character vector of sequence names to read, such
#'   as \code{"chr1"} or \code{c("chr1", "chrX")}. When supplied, matching
#'   regions are selected before the attribute column is processed.
#' @param features An optional character vector of feature types to retain. Each
#'   value must occur in the file after sequence filtering. Common values include
#'   \code{"gene"}, \code{"transcript"}, \code{"exon"}, \code{"CDS"},
#'   \code{"start_codon"}, and \code{"stop_codon"}.
#' @param gene_names An optional character vector used to filter the raw
#'   \code{attribute} field before it is parsed. Matching is case-insensitive.
#' @param gene_names_full_match Logical. If \code{TRUE}, each value in
#'   \code{gene_names} is matched as a complete quoted GTF attribute value. If
#'   \code{FALSE}, values are treated as case-insensitive regular expressions
#'   and may match substrings.
#' @param process_attr_col Logical. If \code{TRUE}, parse the GTF attribute field
#'   with \code{process_gtf_attribute_col()}; if \code{FALSE}, retain the raw
#'   \code{attribute} column and return \code{NULL} for the parsed attributes.
#' @param process_attr_col_args A named list of arguments passed to
#'   \code{process_gtf_attribute_col()}. The default keeps common gene,
#'   transcript, and exon identifiers.
#' @param process_attr_col_args_repl A named list of values that replace entries
#'   in \code{process_attr_col_args}. Use this to change selected defaults
#'   without repeating the complete argument list.
#' @param col_names A character vector of column names assigned to the nine GTF
#'   fields. The default follows the GTF specification. Internal processing
#'   requires columns named \code{seqname}, \code{feature}, and
#'   \code{attribute}; changing those names will break the function.
#'
#' @details
#' GTF coordinates are one-based and inclusive. For records on either strand,
#' \code{start} is numerically less than or equal to \code{end}. Biological
#' transcription therefore proceeds from larger to smaller coordinates for a
#' feature on the minus strand.
#'
#' Duplicate rows are removed before the result is returned. If attribute
#' processing and its requested transformations remove all rows, the function
#' may return \code{NULL}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{gtf}}{A data frame containing the selected GTF records. If
#'     \code{process_attr_col = TRUE}, retained attributes are joined as
#'     separate columns; otherwise the original \code{attribute} column is
#'     preserved.}
#'   \item{\code{attr}}{A wide data frame of parsed, retained attributes, or
#'     \code{NULL} when \code{process_attr_col = FALSE}.}
#' }
#'
#' @references
#' Ensembl, "GFF/GTF File Format":
#' \url{https://www.ensembl.org/info/website/upload/gff.html}
#' \url{https://www.ensembl.org/info/website/upload/gff.html?redirect=no}
#' \url{https://www.ncbi.nlm.nih.gov/datasets/genome/}
#' \url{https://www.ensembl.org/index.html}
#'
#' @seealso \code{process_gtf_attribute_col()}
#'
#' @export
#'
#' @importFrom rlang :=
#'
#' @examples
#' \dontrun{
#' # Read selected feature types and parse common attributes.
#' x <- read_gtf(
#'   "genes.gtf.gz",
#'   features = c("gene", "exon")
#' )
#' head(x$gtf)
#' head(x$attr)
#'
#' # Read one sequence while preserving the raw attribute field.
#' chr1 <- read_gtf(
#'   "genes.gtf",
#'   seqnames = "chr1",
#'   process_attr_col = FALSE
#' )
#' head(chr1$gtf$attribute)
#'
#' # Override selected attribute-processing defaults.
#' exons <- read_gtf(
#'   "genes.gtf",
#'   process_attr_col_args_repl = list(
#'     exons_only = TRUE,
#'     attr_keep = c("gene_id", "gene_name", "transcript_id", "exon_number")
#'   )
#' )
#' }
read_gtf <- function(
    file_path,
    seqnames = NULL,
    features = NULL,
    gene_names = NULL,
    gene_names_full_match = T,
    process_attr_col = T,
    process_attr_col_args = list(
      attr_keep = c(
        "gene_id",
        "gene_name",
        "transcript_id",
        "transcript_name",
        "exon_number"
      ),
      attr_rename = NULL,
      gene_name_prefix = "",
      gene_name_force = NULL,
      rm_index = T,
      use_fun = "r",
      rm_exon = F,
      exons_only = F,
      features_to_exon = c(""),
      aggregate_exons = F,
      aggregate_overlapping_exon_ranges = F,
      check_for_rotation = F,
      genome_length = NULL,
      fill_na = F
    ),
    process_attr_col_args_repl = list(),

    col_names = c(
      "seqname",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attribute"
    )
) {

  ## add fst folder option


  if (missing(file_path)) {
    stop("Provide file_path or gtf data frame.")
  }
  if (!file.exists(file_path)) {
    stop(file_path, " not found.")
  }

  if (tools::file_ext(file_path) == "gz") {
    unpack_fun <- gzfile
  } else {
    unpack_fun <- function(description) {description}
  }

  if (!is.null(seqnames)) {
    gtf <- do.call(rbind, lapply(
      seqnames,
      vroom_gtf,
      file_path = file_path,
      col_names = col_names,
      unpack_fun = unpack_fun
    ))
  } else {
    gtf <- vroom::vroom(
      file = do.call(unpack_fun, args = list(description = file_path)),
      col_names = col_names,
      comment = "#",
      delim = "\t",
      show_col_types = F,
      progress = F
    )
  }

  if (!is.null(features)) {
    features <- rlang::arg_match(features, values = unique(gtf$feature), multiple = T)
    gtf <- gtf[which(gtf$feature %in% features),]
  }

  if (!is.null(gene_names)) {
    if (gene_names_full_match) {
      gene_names <- paste0("\"", gene_names, "\"")
    }
    pattern <- paste(gene_names, collapse="|")
    gtf <- gtf[stringi::stri_detect_regex(gtf$attribute, pattern, case_insensitive = TRUE), ]

    # data.table::setDT(gtf)
    # # Extract gene_id once (fast C-level regex)
    # gtf[, gene_id := stringi::stri_match_first_regex(
    #   attribute,
    #   'gene_id "([^"]+)"'
    # )[,2]]
    # # Exact matching (very fast hash match)
    # gtf2 <- gtf[gene_id %in% gene_names]

  }

  if (nrow(gtf) == 0) {
    stop("No rows left in gtf, check features argument.")
  }

  if (process_attr_col) {
    # message("processing the attribute column.")
    # use waldo::compare to compare results
    for (i in names(process_attr_col_args_repl)) {
      process_attr_col_args[[i]] <-  process_attr_col_args_repl[[i]]
    }

    ret_list <- Gmisc::fastDoCall(what = process_gtf_attribute_col,
                                  args = c(list(gtf = gtf),
                                           process_attr_col_args))

    if (is.null(ret_list)) {
      return(NULL)
    }
    if (anyDuplicated(ret_list$gtf)) {
      ret_list$gtf <- ret_list$gtf |> dplyr::distinct()
      message("duplicate rows in gtf. made unique.")
    }
    return(ret_list)
  } else {
    if (anyDuplicated(gtf)) {
      gtf <- gtf |> dplyr::distinct()
      message("duplicate rows in gtf. made unique.")
    }
    return(list(gtf = gtf, attr = NULL))
  }

}

#' Parse and transform a GTF attribute column
#'
#' Parse the ninth GTF field into separate attributes and optionally normalize
#' identifiers, reshape exon records, merge overlapping ranges, or rotate
#' coordinates for a circular reference. The result can retain attributes as
#' columns or reconstruct a standard GTF key-value attribute field.
#'
#' @param gtf A data frame containing the nine standard GTF fields, including
#'   \code{seqname}, \code{feature}, \code{start}, \code{end},
#'   \code{strand}, and \code{attribute}. This is typically
#'   \code{read_gtf(..., process_attr_col = FALSE)$gtf}.
#' @param attr_keep A character vector of attribute names to retain, for example
#'   \code{c("gene_id", "transcript_id", "gene_name", "gene_type")}.
#'   Use \code{NULL} to retain every parsed attribute. Destination names in
#'   \code{attr_rename} are retained automatically.
#' @param attr_rename An optional named character vector mapping source attribute
#'   names to destination names, for example
#'   \code{c(product = "gene_name", gene_biotype = "gene_type")}.
#' @param rename_replace Logical. If \code{FALSE}, mapped attributes are copied
#'   under their destination names and the source attributes are retained. If
#'   \code{TRUE}, source names are replaced. Pre-existing destination
#'   attributes are removed before the mapping is applied.
#' @param attr_as Output representation for attributes. \code{"cols"} adds one
#'   column per retained attribute to \code{gtf}; \code{"kv"} reconstructs
#'   the GTF \code{attribute} field as quoted key-value pairs.
#' @param gene_name_prefix A character string prepended to non-missing
#'   \code{gene_name} values. It is ignored when \code{gene_name_force = NULL}.
#' @param gene_name_force The name of an attribute used when \code{gene_name}
#'   is missing. The named attribute must be retained and present. Set to
#'   \code{NULL} to disable fallback; a non-empty \code{gene_name_prefix} is
#'   then reset to an empty string.
#' @param gene_name_fix Logical. If \code{TRUE}, replace each run of
#'   non-alphanumeric characters in \code{gene_name} with \code{"-"}, then
#'   remove leading and trailing hyphens.
#' @param gene_name_replace An optional named character vector or list of
#'   regular-expression replacements applied sequentially to \code{gene_name}.
#'   Names are patterns and values are replacements, for example
#'   \code{c("hypothetical protein" = "HYPPROT")}.
#' @param use_fun Attribute parser to use. Only \code{"r"} is currently
#'   supported.
#' @param rm_index Logical. If \code{TRUE}, remove the temporary row-linkage
#'   \code{index} column from the returned \code{gtf} component.
#' @param rm_exon Logical. If \code{TRUE}, remove existing exon records before
#'   applying \code{features_to_exon}.
#' @param features_to_exon A character vector of feature types to relabel as
#'   \code{"exon"}, such as \code{"CDS"}.
#' @param exons_only Logical. If \code{TRUE}, retain only exon records after
#'   removal and relabeling. The function returns \code{NULL} if none remain.
#' @param aggregate_exons Logical. If \code{TRUE}, collapse all exons belonging
#'   to each \code{transcript_id} to its minimum start and maximum end. This
#'   discards intron boundaries and sets \code{exon_number} to 1.
#' @param aggregate_overlapping_exon_ranges Logical. If \code{TRUE}, merge
#'   connected overlapping exon ranges within each sequence. Same-strand ranges
#'   retain their strand; a merged range containing opposing strands is assigned
#'   \code{"+"}. Using \code{aggregate_exons = TRUE} first is strongly
#'   recommended so individual exons from different genes are not joined into
#'   an invalid annotation. Requires the \pkg{igraph} package.
#' @param check_for_rotation Logical. For a circular reference, detect a gene
#'   spanning the artificial origin, choose a cut in an intergenic gap, and
#'   rotate GTF coordinates. The corresponding genome sequence must be rotated
#'   separately using the cut position reported in the message.
#' @param genome_length A positive numeric value giving the reference sequence
#'   length. Required when \code{check_for_rotation = TRUE}.
#' @param rm_entries_wo_matching_exon Logical. If \code{TRUE}, retain records
#'   whose \code{transcript_id} occurs in an exon record, plus gene records with
#'   a matching \code{gene_name}. This can be used to satisfy Cell Ranger
#'   \code{mkref} requirements.
#' @param verbose Logical. If \code{TRUE}, emit informational messages where
#'   supported.
#' @param fill_na Logical. If \code{TRUE}, fill missing \code{transcript_id}
#'   from \code{gene_id}, missing \code{transcript_name} from
#'   \code{gene_name}, and missing retained \code{exon_number} values with 0.
#'
#' @details
#' Transformations occur in the following order: existing exons are optionally
#' removed; selected feature types are relabeled as exons; non-exons are
#' optionally discarded; attributes are parsed, selected, and renamed; gene
#' names and duplicate identifiers are normalized; circular coordinates are
#' optionally rotated; exons and overlapping ranges are optionally aggregated;
#' entries without matching exons are optionally removed; and the requested
#' attribute representation is produced.
#'
#' The parser expects conventional GTF attributes of the form
#' \code{key "value";} separated by a semicolon and a space. Repeated values
#' for the same attribute and input row are joined with \code{"__"}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{gtf}}{The processed annotation data frame. Attributes are
#'     returned as individual columns when \code{attr_as = "cols"}, or in a
#'     reconstructed \code{attribute} column when \code{attr_as = "kv"}.}
#'   \item{\code{attr}}{A wide data frame of the parsed attributes selected by
#'     \code{attr_keep}, keyed by the temporary input-row index. This component
#'     reflects parsing and renaming but not later coordinate or exon
#'     transformations.}
#' }
#' The function returns \code{NULL} instead if \code{exons_only = TRUE} and
#' no exon records remain.
#'
#' @seealso \code{read_gtf()}, \code{make_kv_attr_col()}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse selected attributes into columns.
#' raw_gtf <- read_gtf("genes.gtf", process_attr_col = FALSE)$gtf
#' parsed <- process_gtf_attribute_col(
#'   raw_gtf,
#'   attr_keep = c("gene_id", "gene_name", "transcript_id", "exon_number")
#' )
#' head(parsed$gtf)
#'
#' # Copy gene_biotype to gene_type and reconstruct a valid attribute field.
#' renamed <- process_gtf_attribute_col(
#'   raw_gtf,
#'   attr_rename = c(gene_biotype = "gene_type"),
#'   rename_replace = FALSE,
#'   attr_as = "kv"
#' )
#'
#' # Build exon-only records from CDS features.
#' cds_exons <- process_gtf_attribute_col(
#'   raw_gtf,
#'   attr_keep = c("gene_id", "gene_name", "transcript_id", "exon_number"),
#'   rm_exon = TRUE,
#'   features_to_exon = "CDS",
#'   exons_only = TRUE,
#'   attr_as = "kv"
#' )
#' }
process_gtf_attribute_col <- function(gtf,
                                      attr_keep = NULL,
                                      attr_rename = NULL,
                                      rename_replace = F,
                                      attr_as = c("cols", "kv"),
                                      gene_name_prefix = "",
                                      gene_name_force = "gene_id",
                                      gene_name_fix = T,
                                      gene_name_replace = list(),
                                      use_fun = "r",
                                      rm_index = T,
                                      rm_exon = F,
                                      exons_only = F,
                                      features_to_exon = c(""),
                                      aggregate_exons = F,
                                      aggregate_overlapping_exon_ranges = F,
                                      check_for_rotation = F,
                                      genome_length = NULL,
                                      rm_entries_wo_matching_exon = F,
                                      verbose = T,
                                      fill_na = F) {

  attr_as <- rlang::arg_match(attr_as) # kv = key value pair
  use_fun <- rlang::arg_match(use_fun)

  if (check_for_rotation && is.null(genome_length)) {
    stop("genome_length needed when check_for_rotation.")
  }

  if (gene_name_prefix != "" && is.null(gene_name_force)) {
    message("setting gene_name_prefix to \"\" as gene_name_force is NULL")
    gene_name_prefix <- ""
  }

  if (rm_exon) {
    gtf <- dplyr::filter(gtf, feature != "exon")
  }
  gtf$feature[which(gtf$feature %in% features_to_exon)] <- "exon"

  if (exons_only) {
    before <- nrow(gtf)
    gtf <- dplyr::filter(gtf, feature == "exon")
    if (nrow(gtf) == 0) {
      if (verbose) {
        message("no rows found with exon.")
      }
      return(NULL)
    } else {
      if (verbose) {
        message("all gtf rows: ", formatC(before, big.mark = ","), ". exons only: ", formatC(nrow(gtf), big.mark = ","), ".")
      }
    }
  }

  if (use_fun == "rcpp") {
    # rcpp fun is currently slower than the r procedure
    attr_ind <- rep(seq_along(gtf$attribute), lengths(stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)))
    attr_col <- unlist(processStrings(gtf$attribute)) #igsc:::
  } else if (use_fun == "r") {
    # optimized for speed
    attr_col <- stringi::stri_replace_last(gtf$attribute, replacement = "", fixed = ";")
    #attr_col <- stringi::stri_split_fixed(attr_col, pattern = "; ", omit_empty = T) # this failed when there were ';' in note
    attr_col <- stringi::stri_split_fixed(attr_col, pattern = "\"; ", omit_empty = T)
    attr_ind <- rep(seq_along(attr_col), lengths(attr_col))
    # n = 2 --> only split at first space; tag name must not have spaces, tag content may have spaces
    attr_col <- unlist(stringi::stri_split_fixed(unlist(attr_col), pattern = " ", omit_empty = T, n = 2)) #pattern = ' \"'
    attr_col[seq(2, length(attr_col), 2)] <- stringi::stri_replace_all(attr_col[seq(2, length(attr_col), 2)], replacement = "", fixed = '"') #stringi::stri_trim_both
  } else if (use_fun == "rust") {
    # rust fun was found slower; so using rcpp fun with proper registration
    stop("rust fun is not up to date. e.g. splitting at first space only missing.")
    # rextendr::rust_source(system.file("extdata/lib.rs", package = "igsc"))
    # attr_ind <- rep(seq_along(gtf$attribute), lengths(stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)))
    # attr_col <- process_attr_col_rust(gtf$attribute) #igsc:::
  }

  gtf$index <- as.character(seq(1, nrow(gtf), 1))

  attr_col <- data.frame(attribute = attr_col[seq(1, length(attr_col), 2)],
                         value = attr_col[seq(2, length(attr_col), 2)],
                         index = as.character(attr_ind))

  if (!is.null(attr_rename)) {
    # remove previous
    attr_rename_exists <- attr_rename[which(attr_rename %in% attr_col$attribute)]
    attr_col <- dplyr::filter(attr_col, !attribute %in% attr_rename_exists)
    if (rename_replace) {
      # attr_col$attribute[attr_col$attribute %in% names(attr_rename)] <- attr_rename[attr_col$attribute[attr_col$attribute %in% names(attr_rename)]]
      attr_col <- dplyr::mutate(attr_col, attribute = dplyr::recode(attribute, !!!attr_rename))
    } else {
      attr_col <- dplyr::bind_rows(
        attr_col,
        attr_col |>
          dplyr::filter(attribute %in% names(attr_rename)) |>
          dplyr::mutate(attribute = dplyr::recode(attribute, !!!attr_rename))
      )
    }
  }
  if (is.null(attr_keep)) {
    attr_keep <- unique(attr_col$attribute)
  }
  attr_keep <- unique(c(attr_keep, unname(attr_rename))) # append attr_keep by attr_rename
  attr_keep <- attr_keep[which(attr_keep %in% unique(attr_col$attribute))]

  attr_col2 <- attr_col |>
    dplyr::filter(attribute %in% attr_keep) |>
    tidyr::pivot_wider(names_from = attribute, values_from = value, values_fn = ~paste(.x, collapse = "__"))

  if ("gene_name" %in% names(attr_col2) || !is.null(gene_name_force)) {
    if (!is.null(gene_name_force) && !gene_name_force %in% names(attr_col2)) {
      stop("gene_name_force column not found.")
    }
    if (!"gene_name" %in% names(attr_col2)) {
      attr_col2$gene_name <- NA
    }
    if (!is.null(gene_name_force)) {
      ## needed when not only exons are present and gene_name is NA in some cases of the same gene_id
      attr_col2 <- fill_gene_name_warn(attr_col2) # fill_gene_name_strict
      # needed when all are NA
      attr_col2$gene_name <- ifelse(is.na(attr_col2$gene_name), attr_col2[[gene_name_force]], attr_col2$gene_name)
    }
    attr_col2$gene_name <- ifelse(is.na(attr_col2$gene_name), NA, paste0(gene_name_prefix, attr_col2$gene_name))

    for (i in seq_along(gene_name_replace)) {
      attr_col2$gene_name <- stringr::str_replace_all(
        attr_col2$gene_name,
        names(gene_name_replace)[i],
        unname(gene_name_replace)[i]
      )
    }

    if (gene_name_fix) {
      attr_col2$gene_name <- gsub("[^A-Za-z0-9]+", "-", attr_col2$gene_name)
      attr_col2$gene_name <- gsub("^-|-$", "", attr_col2$gene_name)
    }


    if (fill_na) {
      attr_col2 <- attr_col2 |>
        dplyr::mutate(transcript_id = ifelse(is.na(transcript_id), gene_id, transcript_id)) |>
        dplyr::mutate(transcript_name = ifelse(is.na(transcript_name), gene_name, transcript_name))
      if ("exon_number" %in% names(attr_col2)) {
        if (any(grepl(";", gtf[["exon_number"]]))) {
          attr_col2 <- dplyr::mutate(attr_col2, exon_number = ifelse(is.na(exon_number), "0; exon_id 0", exon_number))
        } else {
          attr_col2 <- dplyr::mutate(attr_col2, exon_number = ifelse(is.na(exon_number), "0", exon_number))
        }
      }
    }

    # check for duplicates of pairs of gene_id, gene_name
    attr_col2 <- fix_duplicates(
      attr_col = attr_col2,
      verbose = verbose
    )

  }

  ## rotate the circular reference genome
  if (check_for_rotation) {
    df <- dplyr::left_join(attr_col2, gtf |> dplyr::select(-attribute), by = "index")
    cut <- pick_best_cut(df, genome_length = genome_length)
    if (!is.null(cut)) {
      df <- rotate_coords(df, cut = cut$cut_position, genome_length = genome_length)
      message("gtf coords have been rotated. use igsc:::rotate_genome_string with cut = ", cut$cut_position, " to rotate the genome. then save to fasta.")
      ## separate gtf cols and attr cols so that code below works
      gtf <- df |> dplyr::select(dplyr::all_of(names(gtf)[which(names(gtf) %in% names(df))]))
      attr_col2 <- df |> dplyr::select(dplyr::all_of(names(df)[which(!names(df) %in% names(gtf))]), index)
    }
  }

  ## highly advisable for when aggregate_overlapping_exon_ranges is done
  # if genes are not aggregated, then single exons from different genes may be joined
  # but then the gtf becomes invalid
  if (aggregate_exons) {
    df <- dplyr::left_join(attr_col2, gtf |>
                             dplyr::select(-dplyr::any_of("attribute")), by = "index")
    df <- split(df, df$feature)
    df[["exon"]] <- split(df[["exon"]], df[["exon"]]$transcript_id)
    ## min max of start and end and summary of other cols
    df[["exon"]] <- purrr::map(df[["exon"]], function(x) {
      minstart <- min(x$start)
      maxend <- max(x$end)
      allindex <- paste(unique(x$index), collapse = "---")
      x |>
        dplyr::select(-dplyr::any_of("exon_number")) |>
        dplyr::summarise(dplyr::across(-c(start, end), most_frequent)) |>
        dplyr::mutate(start = minstart, end = maxend, index = allindex, exon_number = 1)
    })
    df[["exon"]] <- dplyr::bind_rows(df[["exon"]])
    df <- dplyr::bind_rows(df)

    ## separate gtf cols and attr cols so that code below works
    gtf <- df |> dplyr::select(dplyr::all_of(names(gtf)[which(names(gtf) %in% names(df))]))
    attr_col2 <- df |> dplyr::select(dplyr::all_of(names(df)[which(!names(df) %in% names(gtf))]), index)
  }

  ### this is specifically for viral genomes with overlapping ranges
  if (aggregate_overlapping_exon_ranges) {
    igsc:::.ensure_package("igraph")

    # attr_col2 then needs to be filtered by whats left in df
    df <- dplyr::left_join(attr_col2, gtf |>
                             dplyr::select(-dplyr::any_of("attribute")), by = "index")

    ## split by feature; only handle exons
    df <- split(df, df$feature)
    ## split by seqname
    df[["exon"]] <- split(df[["exon"]], df[["exon"]]$seqname)
    df[["exon"]] <- purrr::map(df[["exon"]], function(x) {
      # just run two rounds each
      x <- aggregate_overlapping_exon_ranges_fun(x, strand = "same")
      x <- aggregate_overlapping_exon_ranges_fun(x, strand = "opposing")
      x <- aggregate_overlapping_exon_ranges_fun(x, strand = "same")
      x <- aggregate_overlapping_exon_ranges_fun(x, strand = "opposing")
      return(x)
    })
    df[["exon"]] <- dplyr::bind_rows(df[["exon"]])
    df <- dplyr::bind_rows(df)

    ## separate gtf cols and attr cols so that code below works
    gtf <- df |> dplyr::select(dplyr::all_of(names(gtf)[which(names(gtf) %in% names(df))]))
    attr_col2 <- df |> dplyr::select(dplyr::all_of(names(df)[which(!names(df) %in% names(gtf))]), index)
  }

  # if (attr_as == "kv") {
  #   # optimized for speed
  #   attr_col2 <- attr_col2 |> tidyr::nest(data = index)
  #   for (i in names(attr_col2)[-ncol(attr_col2)]) {
  #     attr_col2[[i]] <- paste0(i, " \"", attr_col2[[i]], "\"; ")
  #   }
  #   attr_col2$attribute <- apply(attr_col2[,-ncol(attr_col2)], 1, paste, collapse = "")
  #   attr_col2 <- attr_col2 |>
  #     dplyr::select(data, attribute) |>
  #     tidyr::unnest(cols = c(data))
  # }

  ## when did i need this: attr_col2 <- attr_col2 |> tidyr::nest(data = index) ??

  gtf <- dplyr::left_join(dplyr::select(gtf, -dplyr::any_of("attribute")), attr_col2, by = "index")

  if (rm_entries_wo_matching_exon) {
    exon_entry_transcript_id <- gtf |>
      dplyr::filter(feature == "exon") |>
      dplyr::distinct(transcript_id) |>
      dplyr::pull(transcript_id)
    transcript_id_gene_names <- gtf |>
      dplyr::filter(transcript_id %in% exon_entry_transcript_id) |>
      dplyr::distinct(gene_name) |>
      dplyr::pull(gene_name)
    gtf <- dplyr::bind_rows(dplyr::filter(gtf, transcript_id %in% exon_entry_transcript_id),
                            gtf |>
                              dplyr::filter(gene_name %in% transcript_id_gene_names) |>
                              dplyr::filter(feature == "gene"))
  }


  if (attr_as == "kv") {
    gtf <- make_kv_attr_col(gtf, keep_index_col = T, verbose = verbose)
  }

  if (rm_index) {
    gtf <- dplyr::select(gtf, -index)
  }

  if ("exon_number" %in% names(gtf)) {
    if (any(grepl(";", gtf[["exon_number"]]))) {
      # assume 1; exon_id ENSE00003967718
      gtf <- gtf |>
        tidyr::separate(exon_number, into = c("exon_number", "exon_id"), sep = ";") |>
        dplyr::mutate(exon_id = trimws(gsub("exon_id", "", exon_id)))
    }
    gtf[["exon_number"]] <- as.numeric(gtf[["exon_number"]])
  }

  # filter for requested only
  attr_col <- attr_col |>
    dplyr::filter(attribute %in% attr_keep) |>
    tidyr::pivot_wider(names_from = attribute, values_from = value, values_fn = ~paste(.x, collapse = "__"))

  return(list(gtf = gtf, attr = attr_col))
}


get_bounds <- function(x, file_path) {
  # here we get the first and last line of a seqname to read with vroom
  # actually though, rg and grep do return the full lines already, not only linenumbers
  # but, so what

  if (grepl("\\.gz$", file_path)) {
    file_path <- brathering::ungunzip(
      file_path,
      out_dir = tempdir(),
      out_file = tools::file_path_sans_ext(basename(file_path))
    )
    message("unpacking file to: ", file_path)
  }

  out <- tryCatch(
    {
      # use ripgrep if possible
      cmd <- paste0("rg '^", x, "\t' -n ", file_path, " | cut -d: -f1")
      system(cmd, intern = T)
    },
    error = function(err) {
      # else use grep which is slower but more common
      cmd <- paste0("grep '^", x, "\t' -n ", file_path, " | cut -d: -f1")
      system(cmd, intern = T)
    }
  )
  if (length(out) == 0) {
    stop("seqname not found in gtf file.")
  }
  out <- as.numeric(out[c(1, length(out))])
  return(list(bounds = out, file = file_path))
}

vroom_gtf <- function(x, file_path, col_names, unpack_fun) {
  bounds <- get_bounds(x, file_path)
  file <- bounds[["file"]]
  bounds <- bounds[["bounds"]]

  y <- vroom::vroom(file = file,
                    col_names = col_names,
                    skip = bounds[1] - 1,
                    n_max = bounds[2] - bounds[1] + 1,
                    comment = "#",
                    progress = F,
                    show_col_types = F)
  return(y)
}



fill_gene_name_warn <- function(df) {
  df |>
    dplyr::group_by(gene_id) |>
    dplyr::group_modify(~ {
      vals <- stats::na.omit(.x$gene_name)

      if (length(vals) == 0) {
        return(.x)
      }

      tab <- sort(table(vals), decreasing = TRUE)
      chosen <- names(tab)[1]

      if (length(tab) > 1) {
        message(
          sprintf(
            "gene_id %s has multiple gene_name values (%s); using most frequent: %s",
            .y$gene_id,
            paste(names(tab), collapse = ", "),
            chosen
          )
        )
      }

      .x$gene_name[is.na(.x$gene_name)] <- chosen
      .x
    }) |>
    dplyr::ungroup()
}

fill_gene_name_strict <- function(df) {
  df |>
    dplyr::group_by(gene_id) |>
    dplyr::group_modify(~ {
      vals <- unique(stats::na.omit(.x$gene_name))

      if (length(vals) > 1) {
        stop(
          sprintf(
            "Conflicting gene_name values for gene_id %s: %s",
            .y$gene_id,
            paste(vals, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      fill <- if (length(vals) == 1) vals else NA_character_
      .x$gene_name[is.na(.x$gene_name)] <- fill
      .x
    }) |>
    dplyr::ungroup()
}


aggregate_overlapping_exon_ranges_fun <- function(df, strand = c("same", "opposing")) {

  strand <- rlang::arg_match(strand)

  ovlp_groups <- get_range_overlap_groups(df = df, strand = strand)
  if (is.null(ovlp_groups)) {
    return(df)
  }
  ## collapse / aggregate
  df2 <- dplyr::left_join(df, ovlp_groups, by = "index") |>
    ## add unique random integers to NA in group
    dplyr::mutate(
      group = {
        na_n <- sum(is.na(group))
        used <- group[!is.na(group)]
        pool <- setdiff(seq_len(max(used) + na_n * 2), used)
        group[is.na(group)] <- sample(pool, na_n)
        group
      }
    ) |>
    ## collapse all columns by group
    dplyr::group_by(group) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ paste(unique(.), collapse = "---")
      ),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    ## the fix selected columns
    dplyr::mutate(start = as.numeric(min(strsplit(start, "---")[[1]])),
                  end = as.numeric(max(strsplit(end, "---")[[1]])),
                  frame = strsplit(frame, "---")[[1]][1]) |>
    dplyr::ungroup() |>
    dplyr::select(-c(group))

  if (strand == "opposing") {
    # just set to plus strand because STAR needs a choice
    df2$strand <- ifelse(grepl("---", df2$strand), "+", df2$strand)
  }
  return(df2)
}

get_range_overlap_groups <- function(df, strand = c("same", "opposing")) {

  strand <- rlang::arg_match(strand)
  same_strand <- strand == "same"
  n <- nrow(df)
  overlap <- outer(
    seq_len(n), seq_len(n),
    Vectorize(function(i, j) {
      df$start[i] <= df$end[j] &&
        df$start[j] <= df$end[i] &&
        (df$strand[i] == df$strand[j]) == same_strand
    })
  )

  diag(overlap) <- FALSE
  idx <- which(overlap, arr.ind = TRUE)

  graphdf <- data.frame(
    gene1 = df$index[idx[,1]],
    gene2 = df$index[idx[,2]],
    stringsAsFactors = F
  )
  if (nrow(graphdf) == 0) {
    return(NULL)
  }
  g <- igraph::graph_from_data_frame(graphdf, directed = FALSE)

  groups <- igraph::components(g)$membership
  groups <- utils::stack(groups)
  names(groups) <- c("group", "index")
  groups$index <- as.character(groups$index)
  return(groups)
}

detect_wrap_and_cut <- function(df, genome_length = NULL, verbose = T) {

  ## better version of pick_best_cut

  df <- as.data.frame(df)

  if (!"exon_number" %in% names(df)) {
    stop("exon_number not in df.")
  }
  # Use only CDS with exon numbers
  cds <- df[df$feature == "CDS" & !is.na(df$exon_number), ]

  if (nrow(cds) == 0) {
    # try exons
    cds <- df[df$feature == "exon" & !is.na(df$exon_number), ] |>
      dplyr::mutate(feature = "CDS")
    if (nrow(cds) == 0) {
      stop("No CDS entries with exon_number found.")
    }
  }

  # Check monotonicity per gene (strand-aware)
  wrap_flags <- tapply(seq_len(nrow(cds)), cds$gene_id, function(idx) {
    g <- cds[idx, ]

    # order by exon_number
    g <- g[order(g$exon_number), ]

    starts <- g$start
    strand <- unique(g$strand)

    # sanity check
    if (length(strand) != 1) {
      return(TRUE)  # inconsistent strand - treat as problematic
    }

    if (strand == "+") {
      # should increase
      any(diff(starts) < 0)
    } else if (strand == "-") {
      # should decrease
      any(diff(starts) > 0)
    } else {
      TRUE  # unknown strand - treat as problematic
    }
  })

  wrapping_genes <- names(wrap_flags)[wrap_flags]

  if (length(wrapping_genes) == 0) {
    if (verbose) {
      message("All genes follow expected strand-specific monotonicity. No rotation needed.")
    }
    return(NULL)
  }
  if (verbose) {
    message("Wrap-around genes detected: ", paste(wrapping_genes, collapse = ", "))
  }
  # Infer genome length if needed
  if (is.null(genome_length)) {
    genome_length <- max(df$end, na.rm = TRUE)
  }

  cut_info <- pick_best_cut_from_positions(df, genome_length = genome_length)

  # overwrite wrapping genes with biologically detected ones
  cut_info$wrapping_genes <- wrapping_genes

  return(cut_info)
}

pick_best_cut_from_positions <- function(df, genome_length = NULL) {

  df <- as.data.frame(df)

  # Infer genome length if not provided
  if (is.null(genome_length)) {
    genome_length <- max(df$end, na.rm = TRUE)
  }

  L <- genome_length

  # Collapse to gene spans
  gene_spans <- stats::aggregate(
    cbind(start, end) ~ gene_id,
    df,
    function(x) c(min = min(x), max = max(x))
  )

  gene_spans$start <- gene_spans$start[, "min"]
  gene_spans$end   <- gene_spans$end[, "max"]

  # Order genes by start
  gene_spans <- gene_spans[order(gene_spans$start), ]
  n <- nrow(gene_spans)

  if (n < 2) {
    stop("Not enough genes to compute gaps.")
  }

  # Compute gaps
  gaps <- numeric(n)
  gap_starts <- numeric(n)
  gap_ends <- numeric(n)

  # Linear gaps
  for (i in seq_len(n - 1)) {
    gaps[i] <- gene_spans$start[i + 1] - gene_spans$end[i] - 1
    gap_starts[i] <- gene_spans$end[i] + 1
    gap_ends[i] <- gene_spans$start[i + 1] - 1
  }

  # Circular gap (last -> first)
  gaps[n] <- (gene_spans$start[1] + L) - gene_spans$end[n] - 1
  gap_starts[n] <- gene_spans$end[n] + 1
  gap_ends[n] <- gene_spans$start[1] - 1

  if (gap_ends[n] <= 0) {
    gap_ends[n] <- gap_ends[n] + L
  }

  # Find largest gap
  best <- which.max(gaps)

  if (gaps[best] <= 0) {
    stop("No intergenic gap found - genome fully covered.")
  }

  # Midpoint of best gap
  cut <- floor((gap_starts[best] + gap_ends[best]) / 2)
browser()
  # Wrap if needed
  if (cut > L) {
    cut <- cut - L
  }

  return(list(
    cut_position = cut,
    gap_size = gaps[best],
    gap_start = gap_starts[best],
    gap_end = gap_ends[best]
  ))
}

pick_best_cut <- function(df, genome_length = NULL, wrap_threshold = 0.5) {

  df <- as.data.frame(df)

  if (is.null(genome_length)) {
    genome_length <- max(df$end, na.rm = TRUE)
  }

  L <- genome_length

  # Collapse to gene spans
  gene_spans <- stats::aggregate(
    cbind(start, end) ~ gene_id,
    df,
    function(x) c(min = min(x), max = max(x))
  )

  gene_spans$start <- gene_spans$start[, "min"]
  gene_spans$end   <- gene_spans$end[, "max"]

  # Detect wrap-around genes
  span_width <- gene_spans$end - gene_spans$start

  wrapping_genes <- gene_spans$gene_id[span_width > (wrap_threshold * L)]

  if (length(wrapping_genes) == 0) {
    message("No wrap-around gene detected. Rotation not needed.")
    return(NULL)
  }

  message("Wrap-around gene(s) detected: ", paste(wrapping_genes, collapse = ", "))

  # Compute largest intergenic gap (same as before)
  gene_spans <- gene_spans[order(gene_spans$start), ]
  n <- nrow(gene_spans)

  gaps <- numeric(n)
  gap_starts <- numeric(n)
  gap_ends <- numeric(n)

  for (i in seq_len(n - 1)) {
    gaps[i] <- gene_spans$start[i + 1] - gene_spans$end[i] - 1
    gap_starts[i] <- gene_spans$end[i] + 1
    gap_ends[i] <- gene_spans$start[i + 1] - 1
  }

  # circular gap
  gaps[n] <- (gene_spans$start[1] + L) - gene_spans$end[n] - 1
  gap_starts[n] <- gene_spans$end[n] + 1
  gap_ends[n] <- gene_spans$start[1] - 1
  if (gap_ends[n] <= 0) gap_ends[n] <- gap_ends[n] + L

  best <- which.max(gaps)

  if (gaps[best] <= 0) {
    stop("No intergenic gap found - genome fully covered.")
  }

  cut <- floor((gap_starts[best] + gap_ends[best]) / 2)
  if (cut > L) cut <- cut - L

  return(list(
    cut_position = cut,
    gap_size = gaps[best],
    wrapping_genes = wrapping_genes
  ))
}

rotate_coords <- function(df, cut, genome_length = NULL) {

  df <- as.data.frame(df)

  if (is.null(genome_length)) {
    genome_length <- max(df$end, na.rm = TRUE)
  }

  L <- genome_length

  # shift coordinates
  shift <- function(x) {
    y <- x - cut + 1
    y[y <= 0] <- y[y <= 0] + L
    y
  }

  df$start <- shift(df$start)
  df$end   <- shift(df$end)

  # split features that cross new origin
  wrap <- df$start > df$end

  if (any(wrap)) {

    df_a <- df[wrap, ]
    df_b <- df[wrap, ]

    df_a$end   <- L
    df_b$start <- 1

    df <- rbind(df[!wrap, ], df_a, df_b)
  }

  rownames(df) <- NULL
  df <- fix_duplicate_rows(df)
  return(df)
}

fix_duplicate_rows <- function(df) {

  ## see example: EBV genome, gene of HHV4tp2_gp01
  df <- df |> dplyr::mutate(row = dplyr::row_number())
  # same entries but multi start/end due to crossing the artificial cut
  dfdup <- df |>
    dplyr::group_by(dplyr::across(-c(start, end, row))) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::summarise(
      start = min(start),
      end = max(end),
      row = min(row),
      .groups = "drop")
  dfunique <- df |>
    dplyr::group_by(dplyr::across(-c(start, end, row))) |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::ungroup()
  df <- dplyr::bind_rows(dfunique, dfdup) |>
    dplyr::arrange(row) |>
    dplyr::select(-row)
  return(df)
}

rotate_genome_string <- function(genome, cut) {

  L <- nchar(genome)

  if (cut < 1 || cut > L) {
    stop("cut must be between 1 and genome length")
  }

  if (cut == 1) {
    return(genome)  # no rotation needed
  }

  part1 <- substr(genome, cut, L)
  part2 <- substr(genome, 1, cut - 1)

  return(stats::setNames(paste0(part1, part2), names(genome)))
}


most_frequent <- function(x) {
  ux <- stats::na.omit(x)
  if (length(ux) == 0) return(NA)
  tab <- table(ux)
  names(tab)[which.max(tab)]
}



#' Make kv pairs in attribute col for saving proper gtf df
#'
#' From a gtf df with processed attr col, make a attribute column with
#' key-value pairs
#'
#' @param gtf_df subsetted gtf df with processed attribute column
#' @param keep_index_col keep index column away from attributes
#' @param verbose print messages?
#'
#' @returns data frame
#' @export
#'
#' @examples
#' \dontrun{
#' gtf <- make_kv_attr_col(gtf)
#' write_gtf(gtf)
#' }
make_kv_attr_col <- function(gtf_df,
                             keep_index_col = F,
                             verbose = T) {
  gtf_cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame")
  if (keep_index_col && "index" %in% names(gtf_df)) {
    gtf_cols <- c(gtf_cols, "index")
  }
  if (any(!gtf_cols %in% names(gtf_df))) {
    stop(paste(gtf_cols, collapse = ","), " columns are needed.")
  }
  attr_cols <- names(gtf_df)[!names(gtf_df) %in% gtf_cols]
  if (!length(attr_cols)) {
    if (verbose) {
      message("no attr cols found")
    }
    return(gtf_df)
  } else {
    if (verbose) {
      message("attr cols: ", paste(attr_cols, collapse = ","))
    }
  }

  gtf <- gtf_df[gtf_cols]
  attr_col <- gtf_df[attr_cols]

  # optimized for speed, maybe
  for (i in names(attr_col)) {
    attr_col[[i]] <- paste0(i, " \"", attr_col[[i]], "\"; ")
  }
  gtf$attribute <- apply(attr_col, 1, paste, collapse = "")

  return(gtf)


  # if ("index" %in% names(gtf) && anyDuplicated(gtf$index)) {
  #
  #   gtf_cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "index")
  #   gtf <- out$gtf_df[gtf_cols]
  #   attr_col <- out$gtf_df[which(!names(out$gtf_df) %in% gtf_cols)]
  #
  #   # optimized for speed, maybe
  #   attr_col <- tidyr::nest(attr_col, data = index)
  #   for (i in names(attr_col)[which(names(attr_col) != "data")]) {
  #     attr_col[[i]] <- paste0(i, " \"", attr_col[[i]], "\"; ")
  #   }
  #   attr_col$attribute <- apply(attr_col[,which(names(attr_col) != "data")], 1, paste, collapse = "")
  #   attr_col <- attr_col |>
  #     dplyr::select(data, attribute) |>
  #     tidyr::unnest(cols = c(data))
  #   gtf <- dplyr::left_join(dplyr::select(gtf, -dplyr::any_of("attribute")), attr_col, by = "index")
  #
  # } else {
  #
  #   gtf_cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame")
  #   gtf <- out$gtf_df[gtf_cols]
  #   attr_col <- out$gtf_df[which(!names(out$gtf_df) %in% gtf_cols)] |>
  #     dplyr::select(-dplyr::any_of("index"))
  #
  #   # optimized for speed, maybe
  #   for (i in names(attr_col)) {
  #     attr_col[[i]] <- paste0(i, " \"", attr_col[[i]], "\"; ")
  #   }
  #   gtf$attribute <- apply(attr_col, 1, paste, collapse = "")
  # }
  #
}

fix_duplicates <- function(attr_col, verbose = T) {

  ## first gene_id vs gene_name
  dups <- attr_col |>
    dplyr::distinct(gene_id, gene_name) |>
    dplyr::add_count(gene_name, name = "n_gene_name") |>
    dplyr::add_count(gene_id, name = "n_gene_id")
  dups1 <- dups |> dplyr::filter(n_gene_name>1) |> dplyr::arrange(gene_name)

  if (nrow(dups1) > 1) {
    if (verbose) {
      message("duplicate gene names made unique:")
      print(dups1)
    }

    gene_name_map <- attr_col |>
      dplyr::distinct(gene_id, gene_name) |>
      dplyr::group_by(gene_name) |>
      dplyr::mutate(
        gene_name_unique = if (dplyr::n() > 1) {
          paste0(gene_name, "--", dplyr::row_number()) # use two dashes to separate from natural numeric extension of some genes
        } else {
          gene_name
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::add_count(gene_id)

    attr_col <- attr_col |>
      dplyr::left_join(
        gene_name_map |> dplyr::distinct(gene_id, gene_name_unique),
        by = "gene_id"
      ) |>
      dplyr::mutate(gene_name = gene_name_unique) |>
      dplyr::select(-gene_name_unique)
    message("UWAGA when many-to-many relationship is shown.")
  }

  dups2 <- dups |> dplyr::filter(n_gene_id>1) |> dplyr::arrange(gene_id)

  if (nrow(dups2) > 1) {
    if (verbose) {
      message("duplicate gene ids made unique:")
      print(dups2)
    }
    gene_id_map <- attr_col |>
      dplyr::distinct(gene_id, gene_name) |>
      dplyr::group_by(gene_id) |>
      dplyr::mutate(
        gene_id_unique = if (dplyr::n() > 1) {
          paste0(gene_id, "--", dplyr::row_number()) # use two dashes to separate from natural numeric extension of some genes
        } else {
          gene_id
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::add_count(gene_name)

    attr_col <- attr_col |>
      dplyr::left_join(
        gene_id_map |> dplyr::distinct(gene_name, gene_id_unique),
        by = "gene_name"
      ) |>
      dplyr::mutate(gene_id = gene_id_unique) |>
      dplyr::select(-gene_id_unique)

    message("UWAGA when many-to-many relationship is shown.")
  }

  ## second gene_id vs transcript_id
  if ("transcript_id" %in% names(attr_col)) {

    dups <- attr_col |>
      dplyr::distinct(gene_id, transcript_id) |>
      dplyr::add_count(transcript_id, name = "n_transcript_id") |>
      dplyr::add_count(gene_id, name = "n_gene_id")

    dups3 <- dups |> dplyr::filter(n_transcript_id>1) |> dplyr::arrange(transcript_id)

    if (nrow(dups3) > 1) {
      if (verbose) {
        message("duplicate transcript ids made unique:")
        print(dups3)
      }
      transcript_id_map <- attr_col |>
        dplyr::distinct(gene_id, transcript_id) |>
        dplyr::group_by(transcript_id) |>
        dplyr::mutate(
          transcript_id_unique = if (dplyr::n() > 1) {
            paste0(transcript_id, "--", dplyr::row_number()) # use two dashes to separate from natural numeric extension of some genes
          } else {
            transcript_id
          }
        ) |>
        dplyr::ungroup() |>
        dplyr::add_count(gene_id)

      attr_col <- attr_col |>
        dplyr::left_join(
          transcript_id_map |> dplyr::distinct(gene_id, transcript_id, transcript_id_unique),
          dplyr::join_by(gene_id, transcript_id)
        ) |>
        dplyr::mutate(transcript_id = transcript_id_unique) |>
        dplyr::select(-transcript_id_unique)
    }

  }

  return(attr_col)
}
