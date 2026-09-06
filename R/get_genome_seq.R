#' Retrieve a genomic subsequence
#'
#' Extract a genomic interval from a multi-record FASTA file, from a directory
#' containing one `fst` file per sequence, or from the UCSC Genome Browser API.
#' Coordinates supplied to this function are one-based and inclusive.
#'
#' The data source is selected in this order: `fasta_file`, `fst_folder`, then
#' `ucsc_api`. In API mode, the start coordinate is converted to UCSC's
#' zero-based, half-open convention before the request is made. Local source
#' names are matched to `chromosome`; if no exact match exists, the name with
#' the smallest edit distance is selected and reported in a message.
#'
#' `query_string` can replace the separate `chromosome`, `start`, and `end`
#' arguments. Accepted forms include `"chr2:1000-2000"`,
#' `"chr2:86,784,610-86,790,913"`,
#' `"chrom=chr2;start=1000;end=2000"`, and `"2,1000,2000"`. The parsed
#' values are printed before retrieval.
#'
#' @param fasta_file Optional path to a multi-record FASTA file. Sequence bounds
#'   are located with [get_fasta_seq_bounds()] and the selected record is read
#'   with [read_fasta()]. When supplied, this source takes precedence over all
#'   folder and API sources.
#' @param fst_folder Optional path to a directory containing per-sequence `.fst`
#'   files. Each filename without its extension is treated as a sequence name,
#'   and the first column must contain sequence characters in positional order.
#' @param fagz_folder Reserved for a directory of per-sequence `.fa.gz` files.
#'   This source is not reachable in the current dispatch implementation; use
#'   `fasta_file` or `fst_folder` instead.
#' @param ucsc_api Base URL for a UCSC sequence endpoint, including the genome
#'   and a trailing `"&"`. Used only when both `fasta_file` and `fst_folder` are
#'   `NULL`. The default requests the human `hg38` assembly through the external
#'   `curl` command.
#' @param query_string Optional character string encoding the sequence name,
#'   start, and end coordinates. When supplied, its parsed values override
#'   `chromosome`, `start`, and `end`.
#' @param chromosome Sequence or chromosome name. Required unless supplied in
#'   `query_string`. API sequence names generally begin with `"chr"`.
#' @param start One-based inclusive start coordinate, or a value coercible to
#'   numeric. Required unless supplied in `query_string`.
#' @param end One-based inclusive end coordinate, or a value coercible to
#'   numeric. Required unless supplied in `query_string`.
#' @param ucsc_toupper Logical; convert API results to uppercase. This does not
#'   alter sequences read from local FASTA or `fst` files.
#'
#' @return A single character string containing the requested genomic sequence.
#' @export
#'
#' @examples
#' \dontrun{
#' # Read from a local multi-record FASTA file.
#' seq <- get_genome_seq(
#'   fasta_file = "reference/genome.fa",
#'   query_string = "chr2:1000-2000"
#' )
#'
#' # Read rows from a per-chromosome fst file.
#' seq <- get_genome_seq(
#'   fst_folder = "reference/fst",
#'   chromosome = "chr2",
#'   start = 1000,
#'   end = 2000
#' )
#'
#' # With no local source, query the default UCSC hg38 endpoint.
#' seq <- get_genome_seq(query_string = "chr2:1000-2000")
#' }
get_genome_seq <- function(fasta_file = NULL,
                           fst_folder = NULL,
                           fagz_folder = NULL,
                           ucsc_api = "http://api.genome.ucsc.edu/getData/sequence?genome=hg38&",
                           query_string = NULL,
                           chromosome = NULL,
                           start = NULL,
                           end = NULL,
                           ucsc_toupper = T) {

  # query_string <- "chr2:86784610-86790913"
  # query_string <- "chr2:86,784,610-86,790,913"
  # query_string <- "chrom=chr1;start=1000000;end=1000100"
  # query_string <- "2,1000,30000"

  if (!is.null(query_string)) {
    query_string <- tolower(query_string)
    query_string <- gsub("chrom=", "", query_string)
    query_string <- gsub(" ", "", query_string)
    if (nchar(query_string) - nchar(gsub(",", "", query_string)) > 2) {
      # we want to allow the query '2,5000,8000'
      # but when there are 1000's seperators, we want to remove them 'chr2:86,784,610-86,790,913'
      query_string <- gsub(",", "", query_string)
    }
    #query_string <- strsplit(query_string, "[^0-9]+")[[1]]
    query_string <- strsplit(query_string, "[^a-z0-9]+")[[1]]
    query_string <- query_string[which(query_string != "")]
    query_string <- c(query_string[1], query_string[which(!is.na(suppressWarnings(as.numeric(query_string))))])

    cat(query_string)
    cat("\n")
    if (length(query_string) != 3) {
      stop("query_string could not be resolved.")
    }
    chromosome <- query_string[1]
    start <- query_string[2]
    end <- query_string[3]
  }


  ## sequence name has to start with "chr"
  if (is.null(fasta_file) && is.null(fst_folder)) {
    # ucsc is 0-based, make start start-1 to match other results
    start <- as.character(as.numeric(start) - 1)
    ucsc_query <- paste0("curl '", ucsc_api, "chrom=", chromosome, "&start=", start, "&end=", end, "'")
    out <- system(ucsc_query, intern = T)[2]
    out <- strsplit(out, " ")[[1]]
    out <- out[length(out)]
    out <- stringr::str_sub(out, 2, -3)
    if (ucsc_toupper) {
      out <- toupper(out)
    }
    #out <- stringi::stri_replace_all(out, replacement = "", fixed = '"')
    return(out)
  }

  if (!is.null(fasta_file)) {
    seq_bounds <- get_fasta_seq_bounds(fasta_file)
    if (!any(chromosome == seq_bounds$name)) {
      chr_name_before <- chromosome
      chromosome <- seq_bounds$name[which.min(utils::adist(chromosome, seq_bounds$name)[1,])]
      message("No exact match for sequence name. Closest match: ", chr_name_before, " --> ", chromosome)
    }
    ind <- which(seq_bounds$name == chromosome)
    refseq <- read_fasta(fasta_file,
                         start_line = seq_bounds[ind, "start_line"],
                         end_line = seq_bounds[ind, "end_line"])
    out <- unname(substr(refseq, start, end))
    return(out)
  }

  if (!is.null(fst_folder)) {
    fst_files <- list.files(fst_folder, pattern = "\\.fst$", full.names = T, ignore.case = T)
    fst_names <- gsub("\\.fst$", "", basename(fst_files))
    if (!any(chromosome == fst_names)) {
      chr_name_before <- chromosome
      chromosome <- fst_names[which.min(utils::adist(chromosome, fst_names)[1,])]
      message("No exact match for sequence name. Closest match: ", chr_name_before, " --> ", chromosome)
    }
    ind <- which(fst_names == chromosome)
    out <- fst::read_fst(path = fst_files[ind], from = as.numeric(start), to = as.numeric(end))
    out <- paste(out[,1,drop=T], collapse = "")
    return(out)
  }

  if (!is.null(fagz_folder)) {
    fagz_files <- list.files(fst_folder, pattern = "\\.fa\\.gz$", full.names = T, ignore.case = T)
    if (length(fagz_files) > 0) {
      fagz_names <- gsub("\\.fa\\.gz$", "", basename(fagz_names))
      if (!any(chromosome == fagz_names)) {
        chr_name_before <- chromosome
        chromosome <- fagz_names[which.min(utils::adist(chromosome, fagz_names)[1,])]
        message("No exact match for sequence name. Closest match: ", chr_name_before, " --> ", chromosome)
      }
      ind <- which(fagz_names == chromosome)
      out <- read_fasta(fagz_files[ind])
      out <- stringr::str_sub(out, start = as.numeric(start), end = as.numeric(end))
      return(out)
    }
  }
}

