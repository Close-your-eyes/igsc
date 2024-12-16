#' Title
#'
#' @param fasta_file
#' @param fst_folder
#' @param ucsc_api
#' @param query_string
#' @param chromosome
#' @param start
#' @param end
#' @param ucsc_toupper
#'
#' @return
#' @export
#'
#' @examples
get_genome_seq <- function(fasta_file = NULL,
                           fst_folder = NULL,
                           ucsc_api = "http://api.genome.ucsc.edu/getData/sequence?genome=hg38&",
                           query_string = NULL,
                           chromosome = NULL,
                           start = NULL,
                           end = NULL,
                           ucsc_toupper = T) {

  if (!is.null(query_string)) {
    # query_string <- "chr2:86784610-86790913"
    # query_string <- "chrom=chr1;start=1000000;end=1000100"
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
      chromosome <- seq_bounds$name[which.min(adist(chromosome, seq_bounds$name)[1,])]
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
      chromosome <- fst_names[which.min(adist(chromosome, fst_names)[1,])]
      message("No exact match for sequence name. Closest match: ", chr_name_before, " --> ", chromosome)
    }
    ind <- which(fst_names == chromosome)
    out <- fst::read_fst(path = fst_files[ind], from = as.numeric(start), to = as.numeric(end))
    out <- paste(out[,1,drop=T], collapse = "")
    return(out)
  }
}
