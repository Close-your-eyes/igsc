#' Read and process a GTF (Gene Transfer Format) file
#'
#' Since GTF files come in an unhandy format, this function may help
#' to easily read them into memory. As a whole or partly by providing
#' seqnames and/or features. Processing the attribute column is computationally
#' costly.
#'
#' See https://www.ensembl.org/info/website/upload/gff.html?redirect=no for
#' explanation of GTF file format.
#' GTF files and genomic fasta files may be downloaded here https://www.ncbi.nlm.nih.gov/datasets/genome/
#' or here https://www.ensembl.org/index.html
#'
#' @param file_path path to the file; file may be gunzipped (ending with .gz)
#' @param gtf data frame from reading gtf file with vroom or similar.
#' useful to pass a subset of rows only.
#' @param attr_col_as_list have the attributes column return as named list (TRUE) or
#' separated by names and values into two columns
#' @param seqnames seqnames to filter the gtf file for; will decrease computation time
#' required for processing the attribute column
#' @param features features to filter the gtf file for; will decrease computation time
#' required for processing the attribute column
#' @param process_attr_col convert the attribute column into separate columns
#' @param attr_keep which attributes to keep from attribute column
#' @param col_names column names to assign to the gtf data frame;
#' changing seqname, feature or attribute will break this function;
#' better leave col_names as it is
#' @param use_fun which function to use for processing the attr_col; rcpp is
#' fastest currently
#' @return a list with (i) entries of the GTF file including the attribute
#' column as list and some attributes as separate columns and
#' (ii) the attributes as long data frame
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#'
#' gtf <- read_gtf(your_path, attr_col_as_list = F)
#' # when attr_col_as_list = F attributes are split into names and values columns
#' # this is how to expand the attributes names and values columns
#' # tidyr unnest over two columns matches the list indices
#' gtf2 <-
#' gtf[["gtf"]] %>%
#' dplyr::mutate(attribute_names = I(strsplit(attribute_names, ",")),
#'                attribute_values = I(strsplit(attribute_values, ","))) %>%
#'                tidyr::unnest(cols = c(attribute_names, attribute_values))
#'
#' # this is how to make a named list from separate names and values columns
#' # this is returned when attr_col_as_list = T
#' gtf2 <-
#' gtf[["gtf"]] %>%
#' dplyr::mutate(attribute_names = I(strsplit(attribute_names, ",")),
#' attribute_values = I(strsplit(attribute_values, ","))) %>%
#' dplyr::rowwise() %>%
#' dplyr::mutate(attr = I(list(setNames(attribute_names, attribute_values))))
#'
#' # you may want to use the fst package to write the data frames to disk
#' # this allows quick reading and random access
#' fst::write_fst(x = gtf,
#' path = "your_path/your_filename_gtf.fst",
#' compress = 10)
#' fst::write_fst(x = attr_col,
#' path = "your_path/your_filename_attr.fst",
#' compress = 10)
#' }
read_gtf <- function(file_path,
                     gtf,
                     attr_col_as_list = F,
                     seqnames = NULL,
                     features = NULL,
                     process_attr_col = T,
                     col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
                     attr_keep = c("gene_id", "gene_name", "transcript_id", "transcript_name", "exon_number"),
                     use_fun = c("rcpp","rust", "r")) {

  if (missing(file_path) && missing(gtf)) {
    stop("Provide file_path or gtf data frame.")
  }

  if (!missing(file_path) && !missing(gtf)) {
    message("file_path and gtf provided. Will ignore file_path and work with the gtf data frame provided.")
  }

  if (grepl("\\.gz$", file_path)) {
    unpack_fun <- gzfile
  } else {
    unpack_fun <- function(description) {description}
  }

  use_fun <- match.arg(use_fun, c("rcpp","rust", "r"))

  if (missing(gtf)) {
    if (!is.null(seqnames)) {
      gtf <- do.call(rbind, lapply(seqnames, vroom_gtf, file_path = file_path, col_names = col_names, unpack_fun = unpack_fun))
    } else {
      gtf <- vroom::vroom(file = do.call(unpack_fun, args = list(description = file_path)),
                          col_names = col_names,
                          comment = "#",
                          show_col_types = F)
    }
  } else {
    if (!"attribute" %in% names(gtf)) {
      stop("attribute column not found in gtf data frame.")
    }
  }

  if (!is.null(features)) {
    gtf <- gtf[which(gtf$feature %in% features),]
  }
  if (nrow(gtf) == 0) {
    stop("No rows left in gtf, check features argument.")
  }

  if (process_attr_col) {
    #out<<- gtf
    #browser()
    #attribue_col <<- gtf$attribute
    #browser()
    message("processing the attribute column.")
    # use waldo::compare to compare results
    if (use_fun == "rcpp") {
      #Rcpp::sourceCpp(system.file("extdata/proc_gtf_attr.cpp", package = "igsc"))
      attr_col <- igsc:::process_attr_col_rcpp(gtf$attribute) #igsc:::
      attr_ind <- rep(seq_along(attr_col), lengths(attr_col)/2)
      attr_col <- unlist(attr_col)
    } else if (use_fun == "r") {
      attr_col <- stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)
      attr_ind <- rep(seq_along(attr_col), lengths(attr_col))
      attr_col <- unlist(stringi::stri_split_fixed(unlist(attr_col), pattern = " ", omit_empty = T)) #pattern = ' \"'
      attr_col <- stringi::stri_trim_both(stringi::stri_replace_all(attr_col, replacement = "", fixed = '"'))
    } else if (use_fun == "rust") {
      # rust fun was found slower; so using rcpp fun with proper registration
      rextendr::rust_source(system.file("extdata/lib.rs", package = "igsc"))
      attr_ind <- rep(seq_along(gtf$attribute), lengths(stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)))
      attr_col <- process_attr_col_rust(gtf$attribute) #igsc:::
    }

    if (attr_col_as_list == T) {
      gtf$attribute <- I(unname(split(stats::setNames(attr_col[seq(2, length(attr_col), 2)],
                                                      attr_col[seq(1, length(attr_col), 2)]),
                                      attr_ind)))
    } else {
      gtf <- gtf[,-which(names(gtf) == "attribute")]
      #gtf$attribute_names <- unname(purrr::map_chr(split(attr_col[seq(1, length(attr_col), 2)], attr_ind), paste, collapse = ","))
      #gtf$attribute_values <- unname(purrr::map_chr(split(attr_col[seq(2, length(attr_col), 2)], attr_ind), paste, collapse = ","))
    }

    attr_col <- data.frame(attribute = attr_col[seq(1, length(attr_col), 2)],
                           value = attr_col[seq(2, length(attr_col), 2)],
                           index = attr_ind)

    gtf$index <- seq(1, nrow(gtf), 1)
    gtf <- dplyr::left_join(gtf,
                            attr_col %>%
                              dplyr::filter(attribute %in% attr_keep) %>%
                              tidyr::pivot_wider(names_from = attribute, values_from = value), by = "index")

    return(list(gtf = gtf, attr = attr_col))
  } else {
    return(list(gtf = gtf, attr = NULL))
  }

}

get_bounds <- function(x, file_path) {
  # here we get the first and last line of a seqname to read with vroom
  # actually though, rg and grep do return the full lines allready, not only linenumbers
  # but, so what
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
  return(out)
}

vroom_gtf <- function(x, file_path, col_names, unpack_fun) {
  bounds <- get_bounds(x, file_path)
  y <- vroom::vroom(file = do.call(unpack_fun, args = list(description = file_path)),
                    col_names = col_names,
                    skip = bounds[1] - 1,
                    n_max = bounds[2] - bounds[1] + 1,
                    comment = "#",
                    show_col_types = F)
  return(y)
}


## too slow
#gtf_attr <- purrr::map(tt3, stringi::stri_split_fixed, pattern = " ", omit_empty = T, .progress = T)
#out <- purrr::map(gtf_attr, function(x) data.frame(attribute = sapply(x, "[", 2), value = sapply(x, "[", 1)), .progress = T)
#library(data.table)
#out2 <- rbindlist(lapply(out, as.data.frame.list), fill=TRUE)
