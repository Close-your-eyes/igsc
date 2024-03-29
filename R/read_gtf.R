#' Read and process a GTF (Gene Transfer Format) file
#'
#' See https://www.ensembl.org/info/website/upload/gff.html?redirect=no
#'
#' @param file_path path to the file
#' @param gtf data frame from reading gtf file with vroom or similar.
#' useful to pass a subset of rows only.
#' @param attr_col_as_list have the attributes column return as named list (TRUE) or
#' separated by names and values into two columns
#' @param ... arguments to vroom::vroom
#'
#' @return a list with (i) entries of the GTF file including the attribute column as list and some attributes as separate columns and (ii) the attributes as long data frame
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#'
#' gtf <- read_gtf(your_path, attr_col_as_list = F)
#' # when attr_col_as_list = F attributes are split into names and values columns
#' # this is how to expand the attributes names and vals columns
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
                     ...) {

  # file_path <- "/Volumes/CMS_SSD_2TB/2023_UriSeq/RNA_bulk_sequencing/GRCh38.p14/Homo_sapiens.GRCh38.111.gtf"

  if (missing(file_path) && missing(gtf)) {
    stop("Provide file_path or gtf data frame.")
  }

  if (!missing(file_path) && !missing(gtf)) {
    message("file_path and gtf provided. Will ignore file_path and work with the gtf data frame provided.")
  }

  col_names <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

  if (missing(gtf)) {
    gtf <- vroom::vroom(file_path,
                        col_names = col_names,
                        comment = "#")
  } else {

    if (!"attribute" %in% names(gtf)) {
      stop("attribute column not found in gtf data frame.")
    }

  }


  message("processing the attribute column.")
  attr_col <- stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)
  attr_ind <- rep(seq_along(attr_col), lengths(attr_col))
  attr_col <- unlist(stringi::stri_split_fixed(unlist(attr_col), pattern = ' \"', omit_empty = T))
  attr_col <- stringi::stri_trim_both(stringi::stri_replace_all(attr_col, replacement = "", fixed = '"'))

  if (attr_col_as_list == T) {
    gtf$attribute <- I(unname(split(stats::setNames(attr_col[seq(2, length(attr_col), 2)],
                                                    attr_col[seq(1, length(attr_col), 2)]),
                                    attr_ind)))
  } else {
    gtf <- gtf[,-which(names(gtf) == "attribute")]
    gtf$attribute_names <- unname(purrr::map_chr(split(attr_col[seq(1, length(attr_col), 2)], attr_ind), paste, collapse = ","))
    gtf$attribute_values <- unname(purrr::map_chr(split(attr_col[seq(2, length(attr_col), 2)], attr_ind), paste, collapse = ","))
  }

  attr_col <- data.frame(attribute = attr_col[seq(1, length(attr_col), 2)],
                         value = attr_col[seq(2, length(attr_col), 2)],
                         index = attr_ind)

  gtf$index <- seq(1, nrow(gtf), 1)
  gtf <- dplyr::left_join(gtf,
                          attr_col %>%
                            dplyr::filter(attribute %in% c("gene_id", "gene_name", "transcript_id", "exon_number")) %>%
                            tidyr::pivot_wider(names_from = attribute, values_from = value), by = "index")

  return(list(gtf = gtf, attr = attr_col))
}




## too slow
#gtf_attr <- purrr::map(tt3, stringi::stri_split_fixed, pattern = " ", omit_empty = T, .progress = T)
#out <- purrr::map(gtf_attr, function(x) data.frame(attribute = sapply(x, "[", 2), value = sapply(x, "[", 1)), .progress = T)
#library(data.table)
#out2 <- rbindlist(lapply(out, as.data.frame.list), fill=TRUE)
