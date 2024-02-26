#' Read and process a GTF (Gene Transfer Format) file
#'
#' See https://www.ensembl.org/info/website/upload/gff.html?redirect=no
#'
#' @param file_path path to the file
#' @param ... arguments to vroom::vroom
#'
#' @return a list with (i) entries of the GTF file including the attribute column as list and some attributes as separate columns and (ii) the attributes as long data frame
#' @export
#'
#' @examples
read_gtf <- function(file_path,
                     ...) {

  gtf <- vroom::vroom(file_path, col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), comment = "#", ...)
  message("processing the attribute column.")
  attr_col <- stringi::stri_split_fixed(gtf$attribute, pattern = ";", omit_empty = T)
  attr_ind <- rep(seq_along(attr_col), lengths(attr_col))
  attr_col <- unlist(stringi::stri_split_fixed(unlist(attr_col), pattern = ' \"', omit_empty = T))
  attr_col <- stringi::stri_trim_both(stringi::stri_replace_all(attr_col, replacement = "", fixed = '"'))
  attr_col <- stats::setNames(attr_col[seq(2, length(attr_col), 2)], attr_col[seq(1, length(attr_col), 2)])
  gtf$attribute <- I(unname(split(attr_col, attr_ind)))
  gtf$index <- seq(1, nrow(gtf), 1)

  attr_col <- data.frame(attribute = names(attr_col), value = attr_col, index = attr_ind)
  gtf <- dplyr::left_join(gtf, attr_col %>% dplyr::filter(attribute %in% c("gene_id", "gene_name", "transcript_id", "exon_number")) %>% tidyr::pivot_wider(names_from = attribute, values_from = value), by = "index")

  return(list(gtf = gtf, attr = attr_col))
}


## too slow
#gtf_attr <- purrr::map(tt3, stringi::stri_split_fixed, pattern = " ", omit_empty = T, .progress = T)
#out <- purrr::map(gtf_attr, function(x) data.frame(attribute = sapply(x, "[", 2), value = sapply(x, "[", 1)), .progress = T)
#library(data.table)
#out2 <- rbindlist(lapply(out, as.data.frame.list), fill=TRUE)
