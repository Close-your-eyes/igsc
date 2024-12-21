#' Title
#'
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
get_fasta_names <- function(file_path) {

  name_lines <- tryCatch(
    {
      cmd <- paste0("rg -n '^>' ", file_path)
      system(cmd, intern = T)
    },
    error = function(err) {
      cmd <- paste0("grep -n '^>' ", file_path)
      system(cmd, intern = T)
    }
  )

  split <- strsplit(name_lines, "\\:>")
  name_lines <- stats::setNames(sapply(split, "[", 1), sapply(split, "[", 2))
  name_lines <- utils::stack(name_lines)
  names(name_lines)[1:2] <- c("start_line", "fastaname")
  name_lines$start_line <- as.numeric(name_lines$start_line)
  name_lines$fastaname <- as.character(name_lines$fastaname)
  name_lines$seqname <- sapply(strsplit(name_lines$fastaname, " "), "[", 1)
  return(name_lines)
}
