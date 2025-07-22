#' Title
#'
#' intended for
#'
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
get_fasta_names <- function(file_path,
                            prep_chr = T) {

  name_lines <- tryCatch(
    {
      if (tools::file_ext(file_path) == "gz") {
        cmd <- paste0("gunzip -c ", file_path, " | rg -n '^>' ")
      } else {
        cmd <- paste0("rg -n '^>' ", file_path)
      }
      system(cmd, intern = T)
    },
    error = function(err) {
      if (tools::file_ext(file_path) == "gz") {
        cmd <- paste0("gunzip -c ", file_path, " | grep -n '^>' ")
      } else {
        cmd <- paste0("grep -n '^>' ", file_path)
      }
      system(cmd, intern = T)
    }
  )

  ## factor names
  split <- strsplit(name_lines, "\\:>")
  name_lines <- stats::setNames(sapply(split, "[", 1), sapply(split, "[", 2))
  name_lines <- utils::stack(name_lines)
  names(name_lines)[1:2] <- c("start_line", "fastaname")
  name_lines$start_line <- as.numeric(name_lines$start_line)
  name_lines$fastaname <- as.character(name_lines$fastaname)

  if (prep_chr) {
    name_lines$seqname <- sapply(strsplit(name_lines$fastaname, " "), "[", 1)

    name <- gsub("chr", "", name_lines$seqname)
    name_num <- suppressWarnings(as.numeric(name))
    name_ord1 <- sort(name_num[which(!is.na(name_num))])
    name_ord2 <- name[which(is.na(name_num))]
    name_lines$fctname <- factor(name, c(name_ord1, name_ord2))

    name_lines$seqname <- factor(name_lines$seqname, levels = name_lines$seqname[order(name_lines$fctname)])
    name_lines$fastaname <- factor(name_lines$fastaname, levels = name_lines$fastaname[order(name_lines$fctname)])
  }

  return(name_lines)
}
