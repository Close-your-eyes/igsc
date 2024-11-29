#' From a genome fasta file, get chromosome sequence boundaries
#'
#' @param file_path path to fasta file
#'
#' @return
#' @export
#'
#' @examples
get_fasta_seq_bounds <- function(file_path) {

  # get width from second line, assuming this is the first line with sequence
  line_width <- nchar(trimws(vroom::vroom_lines(file = file_path,
                                                skip_empty_rows = T,
                                                skip = 1,
                                                n_max = 1)))
  total_lines <- tryCatch(
    {
      cmd <- paste0("rg --count '^' ", file_path)
      as.numeric(trimws(system(cmd, intern = T)))
    },
    error = function(err) {
      cmd <- paste0("grep --count '^' ", file_path)
      as.numeric(trimws(system(cmd, intern = T)))
    }
  )

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

  #name_lines <- sapply(strsplit(name_lines, " "), "[", 1)
  name_lines <- stats::setNames(sapply(strsplit(name_lines, "\\:>"), "[", 1), sapply(strsplit(name_lines, "\\:>"), "[", 2))
  name_lines <- utils::stack(name_lines)
  name_lines$values <- as.numeric(name_lines$values)
  name_lines$end_line <- c(name_lines$values[-1]-1, total_lines)
  names(name_lines)[1:2] <- c("start_line", "name")
  name_lines$name <- as.character(name_lines$name)
  name_lines$n_seq_lines <- name_lines$end_line - name_lines$start_line # exclude name line
  name_lines$bp_approx <- line_width*name_lines$n_seq_lines
  name_lines <- name_lines[,c(2,1,3,4,5)]

  return(name_lines)
}
# paste0(name_lines$end_line, "p")
# str(name_lines)
# cmd3 <- paste0("rg '^' -N ", file_path, " | sed -n '", paste(paste0(name_lines$end_line, "p"), collapse = "; "), "' | wc -m")
# cmd3 <- paste0("rg '^' -N ", file_path, " | awk 'NR==4149275 || NR==6379234 {print length($0)}'")
# tt <- system(cmd3, intern = T)
#
# cmd4 <- paste0("sed '4149275!d; 51662791!d' ", file_path)
# system(cmd4, intern = T)


