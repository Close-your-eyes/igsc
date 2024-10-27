#' From a genome fasta file, get chromosome sequence boundaries
#'
#' @param file_path path to fasta file
#'
#' @return
#' @export
#'
#' @examples
get_fasta_seq_bounds <- function(file_path) {
  file_path <- "/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genome.fa"
  #cmd1 <- paste0("wc -l ", file_path) # rg --count '^' filename.txt
  cmd1 <- paste0("rg --count '^' ", file_path) # rg --count '^' filename.txt
  total_lines <- as.numeric(trimws(system(cmd1, intern = T)))
  #total_lines <- strsplit(total_lines, " ")[[1]][1]
  #cmd2 <- paste0("grep -n '^>' ", file_path) # paste0("rg -n '^>' ", file_path)
  cmd2 <- paste0("rg -n '^>' ", file_path) # paste0("rg -n '^>' ", file_path)
  out <- system(cmd2, intern = T)
  out <- sapply(strsplit(out, " "), "[", 1)
  out <- stats::setNames(sapply(strsplit(out, "\\:>"), "[", 1), sapply(strsplit(out, "\\:>"), "[", 2))
  out <- utils::stack(out)
  out$values <- as.numeric(out$values)
  out$end_line <- c(out$values[-1]-1, total_lines)
  names(out)[1:2] <- c("start_line", "name")
  out$name <- as.character(out$name)
  out$n_seq_lines <- out$end_line - out$start_line # exclude name line
  out$bp_approx <- 60*out$n_seq_lines
  out <- out[,c(2,1,3,4,5)]

  return(out)
}
# paste0(out$end_line, "p")
# str(out)
# cmd3 <- paste0("rg '^' -N ", file_path, " | sed -n '", paste(paste0(out$end_line, "p"), collapse = "; "), "' | wc -m")
# cmd3 <- paste0("rg '^' -N ", file_path, " | awk 'NR==4149275 || NR==6379234 {print length($0)}'")
# tt <- system(cmd3, intern = T)
#
# cmd4 <- paste0("sed '4149275!d; 51662791!d' ", file_path)
# system(cmd4, intern = T)


