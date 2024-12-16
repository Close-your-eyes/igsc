#' From a genome fasta file, get chromosome sequence boundaries
#'
#' @param file_path path to fasta file or fst file
#' @param save_fst
#' @param check_fst
#'
#' @return
#' @export
#'
#' @examples
get_fasta_seq_bounds <- function(file_path,
                                 save_fst = T,
                                 check_fst = T) {

  if (tools::file_ext(basename(file_path)) == "fst") {
    return(fst::read_fst(file_path))
  }
  if (check_fst) {
    fst_file <- file.path(dirname(file_path), paste0(basename(file_path), "_seq_bounds.fst"))
    if (file.exists(fst_file)) {
      message("reading fst file.")
      return(fst::read_fst(fst_file))
    }
  }

  # get width from second line, assuming this is the first line with sequence
  line_width <- nchar(trimws(vroom::vroom_lines(file = file_path,
                                                skip_empty_rows = T,
                                                skip = 1,
                                                n_max = 1,
                                                progress = F)))

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

  name_lines <- get_fasta_names(file_path)
  name_lines$end_line <- c(name_lines$start_line[-1]-1, total_lines)
  name_lines$n_seq_lines <- name_lines$end_line - name_lines$start_line # exclude name line
  name_lines$bp_approx <- line_width*name_lines$n_seq_lines
  name_lines <- name_lines[,c(2,3,1,4:6)]

  if (save_fst) {
    fst::write_fst(name_lines, path = file.path(dirname(file_path), paste0(basename(file_path), "_seq_bounds.fst")))
  }
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


