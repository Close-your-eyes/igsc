#' Read sequences from a fasta-formatted file
#'
#' vroom::vroom_lines is used to quickly read lines of a text file (as a fasta file is one).
#'
#'
#' @param file path to a file with fasta-formatted sequences; file may
#' be gunzipped (ending with .gz)
#' @param rm_comments remove comment lines? this is very slow large files, currently.
#' @param comment_indicator first characters of a line which indicate a comment
#' @param trimws trim leading and trailing whitespaces? disable in case of trustworthy files. may take time for large fasta files
#' @param rm_leading_arrow remove leading arrow in fasta names?
#' @param make_names apply make_names to fasta names?
#' @param start_line passed to vroom::vroom_lines
#' @param end_line passed to vroom::vroom_lines to determine n_max
#' @param make_names_fun
#' @param make_names_args
#' @param concat
#' @param mc.cores
#'
#' @return a named character vector of sequences
#' @export
#'
#' @examples
#' \dontrun{
#' seqs <- igsc::read_fasta(file = "my/fasta/file.fasta")
#' seqs_zipped <- igsc::read_fasta(file = "my/fasta/file.fa.gz")
#' }
read_fasta <- function(file,
                       trimws = F,
                       rm_comments = F,
                       comment_indicator = c(";", "#"),
                       rm_leading_arrow = T,
                       make_names = F, # option for make_names fun e.g. from janitor
                       make_names_fun = janitor::make_clean_names,
                       make_names_args = list(),
                       start_line = 1,
                       end_line = Inf,
                       concat = T,
                       mc.cores = 1) {

  if (missing(file)) {
    stop("Please provide a path to a file in 'file'.")
  }

  if (!is.character(file)) {
    stop("file has to be the path to a fasta-formatted file.")
  }

  if (!file.exists(file)) {
    stop("file not found.")
  }

  if (!is.numeric(start_line) || !is.numeric(end_line) || length(start_line) == 0 || length(end_line) == 0) {
    stop("start_line or end_line is either not numeric or has length zero.")
  }

  make_names_fun <- match.fun(make_names_fun)

  if (grepl("\\.gz$", file)) {
    unpack_fun <- gzfile
  } else {
    unpack_fun <- function(description) {description}
  }
  lines <- vroom::vroom_lines(file = do.call(unpack_fun, args = list(description = file)),
                              skip = start_line - 1,
                              n_max = end_line - start_line + 1,
                              skip_empty_rows = T,
                              progress = F)
  if (trimws) {
    lines <- stringi::stri_trim_both(lines)
  }

  # this is slow !!
  if (rm_comments) {
    message("looking for comment lines.")
    rm_inds <- unique(unlist(purrr::map(comment_indicator, function(x) {
      which(stringi::stri_startswith_fixed(pattern = x, from = 1, str = lines))
    })))
    if (length(rm_inds) > 0) {
      message(length(rm_inds), " comment lines removed")
      lines <- lines[!which(seq_along(lines) %in% rm_inds)]
    }
  }

  ind <- which(stringi::stri_startswith_fixed(pattern = ">", from = 1, str = lines))
  if (length(ind) == 0) {
    message("fasta-formated sequences (names starting with '>') not found. Using the file name.")
    lines <- c(lines, paste0(">", gsub(paste0("\\.", tools::file_ext(file), "$"), "", basename(file))))
    ind <- c(1, ind+1)
  }
  # if length(ind) was 0 above, now there is at least ind = 1
  if (!1 %in% ind) {
    # this elegantly removes leading comment lines
    message("first row should be a sequence name. removing all lines before first fasta name.")
    lines <- lines[-c(1:(ind[1]-1))]
  }

  # this is done to compensate linebreaks in sequences
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  seqnames <- lines[ind]

  if (concat) {
    if (mc.cores == 1) {
      lines <- purrr::map_chr(igsc:::seq2(start, end), function(x) stringi::stri_paste(lines[x], collapse = ""))
    } else {
      chunk_sizes <- ceiling(lengths(igsc:::seq2(start, end)) / mc.cores)
      lines <- parallel::mcmapply(x = igsc:::seq2(start, end),
                                  y = chunk_sizes, function(x,y) split(lines[x], ceiling(seq_along(lines[x]) / y)),
                                  SIMPLIFY = F,
                                  mc.cores = mc.cores)
      lines <- lapply(lines, function(x) paste(unlist(parallel::mclapply(X = x,
                                                                         FUN = stringi::stri_paste,
                                                                         collapse = "",
                                                                         mc.cores = mc.cores)),
                                               collapse = ""))
    }
  } else {
    lines <- split(lines, rep(seqnames, end-start+2))
  }
  names(lines) <- gsub("^>", "", seqnames)
  if (make_names) {
    names(lines) <- do.call(make_names_fun, args = c(list(string = names(lines)), make_names_args))
  }
  return(lines)
}

#test <- read_fasta("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genome.fa", start_line = 1, end_line = 4149276+20)
