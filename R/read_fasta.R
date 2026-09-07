#' Read sequences from a FASTA file
#'
#' Read a plain or gzip-compressed FASTA file with [vroom::vroom_lines()] and
#' return sequences named by their headers. Sequence lines belonging to the same
#' record are concatenated by default, so wrapped FASTA records are represented
#' as single strings.
#'
#' A subset of the file can be selected with inclusive `start_line` and
#' `end_line` bounds. Alternatively, `seqname` uses [get_fasta_names()] to locate
#' one record and adjusts those bounds. Files ending in `.tar.gz` are not
#' supported. If direct reading fails, the function attempts to decompress the
#' input into the session's temporary directory before reading it again.
#'
#' @param file A single path to a FASTA-formatted text file. Gzip-compressed
#'   files with a `.gz` extension are supported, but `.tar.gz` archives are not.
#' @param trimws Logical; trim leading and trailing whitespace from every input
#'   line before parsing. Disabling this can improve performance for trusted
#'   files.
#' @param rm_comments Logical; remove lines beginning with any value in
#'   `comment_indicator`. This can be slow for large files.
#' @param comment_indicator Character vector of prefixes identifying comment
#'   lines. Used only when `rm_comments = TRUE`.
#' @param rm_leading_arrow Logical retained for compatibility. It is currently
#'   ignored: the leading `">"` is always removed from returned record names.
#' @param make_names Logical; transform record names with `make_names_fun`.
#' @param make_names_fun Function, or its name, used when `make_names = TRUE`.
#'   It must accept the record names through an argument named `string`; the
#'   default is [janitor::make_clean_names()].
#' @param make_names_args Named list of additional arguments passed to
#'   `make_names_fun`.
#' @param start_line Numeric first file line to read, using one-based indexing.
#' @param end_line Numeric last file line to read, inclusive. Use `Inf` to read
#'   through the end of the file.
#' @param concat Logical; concatenate all sequence lines within each FASTA
#'   record. If `FALSE`, retain each record's sequence lines as a character
#'   vector in a list.
#' @param progress Logical; display the progress indicator from
#'   [vroom::vroom_lines()].
#' @param seqname Optional single record name to read. It may match a complete
#'   FASTA header, its first whitespace-delimited field, or the chromosome-style
#'   name prepared by [get_fasta_names()]. When supplied, it sets `start_line`
#'   to the selected header and, unless it is the final record, sets `end_line`
#'   to the line before the next header.
#'
#' @return If `concat = TRUE`, a named character vector with one concatenated
#'   sequence per FASTA record. If `concat = FALSE`, a named list of character
#'   vectors containing the original sequence lines. Names are FASTA headers
#'   without their leading `">"`, optionally transformed by `make_names_fun`.
#' @export
#'
#' @examples
#' \dontrun{
#' fasta <- tempfile(fileext = ".fa")
#' writeLines(
#'   c(">sequence 1", "ACGT", "TGCA", ">sequence_2", "NNNN"),
#'   fasta
#' )
#'
#' read_fasta(fasta, progress = FALSE)
#' read_fasta(fasta, concat = FALSE, progress = FALSE)
#' read_fasta(fasta, seqname = "sequence_2", progress = FALSE)
#' read_fasta(fasta, make_names = TRUE, progress = FALSE)
#'
#' unlink(fasta)
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
                       progress = T,
                       seqname = NULL) {

  igsc:::.ensure_packages(c("brathering", "janitor"))

  if (missing(file)) {
    stop("Please provide a path to a file in 'file'.")
  }

  if (!is.character(file)) {
    stop("file has to be the path to a fasta-formatted file.")
  }

  if (!file.exists(file)) {
    stop("file not found.")
  }
  if (grepl("\\.tar.gz$", file)) {
    stop("tar.gz files are not handled well. please untar or provide .gz file.")
  }

  if (!is.numeric(start_line) || !is.numeric(end_line) || length(start_line) == 0 || length(end_line) == 0) {
    stop("start_line or end_line is either not numeric or has length zero.")
  }

  make_names_fun <- match.fun(make_names_fun)

  if (!is.null(seqname)) {
    if (length(seqname) > 1) {
      stop("seqname can only be length 1.")
    }
    seqnames <- igsc::get_fasta_names(file)
    if (seqname %in% as.character(unlist(seqnames[,-1]))) {
      if (seqname %in% seqnames$fastaname) {
        ind <- which(seqnames$fastaname == seqname)
      } else if (seqname %in% seqnames$seqname) {
        ind <- which(seqnames$seqname == seqname)
      } else if (seqname %in% seqnames$fctname) {
        ind <- which(seqnames$fctname == seqname)
      }
      start_line <- seqnames[ind, "start_line"]
      if (ind < nrow(seqnames)) {
        end_line <- seqnames[ind+1, "start_line"]-1
      }
      message("set start_line: ", start_line, ", end_line: ", end_line)
    } else {
      stop("seqname not found.")
    }
  }


  # zip and gz files are handled well but .tar.gz not. there is a problem with connection size then
  lines <- tryCatch(
    expr = {
      vroom::vroom_lines(file = file,
                         skip = start_line - 1,
                         n_max = end_line - start_line + 1,
                         skip_empty_rows = T,
                         progress = progress)
    },
    error = function(err) {
      print(err)
      out <- brathering::ungunzip(
        file,
        out_dir = tempdir(),
        out_file = tools::file_path_sans_ext(basename(file))
      )
      message("unpacking file to: ", out)
      lines <- vroom::vroom_lines(file = out,
                                  skip = start_line - 1,
                                  n_max = end_line - start_line + 1,
                                  skip_empty_rows = T,
                                  progress = progress)
      return(lines)
    }
  )

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


  #ind <- which(stringi::stri_startswith_fixed(pattern = ">", from = 1, str = lines))
  ind <- which(startsWith(lines, ">"))

  if (!length(ind)) {
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
    lines <- purrr::map_chr(brathering::seq2(start, end), ~stringi::stri_paste(lines[.x], collapse = ""))
  } else {
    lines <- split(lines, rep(seqnames, end-start+2))
  }
  names(lines) <- gsub("^>", "", seqnames)
  if (make_names) {
    names(lines) <- do.call(make_names_fun, args = c(list(string = names(lines)), make_names_args))
  }
  return(lines)
}

#test <- read_fasta("/Users/chris/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genome.fa.gz")#, start_line = 1, end_line = 4149276+20)
