#' Convert cigar string to data frame with matched postions and elements from sequence
#'
#' @param cigar
#' @param start
#' @param seq
#' @param name
#' @param name_col
#' @param rm_clipped
#' @param skip_as
#'
#' @return
#' @export
#'
#' @examples
cigar_to_position <- function(cigar,
                              start,
                              seq,
                              name = NULL,
                              name_col = "seq.name",
                              rm_clipped = F,
                              skip_as = NA) {
  # https://davetang.org/wiki/tiki-index.php?page=SAM
  # https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/cigar.md

  if (missing(cigar)) {
    stop("cigar missing.")
  }
  if (missing(start)) {
    stop("start position for cigar conversion missing.")
  }
  if (missing(seq)) {
    stop("sequence (seq) for cigar conversion missing.")
  }

  # if val 1 is omited in cigar
  if (grepl("[[:alpha:]]{2}", cigar)) {
    message("Omitted 1 in cigar string suspected. Trying to insert the 1.")
    message(cigar)
    cigar <- add_1_between_consecutive_letters(cigar)
    message(cigar)
  }
  cigar_split <- strsplit(cigar, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl=TRUE)[[1]]
  val <- c(0,as.integer(cigar_split[grep("\\d", cigar_split, perl=TRUE)]))
  op <- cigar_split[grep("\\D", cigar_split, perl=TRUE)]

  if (any(!op %in% c("D", "F", "H", "I", "M", "N", "P", "R", "S", "X", "="))) {
    stop("Unknown operations in cigar string: ", paste(unique(op[which(!op %in% c("D", "F", "H", "I", "M", "N", "P", "R", "S", "X", "="))]), collapse = ","))
  }
  val_cum <- cumsum(val)

  # start position is affected by clipping at the beginning. change start respectively
  if (op[1] %in% c("S", "H")) {
    start <- start - val[2] # index two because 0 was added at index 1
  }

  # these are the positions in seq
  val_seq <- val[-1][which(op %in% c("S", "M", "=", "X"))]
  val_seq_cum <- c(0,cumsum(val_seq))

  seq_df <- data.frame(seq = character(val_cum[length(val_cum)]),
                       position = seq(start, start+val_cum[length(val_cum)]-1))
  # i is counter for op
  # j is counter for val_seq_cum
  j <- 1
  for (i in seq_along(op)) {
    if (op[i] %in% c("M", "X", "=")) {
      seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- strsplit(substr(seq, val_seq_cum[j]+1, val_seq_cum[j+1]), "")[[1]]
      j <- j + 1
    }

    # clipping means that respective bases were not uses for the alignment, but are retained in the output
    # only used at the end of reads
    # maybe because sequencing errors are more likely towards the ends?! or whatever ...
    if (op[i] == "S" && !rm_clipped) {
      seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- strsplit(substr(seq, val_seq_cum[j]+1, val_seq_cum[j+1]), "")[[1]]
      j <- j + 1
    } else if (op[i] == "S" && rm_clipped) {
      seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- rep(skip_as, val_seq_cum[j+1] - val_seq_cum[j])
      j <- j + 1
    }
    if (op[i] == "N") {
      seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- rep(skip_as, val[i+1])
      # j remains the same
    }
    if (op[i] %in% c("I", "D", "H", "F", "R")) {
      message("New operation found in cigar string. index: ", i)
      stop("New operation found in cigar string.")
    }
  }
  if (!is.null(name)) {
    seq_df[,name_col] <- name
  }
  return(seq_df)
}

add_1_between_consecutive_letters <- function(input_string) {
  result <- gsub("([A-Za-z])([A-Za-z])", "\\11\\2", input_string, perl = TRUE)
  return(result)
}
