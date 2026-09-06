#' Expand a CIGAR string into aligned sequence positions
#'
#' Convert a CIGAR string and its query sequence into a position-by-position data
#' frame. Matched (`M`), sequence-match (`=`), and sequence-mismatch (`X`)
#' operations copy bases from `seq`; deletions (`D`) are represented by `"-"`;
#' skipped regions (`N`) are represented by `skip_as`; and insertions (`I`) are
#' omitted from the returned reference-coordinate rows.
#'
#' Soft-clipped (`S`) bases are retained by default. An initial soft clip shifts
#' the first returned position upstream by its length so that clipped bases are
#' included before `start`. Set `rm_clipped = TRUE` to replace soft-clipped bases
#' with `skip_as` instead.
#'
#' If two CIGAR operations occur consecutively without an explicit length, the
#' function assumes the omitted length is one and inserts it before parsing.
#'
#' @param cigar A single CIGAR string. Operation lengths must be integers;
#'   omitted lengths of one between consecutive operations are inferred.
#' @param start Integer genomic start position of the alignment. When `cigar`
#'   begins with soft clipping, the returned positions begin before this value
#'   so that the clipped bases can be represented.
#' @param seq A single character string containing the query sequence described
#'   by `cigar`.
#' @param name Optional sequence or read identifier. When non-`NULL`, it is
#'   repeated in an additional output column named by `name_col`.
#' @param name_col Name of the output column containing `name`.
#' @param rm_clipped Logical; replace soft-clipped bases with `skip_as` rather
#'   than retaining their values from `seq`.
#' @param skip_as Scalar value used for skipped (`N`) positions and, when
#'   `rm_clipped = TRUE`, soft-clipped (`S`) positions. The default is `NA`.
#'
#' @return A data frame with columns `seq` and `position`, containing one row per
#'   returned reference-coordinate position. If `name` is supplied, the result
#'   also contains a column named by `name_col`.
#' @export
#'
#' @examples
#' # Three matches, two skipped reference positions, then two matches.
#' cigar_to_position(
#'   cigar = "3M2N2M",
#'   start = 100,
#'   seq = "ACGTT",
#'   name = "read_1"
#' )
#'
#' # Soft-clipped bases are retained and positioned before the alignment start.
#' cigar_to_position(cigar = "2S3M", start = 100, seq = "AACGT")
#'
#' # They can instead be replaced with the value supplied to skip_as.
#' cigar_to_position(
#'   cigar = "2S3M",
#'   start = 100,
#'   seq = "AACGT",
#'   rm_clipped = TRUE,
#'   skip_as = "N"
#' )
cigar_to_position <- function(cigar,
                              start,
                              seq,
                              name = NULL,
                              name_col = "seq.name",
                              rm_clipped = F,
                              skip_as = NA) {
                              #ref_start_pos = NULL) {
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

  # if val 1 is omitted in cigar
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
  val_seq <- val[-1][which(op %in% c("S", "M", "=", "X", "I", "D"))]
  val_seq_cum <- c(0,cumsum(val_seq))

  seq_df <- data.frame(seq = character(val_cum[length(val_cum)]),
                       position = seq(start, start+val_cum[length(val_cum)]-1))
  # i is counter for op
  # j is counter for val_seq_cum
  j <- 1
  for (i in seq_along(op)) {
    if (op[i] %in% c("M", "X", "=")) {
      #browser()
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
    if (op[i] == "D") {
      #browser()
      seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- rep("-", val[i+1])
      #val_cum[(i+1):length(val_cum)] <- val_cum[(i+1):length(val_cum)] - val[i+1]
      val_seq_cum[(j+1):length(val_seq_cum)] <- val_seq_cum[(j+1):length(val_seq_cum)] - val_seq[j]
      j <- j + 1
    }
    if (op[i] == "I") {
      # to allow for insertion would require to pass the ref seq and totally rewrite the function
      message("insertion skipped.")
      #seq_df$seq[(val_cum[i]+1):val_cum[i+1]] <- substr(seq, val_seq_cum[j]+1, val_seq_cum[j+1])
      val_cum[(i+1):length(val_cum)] <- val_cum[(i+1):length(val_cum)] - val[i+1]
      #val_seq_cum[(j+1):length(val_seq_cum)] <- val_seq_cum[(j+1):length(val_seq_cum)] - val_seq[j]
      seq_df <- seq_df[-c((nrow(seq_df)-val[i+1]+1):nrow(seq_df)),]
      j <- j + 1
    }
    if (op[i] %in% c("H", "F", "R")) { # D
      message("New operation found in cigar string. index: ", i)
      stop("New operation found in cigar string.")
    }
  }
  if (!is.null(name)) {
    seq_df[,name_col] <- name
  }
  # seq_df$position_1 <- seq_df$position - min(seq_df$position) + 1
  # if (!is.null(ref_start_pos)) {
  #   seq_df$position_rel_start <- seq_df$position - ref_start_pos + 1
  # }
  return(seq_df)
}

add_1_between_consecutive_letters <- function(input_string) {
  result <- gsub("([A-Za-z])([A-Za-z])", "\\11\\2", input_string, perl = TRUE)
  return(result)
}
