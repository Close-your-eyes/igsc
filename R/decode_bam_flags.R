#' Decode SAM/BAM alignment flags
#'
#' Convert SAM/BAM flag values into integer indicator columns describing
#' alignment properties. Multiple columns can be active for each alignment.
#'
#' @param flag An integer vector of SAM/BAM flag values, or a numeric vector
#'   coercible to integers. Valid flags are integers from 0 through 65535.
#'   Missing values produce missing indicators. Input values are not validated.
#' @param drop_constant Logical; if `TRUE`, remove columns with only one
#'   distinct value across reads. Missing values count as a distinct value.
#'   Empty results retain all columns.
#'
#' @return A data frame with one row per input value. Each column contains
#'   `1L` when the property applies, `0L` otherwise, or `NA_integer_` for missing
#'   input values:
#'   \describe{
#'     \item{paired}{The read belongs to a paired-end or other multi-segment template.}
#'     \item{proper_pair}{Segments are properly aligned according to the aligner.}
#'     \item{unmapped}{The read is unmapped.}
#'     \item{mate_unmapped}{The mate or next segment is unmapped.}
#'     \item{reverse_complemented}{The stored read sequence is reverse-complemented.}
#'     \item{mate_reverse_complemented}{The stored mate or next-segment sequence
#'       is reverse-complemented.}
#'     \item{first_segment}{The first segment in the template, usually R1.}
#'     \item{last_segment}{The last segment in the template, usually R2.}
#'     \item{secondary}{A secondary alignment, typically an alternative mapping.}
#'     \item{qc_failed}{The read failed platform or vendor quality checks.}
#'     \item{duplicate}{The read is marked as a PCR or optical duplicate.}
#'     \item{supplementary}{A supplementary alignment, representing another part
#'       of a split or chimeric alignment.}
#'     \item{mapped}{The unmapped flag is absent.}
#'     \item{primary}{Neither the secondary nor supplementary flag is set.}
#'     \item{plus_strand}{The read is mapped and is not reverse-complemented.}
#'     \item{minus_strand}{The read is mapped and is reverse-complemented.}
#'   }
#'
#' @details
#' The first twelve columns decode individual SAM/BAM flag bits. The remaining
#' four columns describe properties derived from those bits.
#'
#' Interpret `proper_pair`, mate-related indicators, `first_segment`, and
#' `last_segment` only when `paired` is active. For paired-end reads,
#' `first_segment` and `last_segment` usually identify R1 and R2; they do not
#' indicate genomic strand.
#'
#' Unmapped reads have both strand indicators set to zero. An unmapped read
#' can still have `primary` set to one. Primary status does not imply unique
#' mapping, high mapping quality, or absence of duplicate marking.
#'
#' Reserved flag bits are ignored. An empty input vector returns a data frame
#' with zero rows and all sixteen columns.
#'
#' @references
#' \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#'
#' @examples
#' # Properly paired R1 on the plus strand, with its mate on the minus strand
#' decode_bam_flags(99L)
#'
#' # Decode multiple alignments
#' decode_bam_flags(c(0L, 16L, 99L, 147L, 256L, 2048L, NA_integer_))
#'
#' # Append indicators to a data frame of alignments
#' reads <- data.frame(qname = c("read1", "read2"), flag = c(99L, 147L))
#' reads <- cbind(reads, decode_bam_flags(reads$flag))
#'
#' @export
decode_bam_flags <- function(flag,
                             drop_constant = T) {

  bits <- c(
    paired = 1L,
    proper_pair = 2L,
    unmapped = 4L,
    mate_unmapped = 8L,
    reverse_complemented = 16L,
    mate_reverse_complemented = 32L,
    first_segment = 64L,
    last_segment = 128L,
    secondary = 256L,
    qc_failed = 512L,
    duplicate = 1024L,
    supplementary = 2048L
  )

  if (!is.numeric(flag) || is.complex(flag) || is.object(flag)) {
    stop(
      "`flag` must be a plain integer or numeric vector.",
      call. = FALSE
    )
  }

  present <- !is.na(flag)
  values <- flag[present]

  if (any(!is.finite(values) |
          values < 0 |
          values > 65535 |
          values != trunc(values))) {
    stop(
      "`flag` must contain whole numbers from 0 through 65535, or NA.",
      call. = FALSE
    )
  }

  flag <- as.integer(flag)

  decoded <- as.data.frame(lapply(bits, function(bit) {
    as.integer(bitwAnd(as.integer(flag), bit) != 0L)
  }))

  decoded$mapped <- as.integer(decoded$unmapped == 0L)

  decoded$primary <- as.integer(
    decoded$secondary == 0L &
      decoded$supplementary == 0L
  )

  decoded$plus_strand <- as.integer(
    decoded$mapped == 1L &
      decoded$reverse_complemented == 0L
  )

  decoded$minus_strand <- as.integer(
    decoded$mapped == 1L &
      decoded$reverse_complemented == 1L
  )

  if (drop_constant && nrow(decoded) > 0L) {
    keep <- vapply(decoded, function(x) length(unique(x)) > 1L, logical(1))

    if (any(!keep)) {
      message(
        "Removed columns with one distinct value across all reads: ",
        paste(names(decoded)[!keep], collapse = ", "),
        "."
      )
    }

    decoded <- decoded[, keep, drop = FALSE]
  }

  return(decoded)
}
