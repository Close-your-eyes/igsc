# R wrappers for the Bash functions in igsc/inst/extdata/compress_fastq/.
# Requires Bash, pigz, awk, and Python 3 on PATH as used by the Bash script.
# For a development checkout, set options(igsc.fastq_script = "/path/to/script").
# Otherwise the script is located in the installed igsc package.

#' Run a FASTQ Bash function
#'
#' Source the FASTQ utility script in a separate Bash process and invoke one of
#' its functions. Paths and arguments are quoted before passing them to Bash.
#'
#' @param fun Character scalar naming the Bash function to invoke.
#' @param input Character scalar giving the input file path.
#' @param output Character scalar giving the output file path.
#' @param args Character vector of additional command-line options and values.
#' @param script Optional path to the Bash utility script. When `NULL`, use
#'   `getOption("igsc.fastq_script")`, then the script in the installed `igsc`
#'   package if the option is unset.
#'
#' @returns `output`, invisibly. An error is raised if the script or Bash cannot
#'   be found, the paths are invalid, or the Bash function returns a nonzero status.
#' @keywords internal
#' @noRd
.run_fastq_bash <- function(fun, input, output, args = character(), script = NULL) {
  if (is.null(script)) script <- getOption("igsc.fastq_script")
  if (is.null(script)) {
    script <- system.file(
      "extdata", "compress_fastq", "fastq_compression.bash",
      package = "igsc"
    )
  }
  if (length(script) != 1L || is.na(script) || !nzchar(script) ||
      !file.exists(script)) {
    stop("Set options(igsc.fastq_script = '/path/to/fastq_compression.bash') ",
         "or supply script= with an existing Bash script.", call. = FALSE)
  }
  for (path in list(input, output)) {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
      stop("input and output must each be one nonempty path.", call. = FALSE)
    }
  }
  bash <- Sys.which("bash")
  if (!nzchar(bash)) stop("Bash is required on PATH.", call. = FALSE)
  script <- normalizePath(script, mustWork = TRUE)
  argv <- c(script, fun, input, output, args)
  status <- system2(
    bash,
    args = c(
      "-c", shQuote('source "$1" && shift && "$@"'),
      "fastq-r", shQuote(argv)
    )
  )
  if (status != 0L) {
    stop(fun, " failed with exit status ", status, ".", call. = FALSE)
  }
  invisible(output)
}

#' Format an optional FASTQ command-line argument
#'
#' Construct an option and value pair without scientific notation, or omit the
#' option when its value is `NULL`.
#'
#' @param option Character scalar giving a command-line option, such as
#'   `"--reads"`.
#' @param value Numeric value to format, or `NULL` to omit the option.
#'
#' @returns A character vector containing the option and its formatted value,
#'   or an empty character vector when `value` is `NULL`.
#' @keywords internal
#' @noRd
.fastq_option <- function(option, value) {
  if (is.null(value)) return(character())
  c(option, format(value, scientific = FALSE, trim = TRUE))
}


#' Recompress FASTQ with optional read names and constant qualities
#'
#' Stream a gzip-compressed FASTQ file through the Bash utility, optionally
#' replacing read headers and quality strings before recompressing with pigz.
#'
#' @param input Character scalar giving the path to a readable gzip-compressed
#'   FASTQ file with four lines per read.
#' @param output Character scalar giving the destination gzip file. It must
#'   differ from `input`. Existing output is overwritten, regardless of suffix.
#' @param rename Logical; replace read headers with `@read1`, `@read2`, and so on
#'   in input order. Numbering restarts at one per call. Default: `FALSE`.
#' @param quality Logical; replace every quality line with repeated `I`
#'   characters matching the first read length. Default: `FALSE`.
#' @param level Integer compression level from 0 to 9. Default: 6.
#' @param reads Positive integer limiting output to the first requested reads,
#'   or `NULL` to retain all reads. Fewer available reads are retained in full.
#' @param threads Positive integer passed to both pigz processes, or `NULL` for
#'   pigz's default thread count. The awk transformation remains serial.
#' @param script Optional character path to `fastq_compression.bash`. When
#'   `NULL`, use `getOption("igsc.fastq_script")` if set; otherwise use
#'   `system.file("extdata", "compress_fastq", "fastq_compression.bash",
#'   package = "igsc")` from the installed package.
#'
#' @details
#' Bash, pigz, and awk must be available on the R session's `PATH`.
#' Input is assumed to be unwrapped ASCII FASTQ. Sequence and separator lines
#' are preserved. Equal read lengths are assumed when replacing qualities;
#' read lengths and the full FASTQ structure are not validated by this function.
#' With both replacement options disabled, the input is only recompressed.
#' Limiting reads can stop decompression before the entire input is checked.
#' Errors from the Bash utility are written to stderr, and failures may leave
#' partial output files.
#'
#' @returns The `output` path, invisibly. Raises an R error if the Bash utility
#'   cannot be run or reports failure.
#' @seealso [compress_fastq_sequences()], [decode_fastq_triplets()]
#' @export
#'
#' @examples
#' \dontrun{
#' relabel_fastq(
#'   "input.fastq.gz", "renamed.fastq.gz",
#'   rename = TRUE, quality = TRUE, level = 6L, threads = 8L
#' )
#' relabel_fastq("input.fastq.gz", "first10000.fastq.gz", reads = 10000L)
#' }
relabel_fastq <- function(input, output, rename = FALSE, quality = FALSE,
                          level = 6L, reads = NULL, threads = NULL,
                          script = NULL) {
  args <- c(
    if (isTRUE(rename)) "--rename",
    if (isTRUE(quality)) "--quality",
    .fastq_option("--level", level),
    .fastq_option("--reads", reads),
    .fastq_option("--threads", threads)
  )
  .run_fastq_bash("relabel_fastq", input, output, args, script)
}


#' Compress concatenated FASTQ sequences with optional triplet encoding
#'
#' Discard FASTQ headers, separator lines, and qualities, then concatenate the
#' sequences in input order without line breaks and gzip-compress the result.
#' Optionally apply the reversible custom triplet encoding before compression.
#'
#' @param input Character scalar giving the path to a readable gzip-compressed
#'   FASTQ file with four lines per read.
#' @param output Character scalar giving the destination gzip file. It must
#'   differ from `input`. Existing output is overwritten, regardless of suffix.
#' @param read_length Positive integer specifying the exact length of every
#'   selected read. Default: 98. Unexpected lengths cause an error.
#' @param triplets Logical; enable the custom triplet encoding before gzip
#'   compression. Default: `FALSE`, which stores the concatenated bases directly.
#' @param level Integer compression level from 0 to 9. Default: 6.
#' @param reads Positive integer limiting output to the first requested reads,
#'   or `NULL` to retain all reads. Fewer available reads are retained in full.
#' @param threads Positive integer passed to both pigz processes, or `NULL` for
#'   pigz's default thread count. The awk transformation remains serial.
#' @param script Optional character path to `fastq_compression.bash`. When
#'   `NULL`, use `getOption("igsc.fastq_script")` if set; otherwise locate the
#'   script under `extdata/compress_fastq` in the installed `igsc` package.
#'
#' @details
#' Bash, pigz, and awk must be available on the R session's `PATH`. Input is
#' assumed to be unwrapped ASCII FASTQ. Selected read lengths and incomplete
#' selected records are checked, but this is not a full FASTQ validator.
#' A trailing carriage return is removed from each sequence line.
#'
#' Triplet encoding maps non-overlapping A/C/G/T triplets in lexical order to
#' `ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/`.
#' For example, `AAA` maps to `A` and `TTT` maps to `/`. Literal escapes preserve
#' triplets with ambiguous bases or lowercase characters, and a final escape
#' preserves one or two trailing bases. Triplet grouping continues across read
#' boundaries. A smaller gzip output is not guaranteed by this encoding.
#'
#' The output does not record its encoding mode or read length. Keep both
#' settings and pass them to [decode_fastq_triplets()] during restoration.
#' Original read names and quality strings are discarded. Limiting reads can
#' stop decompression early; failures may leave partial output files.
#'
#' @returns The `output` path, invisibly. Raises an R error if the Bash utility
#'   cannot be run or reports failure.
#' @seealso [relabel_fastq()], [decode_fastq_triplets()]
#' @export
#'
#' @examples
#' \dontrun{
#' compress_fastq_sequences(
#'   "input.fastq.gz", "reads.L98.seq.gz", read_length = 98L, threads = 8L
#' )
#' compress_fastq_sequences(
#'   "input.fastq.gz", "reads.L98.seq3.gz",
#'   read_length = 98L, triplets = TRUE, threads = 8L
#' )
#' }
compress_fastq_sequences <- function(input, output, read_length = 98L,
                                     triplets = FALSE, level = 6L, reads = NULL,
                                     threads = NULL, script = NULL) {
  args <- c(
    .fastq_option("--read-length", read_length),
    if (isTRUE(triplets)) "--triplets",
    .fastq_option("--level", level),
    .fastq_option("--reads", reads),
    .fastq_option("--threads", threads)
  )
  .run_fastq_bash("compress_fastq_sequences", input, output, args, script)
}


#' Restore FASTQ from compressed concatenated sequences
#'
#' Decompress a sequence-only file, optionally decode its triplets, and split
#' the bases into reads of a fixed length. Generate sequential read names,
#' literal `+` separator lines, and constant `I` quality strings.
#'
#' @param input Character scalar giving the path to a gzip file created by
#'   [compress_fastq_sequences()], with or without triplet encoding.
#' @param output Character scalar giving the destination FASTQ path. It must
#'   differ from `input`. A suffix of `.gz` selects gzip output at level 6;
#'   other filenames receive plain FASTQ. Existing output is overwritten.
#' @param read_length Positive integer giving the original fixed read length.
#'   Default: 98. Must match the value used during compression.
#' @param triplets Logical; decode the custom triplet encoding before splitting
#'   reads. Default: `FALSE`. Must match the compression setting.
#' @param script Optional character path to `fastq_compression.bash`. When
#'   `NULL`, use `getOption("igsc.fastq_script")` if set; otherwise locate the
#'   script under `extdata/compress_fastq` in the installed `igsc` package.
#'
#' @details
#' Bash and Python 3 must be available on the R session's `PATH`. The Python
#' decoder uses only standard-library modules. Decoding, formatting, and optional
#' output compression run serially with bounded processing chunks.
#'
#' Each recovered record has a header such as `@read1`, the recovered sequence,
#' a literal `+` line, and an `I` quality string matching `read_length`.
#' Names restart at one on every call. Sequence order and case are retained;
#' original headers and quality scores cannot be recovered from this format.
#'
#' Encoding mode and read length are not detected automatically. Incorrect
#' settings may produce incorrect sequences or read boundaries without an error.
#' Invalid triplet symbols or escapes, line breaks in sequence data, and a final
#' incomplete read cause an error. Empty input produces an empty FASTQ file.
#' Errors are written to stderr, and failures may leave partial output files.
#'
#' @returns The `output` path, invisibly. Raises an R error if the Bash utility
#'   cannot be run or reports failure.
#' @seealso [compress_fastq_sequences()], [relabel_fastq()]
#' @export
#'
#' @examples
#' \dontrun{
#' decode_fastq_triplets(
#'   "reads.L98.seq.gz", "recovered.fastq", read_length = 98L
#' )
#' decode_fastq_triplets(
#'   "reads.L98.seq3.gz", "recovered.fastq.gz",
#'   read_length = 98L, triplets = TRUE
#' )
#' }
decode_fastq_triplets <- function(input, output, read_length = 98L,
                                  triplets = FALSE, script = NULL) {
  args <- c(
    .fastq_option("--read-length", read_length),
    if (isTRUE(triplets)) "--triplets" else "--no-triplets"
  )
  .run_fastq_bash("decode_fastq_triplets", input, output, args, script)
}
