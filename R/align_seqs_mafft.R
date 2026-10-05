#' Align a sequence set with MAFFT
#'
#' Run the external MAFFT executable and return aligned Biostrings sequences.
#' The input, processor, and verbosity arguments follow [DECIPHER::AlignSeqs()].
#' MAFFT has different alignment algorithms and scoring conventions, so this
#' is not a replacement for DECIPHER-specific tuning arguments.
#'
#' @param myXStringSet A non-empty `DNAStringSet`, `RNAStringSet`, or
#'   `AAStringSet`. Existing gap characters (`-`) are removed before alignment.
#'   Every sequence must contain at least one non-gap character.
#' @param processors Positive integer thread count, or `NULL` to let MAFFT
#'   detect and use all available processors.
#' @param verbose Logical. Print a starting message and MAFFT's diagnostic log
#'   after the process finishes. Failure diagnostics are included in errors
#'   regardless of this setting.
#' @param ... Reserved for detecting unsupported arguments. DECIPHER arguments
#'   such as `guideTree`, `iterations`, `refinements`, and `gapOpening` are not
#'   silently translated or ignored; supplying them raises an error.
#' @param strategy MAFFT strategy: `"auto"`, `"linsi"`, `"ginsi"`, `"einsi"`,
#'   or `"fftns2"`. The default lets MAFFT choose a strategy.
#' @param mafft Optional executable path or command name. By default, search
#'   `PATH`, then common MacPorts and Homebrew locations, including
#'   `/opt/local/bin/mafft`. MAFFT must already be installed.
#' @param mafft_args Character vector of additional command-line arguments,
#'   with flags and their values as separate entries, for example
#'   `c("--op", "2", "--ep", "0.1")`. Use MAFFT's own scoring units.
#'   Options that change the input set, sequence direction, sequence type,
#'   thread count, or output format are not supported. Output order is restored
#'   even if `"--reorder"` is supplied. See the sections below for common options,
#'   strategy selection, and arguments managed by the wrapper.
#'
#' @return An aligned sequence set of the same type as `myXStringSet`, with
#'   equal sequence widths and the original sequence order and names, including
#'   duplicate names or absent names. DNA returns a `DNAStringSet`, protein
#'   returns an `AAStringSet`, and RNA returns an `RNAStringSet`.
#' @details Temporary FASTA files use generated identifiers so that spaces,
#'   punctuation, and duplicate sequence names are preserved on return. All
#'   temporary files are removed on success or failure. A single sequence is
#'   returned without running MAFFT, after removing existing gaps.
#'
#'   This performs a new multiple-sequence alignment of all input sequences.
#'   It does not implement reference-guided `--addfragments` alignment.
#'
#' @section Choosing a strategy:
#' - `"auto"`: let MAFFT choose based on the data size; a useful starting point.
#' - `"linsi"`: local pairwise information plus up to 1000 refinement cycles.
#'   Useful when conserved regions are separated by variable regions.
#' - `"ginsi"`: global pairwise information plus up to 1000 refinement cycles.
#'   Useful for sequences of similar lengths that align across their full length.
#' - `"einsi"`: generalized-affine pairwise alignment, up to 1000 refinement
#'   cycles, and `--ep 0`. Intended for sequences with large unalignable regions.
#' - `"fftns2"`: fast progressive alignment, two guide-tree calculations, and
#'   no iterative refinement. Useful when speed is more important.
#'
#' The three refinement strategies can be much slower on large sequence sets.
#' These names select MAFFT algorithms; they do not reproduce DECIPHER results.
#'
#' @section Common MAFFT arguments:
#' Supply each option and value separately in `mafft_args`. Values are strings;
#' options without values are single entries, for example `"--memsave"`.
#'
#' - `--op`: gap-opening penalty. Larger positive values discourage new gaps;
#'   try `c("--op", "3")` to increase the usual value of 1.53.
#' - `--ep`: an offset that behaves like a gap-extension penalty. Increasing it
#'   generally discourages extending gaps, but it is not a conventional
#'   per-residue extension cost. For example, `c("--ep", "0.1")`.
#' - `--maxiterate`: maximum iterative-refinement cycles. More cycles can
#'   increase runtime; zero disables refinement. For example,
#'   `c("--maxiterate", "100")`.
#' - `--retree`: number of guide-tree calculations in progressive alignment.
#'   For the `"fftns2"` strategy, `c("--retree", "1")` provides a faster,
#'   less thoroughly refined guide tree than the preset value of two.
#' - `--bl`: protein BLOSUM matrix number, for example `c("--bl", "62")`.
#'   Applies to amino-acid sequences, not nucleotide scoring.
#' - `--memsave`: use MAFFT's memory-saving alignment algorithm for long
#'   sequences; this can trade speed for lower memory use.
#' - `--nofft`: disable the FFT approximation in group-to-group alignment.
#'
#' Defaults can differ between MAFFT versions and strategies. For example,
#' MAFFT 7.526 reports `--ep 0.0` in its command-line help, while older manuals
#' list 0.123. Consult the installed executable's `--help` output.
#'
#' Extra arguments follow the strategy preset on the command line. For manual
#' control of `--maxiterate` or `--retree`, choose a non-`"auto"` strategy so
#' automatic strategy selection does not decide those settings for you.
#' Gap costs use MAFFT's own scoring scale, not DECIPHER's negative gap scores.
#' Stronger gap penalties may favor mismatches and do not guarantee a gap-free
#' alignment or subject sequence.
#'
#' @section Arguments managed by this wrapper:
#' Set `processors` instead of passing `--thread`; `processors = NULL` passes
#' `--thread -1`. Set `verbose = FALSE` to add `--quiet`. The input class selects
#' `--nuc` or `--amino`, and `--anysymbol` is supplied automatically.
#' `--reorder` cannot change the returned sequence order: original order is
#' restored before returning. FASTA output is required. Options such as
#' `--addfragments`, `--keeplength`, `--adjustdirection`, and `--clustalout`
#' are rejected because they conflict with this wrapper's input/output contract.
#'
#' @seealso [DECIPHER::AlignSeqs()], [aln_plot()],
#'   \url{https://mafft.cbrc.jp/alignment/software/},
#'   \url{https://mafft.cbrc.jp/alignment/software/manual/manual.html}
#' @export
#' @examples
#' \dontrun{
#' dna <- Biostrings::DNAStringSet(c(
#'   first = "AAAACCCCGGGG", second = "AAAATCCCCGGGG", third = "AAAACCCGGGG"
#' ))
#' aligned <- align_seqs_mafft(dna, processors = 4, verbose = FALSE)
#' aln_plot(aligned)
#'
#' # Make new gaps more costly; also discourage gap extension.
#' fewer_gaps <- align_seqs_mafft(
#'   dna, verbose = FALSE, mafft_args = c("--op", "3", "--ep", "0.1")
#' )
#'
#' # Start from L-INS-i but reduce its refinement limit from 1000 to 100.
#' limited_refinement <- align_seqs_mafft(
#'   dna, strategy = "linsi", verbose = FALSE,
#'   mafft_args = c("--maxiterate", "100")
#' )
#'
#' proteins <- Biostrings::AAStringSet(c(first = "MKVLW", second = "MKVLAW"))
#' aligned_aa <- align_seqs_mafft(proteins, strategy = "linsi",
#'                                mafft = "/opt/local/bin/mafft", verbose = FALSE)
#'
#' # Choose a protein substitution matrix in MAFFT's own scoring system.
#' aligned_blosum62 <- align_seqs_mafft(
#'   proteins, strategy = "ginsi", verbose = FALSE,
#'   mafft_args = c("--bl", "62")
#' )
#' }
align_seqs_mafft <- function(myXStringSet,
                              processors = 1,
                              verbose = TRUE,
                              ...,
                              strategy = c("auto", "linsi", "ginsi", "einsi", "fftns2"),
                              mafft = NULL,
                              mafft_args = character()) {
  extra <- as.list(substitute(list(...)))[-1L]
  if (length(extra)) {
    labels <- names(extra)
    if (is.null(labels)) labels <- rep("", length(extra))
    labels[is.na(labels) | !nzchar(labels)] <- "<unnamed>"
    stop("Unsupported argument(s): ", paste(labels, collapse = ", "),
         ". DECIPHER-specific options cannot be translated directly; ",
         "use strategy and mafft_args instead.", call. = FALSE)
  }
  strategy <- match.arg(strategy)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(processors) &&
      (!is.numeric(processors) || length(processors) != 1L ||
       !is.finite(processors) || processors < 1 || processors != trunc(processors) ||
       processors > .Machine$integer.max)) {
    stop("processors must be a positive integer or NULL.", call. = FALSE)
  }
  if (!is.character(mafft_args) || anyNA(mafft_args) || any(!nzchar(mafft_args))) {
    stop("mafft_args must be a character vector of separate flags and values.", call. = FALSE)
  }
  flags <- sub("=.*$", "", mafft_args)
  unsupported <- flags[grepl("^--(add|seed|adjustdirection|thread|out|clustal|phylip|nuc$|amino$|keeplength$|mapout$|compactmapout$|dash$)", flags)]
  if (length(unsupported)) {
    stop("Unsupported mafft_args option(s): ", paste(unique(unsupported), collapse = ", "),
         ". This wrapper preserves the input sequences and returns FASTA-based alignment data.",
         call. = FALSE)
  }
  igsc:::.ensure_package("Biostrings")
  type <- if (methods::is(myXStringSet, "DNAStringSet")) "DNA" else
    if (methods::is(myXStringSet, "RNAStringSet")) "RNA" else
      if (methods::is(myXStringSet, "AAStringSet")) "AA" else NULL
  if (is.null(type)) {
    stop("myXStringSet must be a DNAStringSet, RNAStringSet, or AAStringSet.", call. = FALSE)
  }
  if (!length(myXStringSet)) stop("myXStringSet must not be empty.", call. = FALSE)
  constructor <- switch(type, DNA = Biostrings::DNAStringSet,
                         RNA = Biostrings::RNAStringSet, AA = Biostrings::AAStringSet)
  reader <- switch(type, DNA = Biostrings::readDNAStringSet,
                    RNA = Biostrings::readRNAStringSet, AA = Biostrings::readAAStringSet)
  original_names <- names(myXStringSet)
  sequences <- unname(gsub("-", "", as.character(myXStringSet), fixed = TRUE))
  if (any(!nzchar(sequences))) {
    stop("Every input sequence must contain at least one non-gap character.", call. = FALSE)
  }
  if (length(sequences) == 1L) {
    result <- constructor(sequences)
    names(result) <- original_names
    return(result)
  }
  executable <- .find_mafft_executable(mafft)
  folder <- tempfile("igsc-mafft-")
  if (!dir.create(folder)) stop("Could not create a MAFFT temporary directory.", call. = FALSE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  input <- file.path(folder, "input.fasta")
  output <- file.path(folder, "aligned.fasta")
  logfile <- file.path(folder, "mafft.log")
  ids <- paste0("igsc_seq_", seq_along(sequences))
  Biostrings::writeXStringSet(constructor(stats::setNames(sequences, ids)), input)
  algorithm <- switch(strategy,
    auto = "--auto",
    linsi = c("--localpair", "--maxiterate", "1000"),
    ginsi = c("--globalpair", "--maxiterate", "1000"),
    einsi = c("--genafpair", "--maxiterate", "1000", "--ep", "0"),
    fftns2 = c("--retree", "2", "--maxiterate", "0")
  )
  args <- c(algorithm, "--thread", if (is.null(processors)) "-1" else as.character(processors),
            if (type == "AA") "--amino" else "--nuc", "--anysymbol", "--inputorder",
            if (!verbose) "--quiet", mafft_args, input)
  if (verbose) message("Aligning ", length(sequences), " sequences with MAFFT (", strategy, ").")
  status <- tryCatch(
    system2(executable, args = shQuote(args), stdout = output, stderr = logfile),
    error = function(error) stop("Could not run MAFFT: ", conditionMessage(error), call. = FALSE)
  )
  diagnostics <- if (file.exists(logfile)) readLines(logfile, warn = FALSE) else character()
  if (status != 0L || !file.exists(output) || file.info(output)$size == 0) {
    stop("MAFFT failed (exit status ", status, ").\n",
         paste(utils::tail(diagnostics, 30L), collapse = "\n"), call. = FALSE)
  }
  if (verbose && length(diagnostics)) message(paste(diagnostics, collapse = "\n"))
  result <- reader(output)
  if (length(result) != length(ids) || anyDuplicated(names(result)) ||
      !setequal(names(result), ids)) {
    stop("MAFFT output does not contain exactly the expected sequence identifiers.", call. = FALSE)
  }
  result <- result[match(ids, names(result))]
  if (length(unique(nchar(as.character(result)))) != 1L ||
      !identical(unname(gsub("-", "", as.character(result), fixed = TRUE)), sequences)) {
    stop("MAFFT output has unequal widths or changed non-gap sequence content.", call. = FALSE)
  }
  names(result) <- original_names
  result
}

.find_mafft_executable <- function(mafft = NULL) {
  if (is.null(mafft)) {
    candidates <- unique(c(unname(Sys.which("mafft")), "/opt/local/bin/mafft",
                           "/opt/homebrew/bin/mafft", "/usr/local/bin/mafft"))
  } else {
    if (!is.character(mafft) || length(mafft) != 1L || is.na(mafft) || !nzchar(mafft)) {
      stop("mafft must be an executable path or command name.", call. = FALSE)
    }
    candidates <- unique(c(path.expand(mafft), unname(Sys.which(mafft))))
  }
  available <- candidates[nzchar(candidates) & file.exists(candidates) &
                            !dir.exists(candidates) & file.access(candidates, 1) == 0]
  if (!length(available)) {
    stop("MAFFT executable not found. Supply mafft = '/path/to/mafft' ",
         "or install MAFFT and add it to PATH.", call. = FALSE)
  }
  normalizePath(available[1], mustWork = TRUE)
}
