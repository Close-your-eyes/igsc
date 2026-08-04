#' Calculate TCRdist3 distance matrices
#'
#' Runs the bundled Python script to calculate a full TCRdist distance
#' matrix and a CDR3-only distance matrix for either T-cell receptor alpha
#' (TRA) or beta (TRB) sequences.
#' https://github.com/kmayerb/tcrdist3
#' pip install tcrdist3 pandas
#'
#' This function works with csv files on disk. See run_tcrdist3_reticulate as
#' alternative.
#'
#' @param input Character scalar giving the path to the input CSV file.
#'   A relative path is interpreted relative to `workdir`.
#' @param script Character scalar giving the path to the Python script that
#'   performs the TCRdist3 calculation. By default, the script bundled with
#'   the \pkg{igsc} package is used.
#' @param db_file Character scalar giving the path to the TCRdist3 reference
#'   database file. A relative path is interpreted relative to `workdir`.
#' @param chain Character scalar specifying the receptor chain to process.
#'   Must be `"TRB"` or `"TRA"`. The default is `"TRB"`.
#' @param species Character scalar specifying the organism. Must be
#'   `"human"` or `"mouse"`. Gene names in `input` must be appropriate for
#'   the selected species.
#' @param suffix Character scalar containing an optional suffix to append to
#'   the output filenames. For example, `suffix = "sample1"` produces a
#'   filename ending in `_sample1.csv`. The default, `""`, adds no suffix.
#' @param workdir Character scalar giving the working directory used to
#'   resolve relative input paths and store output files. The directory must
#'   already exist.
#' @param python Character scalar giving the path to a Python 3 executable.
#'   The Python environment must contain the `pandas` and `tcrdist3`
#'   packages.
#'
#' @details
#' The input must be a comma-separated text file with a header row. Required
#' columns depend on `chain`:
#'
#' \tabular{llll}{
#'   Chain \tab CDR3 amino acid sequence \tab V gene \tab J gene \cr
#'   `"TRA"` \tab `cdr3_a_aa` \tab `v_a_gene` \tab `j_a_gene` \cr
#'   `"TRB"` \tab `cdr3_b_aa` \tab `v_b_gene` \tab `j_b_gene`
#' }
#'
#' Rows with missing values in any required column are removed. CDR3 amino
#' acid sequences are converted to uppercase, and surrounding whitespace is
#' removed from the CDR3, V-gene, and J-gene columns.
#'
#' The input may optionally contain a `count` column giving the abundance of
#' each receptor. If it is absent, every row is assigned a count of one.
#' Additional columns are allowed and are passed to TCRdist3 unchanged.
#'
#' V- and J-gene names must use nomenclature recognized by TCRdist3 and must
#' be compatible with `species` and `db_file`.
#'
#' Two CSV files are written to `workdir`:
#'
#' * `tcrdist3_full_distmat_alpha[_suffix].csv` or
#'   `tcrdist3_full_distmat_beta[_suffix].csv`
#' * `tcrdist3_cdr3_distmat_alpha[_suffix].csv` or
#'   `tcrdist3_cdr3_distmat_beta[_suffix].csv`
#'
#' Existing files with the same names are overwritten.
#'
#' @return Invisibly returns a named character vector containing the paths
#'   to the two output files. Element `full` is the full TCRdist matrix and
#'   element `cdr3` is the CDR3-only distance matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process a TRB repertoire using the default filenames.
#' files <- run_tcrdist3(
#'   workdir = "/path/to/analysis",
#'   python = "/opt/local/bin/python3"
#' )
#'
#' # Process a TRA repertoire and add a suffix to the output files.
#' files <- run_tcrdist3(
#'   chain = "TRA",
#'   input = "tra_repertoire.csv",
#'   species = "human",
#'   suffix = "sample1",
#'   workdir = "/path/to/analysis",
#'   python = "/opt/local/bin/python3"
#' )
#'
#' # Read the generated matrices into R.
#' full_distance <- read.csv(
#'   files["full"],
#'   row.names = 1,
#'   check.names = FALSE
#' )
#'
#' cdr3_distance <- read.csv(
#'   files["cdr3"],
#'   row.names = 1,
#'   check.names = FALSE
#' )
run_tcrdist3 <- function(
    input,
    script = system.file("extdata", "run_tcrdist3.py", package = "igsc"),
    db_file = "alphabeta_gammadelta_db.tsv",
    chain = c("TRB", "TRA"),
    species = c("human", "mouse"),
    suffix = "",
    workdir = getwd(),
    python = "/opt/local/bin/python3"
) {
  chain <- match.arg(chain)
  species <- match.arg(species)

  # Normalize paths that must already exist.
  script <- normalizePath(
    script,
    mustWork = TRUE
  )

  workdir <- normalizePath(
    workdir,
    mustWork = TRUE
  )

  if (!file.exists(python)) {
    stop("Python executable does not exist: ", python)
  }

  if (!dir.exists(workdir)) {
    stop("Working directory does not exist: ", workdir)
  }

  args <- c(
    script,
    "--db-file", db_file,
    "--chain", chain,
    "--input", input,
    "--species", species,
    "--suffix", suffix,
    "--workdir", workdir
  )

  # Display the command before running it.
  command_text <- paste(
    shQuote(python),
    paste(shQuote(args), collapse = " ")
  )

  message("Running:\n", command_text)

  status <- system2(
    command = python,
    args = shQuote(args),
    stdout = "",
    stderr = ""
  )

  if (status != 0L) {
    stop("The tcrdist3 script failed with exit status ", status)
  }

  output_prefix <- if (chain == "TRA") "alpha" else "beta"
  suffix_part <- if (nzchar(suffix)) paste0("_", suffix) else ""

  output_files <- c(
    full = file.path(
      workdir,
      paste0(
        "tcrdist3_full_distmat_",
        output_prefix,
        suffix_part,
        ".csv"
      )
    ),
    cdr3 = file.path(
      workdir,
      paste0(
        "tcrdist3_cdr3_distmat_",
        output_prefix,
        suffix_part,
        ".csv"
      )
    )
  )

  message(
    "Finished successfully.\nOutput files:\n",
    paste(output_files, collapse = "\n")
  )

  invisible(output_files)
}




#' Calculate T-cell receptor distance matrices
#'
#' Selects the standard Python implementation for data sets containing
#' at most 10,000 rows and the explicit dense implementation for larger
#' data sets.
#'
#' Run these lines before calling the function:
#' library(reticulate)
#' use_python("/opt/local/bin/python3",required = TRUE)
#' library(igsc)
#'
#' @param data A data frame containing TCR sequences and gene calls.
#' @param chain TCR chain: `"TRA"` or `"TRB"`.
#' @param species Organism: `"human"` or `"mouse"`.
#' @param db_file Name or path of the tcrdist3 reference database.
#' @param cpus Number of CPUs used for the explicit dense calculation.
#' @param dense_threshold Number of rows above which the explicit dense
#'   Python function is used.
#'
#' @return A `tcrdist_result` object.
#' @seealso
#' [run_tcrdist3()]
#'
#' @export
run_tcrdist3_reticulate <- function(
    data,
    chain = c("TRB", "TRA"),
    species = c("human", "mouse"),
    db_file = c("alphabeta_gammadelta_db.tsv", "alphabeta_db.tsv", "gammadelta_db.tsv"),
    cpus = 4L,
    dense_threshold = 10000L
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  chain <- match.arg(chain)
  species <- match.arg(species)
  db_file <- match.arg(db_file)

  if (
    !is.character(db_file) ||
    length(db_file) != 1L ||
    is.na(db_file)
  ) {
    stop(
      "`db_file` must be a single non-missing character value.",
      call. = FALSE
    )
  }

  cpus <- as.integer(cpus)

  if (
    length(cpus) != 1L ||
    is.na(cpus) ||
    cpus < 1L
  ) {
    stop("`cpus` must be a positive integer.", call. = FALSE)
  }

  dense_threshold <- as.integer(dense_threshold)

  if (
    length(dense_threshold) != 1L ||
    is.na(dense_threshold) ||
    dense_threshold < 1L
  ) {
    stop(
      "`dense_threshold` must be a positive integer.",
      call. = FALSE
    )
  }

  module <- get_tcrdist_module()
  number_of_rows <- nrow(data)

  if (number_of_rows > dense_threshold) {
    message(
      sprintf(
        paste0(
          "Input has %s rows, exceeding the threshold of %s. ",
          "Using Python calculate_tcr_distances2() with an ",
          "explicit dense distance calculation."
        ),
        format(number_of_rows, big.mark = ","),
        format(dense_threshold, big.mark = ",")
      )
    )

    result <- module$calculate_tcr_distances2(
      data = data,
      chain = chain,
      species = species,
      db_file = db_file,
      cpus = cpus
    )
  } else {
    message(
      sprintf(
        paste0(
          "Input has %s rows, within the threshold of %s. ",
          "Using Python calculate_tcr_distances()."
        ),
        format(number_of_rows, big.mark = ","),
        format(dense_threshold, big.mark = ",")
      )
    )

    result <- module$calculate_tcr_distances(
      data = data,
      chain = chain,
      species = species,
      db_file = db_file
    )
  }

  labels <- as.character(unlist(
    result$labels,
    use.names = FALSE
  ))

  full_distance <- as.matrix(result$full_distance)
  cdr3_distance <- as.matrix(result$cdr3_distance)

  expected_dimension <- c(
    length(labels),
    length(labels)
  )

  if (!identical(dim(full_distance), expected_dimension)) {
    stop(
      "Python returned a full-distance matrix with unexpected dimensions.",
      call. = FALSE
    )
  }

  if (!identical(dim(cdr3_distance), expected_dimension)) {
    stop(
      "Python returned a CDR3-distance matrix with unexpected dimensions.",
      call. = FALSE
    )
  }

  dimnames(full_distance) <- list(labels, labels)
  dimnames(cdr3_distance) <- list(labels, labels)

  structure(
    list(
      full_distance = full_distance,
      cdr3_distance = cdr3_distance,
      labels = labels,
      clone_data = as.data.frame(result$clone_data),
      chain = chain,
      species = species,
      cpus = cpus,
      python_function = if (number_of_rows > dense_threshold) {
        "calculate_tcr_distances2"
      } else {
        "calculate_tcr_distances"
      }
    ),
    class = "tcrdist_result"
  )
}
