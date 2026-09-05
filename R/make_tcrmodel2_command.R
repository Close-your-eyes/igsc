#' Create a TCRmodel2 command
#'
#' Creates a shell command for running TCRmodel2 with either an MHC class I
#' or MHC class II complex. The MHC class is inferred from the names present
#' in `seqs`.
#'
#' @param seqs A named character vector containing amino-acid sequences.
#'   All inputs require elements named `TCRa`, `TCRb`, and `pep`.
#'   For MHC class I, supply an element named `MHC`. For MHC class II,
#'   supply elements named `MHC1` and `MHC2`, representing the MHC alpha
#'   and beta chains, respectively.
#' @param job_id A character scalar giving the TCRmodel2 job identifier.
#' @param output_dir A character scalar giving the output directory.
#'   Defaults to `"experiments/"`.
#' @param ori_db A character scalar giving the path to the AlphaFold
#'   database directory containing the PDB mmCIF files and model parameters.
#' @param python A character scalar specifying the Python executable.
#'   Defaults to `"python"`.
#' @param script A character scalar specifying the path to
#'   `run_tcrmodel2.py`.
#'
#' @details
#' MHC class I input is translated to the command-line argument
#' `--mhca_seq`.
#'
#' For MHC class II, `MHC1` is translated to `--mhca_seq` and `MHC2`
#' is translated to `--mhcb_seq`.
#'
#' Sequence letters are converted to uppercase, and whitespace within
#' sequences is removed. The function validates the presence and format of
#' the sequences but does not verify that they contain the domains required
#' by TCRmodel2.
#'
#' @return A single character string containing a multiline shell command.
#'
#' @examples
#' class_I <- c(
#'   TCRa = "CAVRDSNYQLIW",
#'   TCRb = "CASSLGQETQYF",
#'   pep  = "GILGFVFTL",
#'   MHC  = "SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQF"
#' )
#'
#' cmd_I <- make_tcrmodel2_command(
#'   class_I,
#'   job_id = "class_I_example",
#'   ori_db = "/path/to/alphafold/data"
#' )
#' cat(cmd_I)
#'
#' class_II <- c(
#'   TCRa = "CALITGGGNKLTF",
#'   TCRb = "CASRLQGWNSPLHF",
#'   pep  = "EVVRHCPHHER",
#'   MHC1 = "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEI",
#'   MHC2 = "TRPRFLEYSTSECHFFNGTERVRFLDRYFHNQ"
#' )
#'
#' cmd_II <- make_tcrmodel2_command(
#'   class_II,
#'   job_id = "class_II_example",
#'   ori_db = "/path/to/alphafold/data"
#' )
#' cat(cmd_II)
#'
#' \dontrun{
#' system(cmd_II)
#' }
#'
#' @export
make_tcrmodel2_command <- function(
    seqs,
    job_id,
    output_dir = "experiments/",
    ori_db = "/sc-scratch/sc-scratch-cc13-scurbi/_bin/tcrmodel2/alphafolddl",
    python = "python",
    script = "run_tcrmodel2.py"
) {
  if (is.null(names(seqs))) {
    stop("`seqs` must be a named vector.")
  }

  if (anyDuplicated(names(seqs))) {
    stop("Sequence names must be unique.")
  }

  has_class_I  <- "MHC" %in% names(seqs)
  has_class_II <- all(c("MHC1", "MHC2") %in% names(seqs))
  has_partial_II <- any(c("MHC1", "MHC2") %in% names(seqs))

  if (has_class_I && has_partial_II) {
    stop("Supply either `MHC` for class I or `MHC1` and `MHC2` for class II.")
  }

  if (has_partial_II && !has_class_II) {
    stop("Class II input requires both `MHC1` and `MHC2`.")
  }

  if (!has_class_I && !has_class_II) {
    stop("Supply `MHC` for class I or `MHC1` and `MHC2` for class II.")
  }

  mhc_class <- if (has_class_II) "II" else "I"

  required <- if (mhc_class == "I") {
    c("TCRa", "TCRb", "pep", "MHC")
  } else {
    c("TCRa", "TCRb", "pep", "MHC1", "MHC2")
  }

  missing <- setdiff(required, names(seqs))

  if (length(missing)) {
    stop("Missing sequence(s): ", paste(missing, collapse = ", "))
  }

  # Normalize case and remove spaces or line breaks
  seqs <- setNames(
    toupper(gsub("\\s+", "", as.character(seqs))),
    names(seqs)
  )

  if (anyNA(seqs[required]) || any(seqs[required] == "")) {
    stop("Sequences cannot be missing or empty.")
  }

  if (any(!grepl("^[A-Z]+$", seqs[required]))) {
    stop("Sequences may contain amino-acid letters only.")
  }

  quote_arg <- function(x) {
    x <- as.character(x)

    if (grepl("^[A-Za-z0-9_./:+-]+$", x)) {
      x
    } else {
      shQuote(x)
    }
  }

  mhc_args <- if (mhc_class == "I") {
    sprintf("--mhca_seq=%s", quote_arg(seqs[["MHC"]]))
  } else {
    c(
      sprintf("--mhca_seq=%s", quote_arg(seqs[["MHC1"]])),
      sprintf("--mhcb_seq=%s", quote_arg(seqs[["MHC2"]]))
    )
  }

  command_parts <- c(
    sprintf("%s %s", quote_arg(python), quote_arg(script)),
    sprintf("--job_id=%s", quote_arg(job_id)),
    sprintf("--output_dir=%s", quote_arg(output_dir)),
    sprintf("--tcra_seq=%s", quote_arg(seqs[["TCRa"]])),
    sprintf("--tcrb_seq=%s", quote_arg(seqs[["TCRb"]])),
    sprintf("--pep_seq=%s", quote_arg(seqs[["pep"]])),
    mhc_args,
    sprintf("--ori_db=%s", quote_arg(ori_db))
  )

  paste0(paste(command_parts, collapse = " \\\n"), "\n")
}
