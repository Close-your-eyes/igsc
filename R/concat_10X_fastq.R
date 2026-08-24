#' Concatenate lane-split 10x FASTQ files
#'
#' Concatenates gzip-compressed FASTQ files belonging to the same sample and
#' read type. Concatenating gzip streams directly is valid: the resulting file
#' is a multi-member gzip archive.
#'
#' Expected input filenames follow this structure:
#' `sample_S1_L001_R1_001.fastq.gz`
#'
#' @param dir Directory searched recursively for input FASTQ files.
#' @param save_path Parent directory for sample-specific output directories.
#' @param exec Logical; execute commands when `TRUE`, otherwise return a plan.
#' @param existing How existing output files are handled: `"error"`, `"skip"`,
#'   or `"overwrite"`.
#' @param create_save_path Logical; create `save_path` when it does not exist.
#' @param verify Logical; after execution, verify that output size equals the
#'   sum of input file sizes.
#'
#' @return A data frame containing output files, input files, commands, and
#'   execution status.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plan <- concat_10X_fastq(
#'   dir = "/path/to/fastq",
#'   save_path = "/path/to/output"
#' )
#'
#' result <- concat_10X_fastq(
#'   dir = "/path/to/fastq",
#'   save_path = "/path/to/output",
#'   exec = TRUE,
#'   existing = "error"
#' )
#' }
concat_10X_fastq <- function(
    dir,
    save_path,
    exec = FALSE,
    existing = c("error", "skip", "overwrite"),
    create_save_path = TRUE,
    verify = TRUE
) {
  existing <- match.arg(existing)

  # Argument validation -------------------------------------------------------
  check_scalar_string <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(sprintf("`%s` must be one non-empty, non-NA string.", name),
           call. = FALSE)
    }
  }

  check_scalar_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }

  check_scalar_string(dir, "dir")
  check_scalar_string(save_path, "save_path")
  check_scalar_logical(exec, "exec")
  check_scalar_logical(create_save_path, "create_save_path")
  check_scalar_logical(verify, "verify")

  dir <- path.expand(dir)
  save_path <- path.expand(save_path)

  if (!dir.exists(dir)) {
    stop("Input directory does not exist: ", dir, call. = FALSE)
  }

  input_info <- file.info(dir)
  if (is.na(input_info$isdir) || !input_info$isdir) {
    stop("`dir` is not a directory: ", dir, call. = FALSE)
  }

  if (file.access(dir, mode = 4L) != 0L) {
    stop("Input directory is not readable: ", dir, call. = FALSE)
  }

  if (!dir.exists(save_path)) {
    if (!create_save_path) {
      stop("Output directory does not exist: ", save_path, call. = FALSE)
    }

    if (!dir.create(save_path, recursive = TRUE, showWarnings = FALSE) &&
        !dir.exists(save_path)) {
      stop("Could not create output directory: ", save_path, call. = FALSE)
    }
  }

  output_info <- file.info(save_path)
  if (is.na(output_info$isdir) || !output_info$isdir) {
    stop("`save_path` is not a directory: ", save_path, call. = FALSE)
  }

  if (file.access(save_path, mode = 2L) != 0L) {
    stop("Output directory is not writable: ", save_path, call. = FALSE)
  }

  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  save_path <- normalizePath(save_path, winslash = "/", mustWork = TRUE)

  # Locate and parse inputs ---------------------------------------------------
  fastq_pattern <- "\\.fastq\\.gz$"

  input_files <- list.files(
    path = dir,
    pattern = fastq_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = FALSE
  )

  if (!length(input_files)) {
    stop("No files ending in `.fastq.gz` were found under: ", dir,
         call. = FALSE)
  }

  input_files <- normalizePath(
    input_files,
    winslash = "/",
    mustWork = TRUE
  )

  input_file_info <- file.info(input_files)

  if (any(input_file_info$isdir %in% TRUE)) {
    stop("One or more detected FASTQ paths are directories.", call. = FALSE)
  }

  if (anyNA(input_file_info$size)) {
    stop("Could not determine the size of one or more input files.",
         call. = FALSE)
  }

  if (any(input_file_info$size == 0)) {
    empty_files <- input_files[input_file_info$size == 0]
    stop(
      "Empty input FASTQ file(s) detected:\n",
      paste(empty_files, collapse = "\n"),
      call. = FALSE
    )
  }

  if (any(file.access(input_files, mode = 4L) != 0L)) {
    unreadable <- input_files[file.access(input_files, mode = 4L) != 0L]
    stop(
      "Unreadable input FASTQ file(s):\n",
      paste(unreadable, collapse = "\n"),
      call. = FALSE
    )
  }

  # Capture:
  # 1: sample prefix, 2: sample number, 3: lane,
  # 4: read type, 5: chunk number
  filename_pattern <-
    "^(.*)_S([0-9]+)_L([0-9]{3})_([IR][12])_([0-9]{3})\\.fastq\\.gz$"

  filenames <- basename(input_files)
  matches <- regexec(filename_pattern, filenames, perl = TRUE)
  parsed <- regmatches(filenames, matches)
  valid <- lengths(parsed) == 6L

  if (!all(valid)) {
    stop(
      "FASTQ filename(s) do not match the expected 10x naming scheme ",
      "`<sample>_S<number>_L<number>_<I1|I2|R1|R2>_<chunk>.fastq.gz`:\n",
      paste(filenames[!valid], collapse = "\n"),
      call. = FALSE
    )
  }

  metadata <- data.frame(
    input = input_files,
    sample = vapply(parsed, `[[`, character(1), 2L),
    sample_number = as.integer(vapply(parsed, `[[`, character(1), 3L)),
    lane = as.integer(vapply(parsed, `[[`, character(1), 4L)),
    type = vapply(parsed, `[[`, character(1), 5L),
    chunk = as.integer(vapply(parsed, `[[`, character(1), 6L)),
    size = input_file_info$size,
    stringsAsFactors = FALSE
  )

  if (any(!nzchar(metadata$sample))) {
    stop("At least one FASTQ file has an empty sample prefix.",
         call. = FALSE)
  }

  # A duplicated sample/lane/type/chunk generally indicates accidentally
  # copied or ambiguously located FASTQ files.
  duplicate_key <- duplicated(
    metadata[c("sample", "sample_number", "lane", "type", "chunk")]
  )

  if (any(duplicate_key)) {
    duplicate_rows <- metadata[
      duplicated(
        metadata[c("sample", "sample_number", "lane", "type", "chunk")],
        fromLast = TRUE
      ) | duplicate_key,
      ,
      drop = FALSE
    ]

    stop(
      "Duplicate sample/lane/read/chunk combinations detected:\n",
      paste(duplicate_rows$input, collapse = "\n"),
      call. = FALSE
    )
  }

  # Ensure deterministic concatenation order.
  metadata <- metadata[
    order(
      metadata$sample,
      metadata$type,
      metadata$sample_number,
      metadata$lane,
      metadata$chunk,
      metadata$input
    ),
    ,
    drop = FALSE
  ]

  groups <- split(
    seq_len(nrow(metadata)),
    interaction(metadata$sample, metadata$type, drop = TRUE, lex.order = TRUE)
  )

  # Build output plan ---------------------------------------------------------
  plan_rows <- lapply(groups, function(index) {
    group <- metadata[index, , drop = FALSE]

    sample <- unique(group$sample)
    type <- unique(group$type)

    if (length(sample) != 1L || length(type) != 1L) {
      stop("Internal grouping error.", call. = FALSE)
    }

    output_dir <- file.path(save_path, sample)
    output <- file.path(
      output_dir,
      sprintf("%s_S1_L001_%s_001.fastq.gz", sample, type)
    )

    if (output %in% group$input) {
      stop(
        "An output file would overwrite one of its own input files: ",
        output,
        call. = FALSE
      )
    }

    command <- paste(
      "cat",
      paste(shQuote(group$input), collapse = " "),
      ">",
      shQuote(output)
    )

    data.frame(
      sample = sample,
      type = type,
      output_dir = output_dir,
      output = output,
      input_count = nrow(group),
      expected_bytes = sum(group$size),
      expected_size = brathering::format_bytes(sum(group$size)),
      inputs = I(list(group$input)),
      command = command,
      stringsAsFactors = FALSE
    )
  })

  plan <- do.call(rbind, plan_rows)
  rownames(plan) <- NULL

  plan$output_exists <- file.exists(plan$output)
  plan$status <- ifelse(plan$output_exists, "exists", "planned")

  if (anyDuplicated(plan$output)) {
    stop("Multiple groups resolve to the same output file.", call. = FALSE)
  }

  if (existing == "error" && any(plan$output_exists)) {
    stop(
      "Output file(s) already exist. Use `existing = \"skip\"` or ",
      "`existing = \"overwrite\"` if appropriate:\n",
      paste(plan$output[plan$output_exists], collapse = "\n"),
      call. = FALSE
    )
  }

  message("total size to write: ", brathering::format_bytes(sum(plan$expected_bytes)))
  if (!exec) {
    return(plan)
  }

  # Execute and verify --------------------------------------------------------
  for (i in seq_len(nrow(plan))) {
    output <- plan$output[i]
    output_dir <- plan$output_dir[i]

    if (file.exists(output) && existing == "skip") {
      plan$status[i] <- "skipped"
      next
    }

    if (!dir.exists(output_dir)) {
      created <- dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      if (!created && !dir.exists(output_dir)) {
        stop("Could not create output directory: ", output_dir,
             call. = FALSE)
      }
    }

    if (file.access(output_dir, mode = 2L) != 0L) {
      stop("Output directory is not writable: ", output_dir,
           call. = FALSE)
    }

    # Write to a temporary file first. This avoids leaving an apparently valid
    # final output when concatenation is interrupted or fails.
    temporary_output <- tempfile(
      pattern = paste0(".", basename(output), "."),
      tmpdir = output_dir
    )

    temporary_command <- paste(
      "cat",
      paste(shQuote(plan$inputs[[i]]), collapse = " "),
      ">",
      shQuote(temporary_output)
    )

    status <- system(temporary_command)

    if (!identical(status, 0L)) {
      if (file.exists(temporary_output)) {
        unlink(temporary_output)
      }

      plan$status[i] <- paste0("failed (exit ", status, ")")
      stop(
        "Concatenation failed for output:\n",
        output,
        "\nCommand exited with status ", status, ".",
        call. = FALSE
      )
    }

    if (!file.exists(temporary_output)) {
      stop("Command succeeded but did not create an output file: ",
           temporary_output, call. = FALSE)
    }

    actual_size <- file.info(temporary_output)$size

    if (is.na(actual_size) || actual_size == 0) {
      unlink(temporary_output)
      stop("Temporary output is missing or empty for: ", output,
           call. = FALSE)
    }

    if (verify && actual_size != plan$expected_bytes[i]) {
      unlink(temporary_output)
      stop(
        "Output-size verification failed for: ", output,
        "\nExpected ", plan$expected_bytes[i],
        " bytes, got ", actual_size, " bytes.",
        call. = FALSE
      )
    }

    # file.rename() does not reliably replace an existing file on every OS.
    if (file.exists(output)) {
      if (existing != "overwrite") {
        unlink(temporary_output)
        stop("Output appeared during execution: ", output, call. = FALSE)
      }

      if (!file.remove(output)) {
        unlink(temporary_output)
        stop("Could not remove existing output: ", output, call. = FALSE)
      }
    }

    if (!file.rename(temporary_output, output)) {
      unlink(temporary_output)
      stop("Could not move temporary output into place: ", output,
           call. = FALSE)
    }

    plan$status[i] <- "created"
    message(output)
  }

  plan$actual_size <- vapply(
    plan$output,
    function(path) {
      if (file.exists(path)) file.info(path)$size else NA_real_
    },
    numeric(1)
  )

  return(plan)
}

# original:
# concat_10X_fastq <- function(dir,
#                              save_path,
#                              exec = F) {
#
#   ff <- list.files(dir,
#                    full.names = T,
#                    recursive = T,
#                    pattern = "\\.fastq\\.gz")
#
#   files <- basename(ff)
#
#   prefixes <- sapply(strsplit(files, "_S[[:digit:]]{1,}_"), "[", 1)
#
#   ff2 <- split(ff, prefixes)
#   ff3 <- purrr::map(ff2, function(x) {
#     files2 <- basename(x)
#     types2 <- stringr::str_extract(files2, "[IR][12]")
#     split(x, types2)
#   })
#
#   save_paths <- file.path(save_path, names(ff3))
#   purrr::map(save_paths, dir.create)
#
#   cmdlst <- purrr::map2(ff3, save_paths, function(x,y) {
#     purrr::map_chr(x, function(z) {
#       prefix <- unique(sapply(strsplit(basename(z), "_S[[:digit:]]{1,}_"), "[", 1))
#       type <- unique(stringr::str_extract(basename(z), "[IR][12]"))
#       filename <- paste0(prefix, "_S1_L001_", type, "_001.fastq.gz")
#       paste0("cat ", paste(z, collapse = " "), " > ", file.path(y, filename))
#     })
#   })
#   cmds <- unlist(cmdlst)
#
#   if (exec) {
#     purrr::map(cmds, system)
#   }
#
#   return(cmds)
# }
