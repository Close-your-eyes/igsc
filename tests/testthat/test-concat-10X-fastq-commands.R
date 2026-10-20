test_that("commands concatenate gzip inputs in deterministic order with quoted paths", {
  root <- tempfile("fastq commands '")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  input_dir <- file.path(root, "input $reads")
  save_path <- file.path(root, "output;reads")
  dir.create(input_dir)
  sample <- "sample_with.dots"
  filenames <- paste0(sample, c(
    "_S2_L002_R1_001.fastq.gz",
    "_S1_L002_R1_002.fastq.gz",
    "_S1_L001_R1_001.fastq.gz",
    "_S1_L002_R1_001.fastq.gz"
  ))
  inputs <- file.path(input_dir, filenames)
  for (i in seq_along(inputs)) {
    connection <- gzfile(inputs[i], "wt")
    writeLines(paste0("read", i), connection)
    close(connection)
  }

  commands <- concat_10X_fastq_commands(inputs, save_path)
  expect_length(commands, 1L)
  expect_null(names(commands))
  expect_identical(commands, concat_10X_fastq_commands(rev(inputs), save_path))
  expect_false(dir.exists(save_path))
  dir.create(file.path(save_path, sample), recursive = TRUE)
  expect_identical(system(commands), 0L)

  output <- file.path(save_path, sample,
                      paste0(sample, "_S1_L001_R1_001.fastq.gz"))
  connection <- gzfile(output, "rt")
  expect_identical(readLines(connection), c("read3", "read4", "read2", "read1"))
  close(connection)
  expect_equal(file.info(output)$size, sum(file.info(inputs)$size))
})

test_that("nonlocal paths are grouped by sample and read type", {
  inputs <- c(
    "remote/sample_S9_L002_R2_001.fastq.gz",
    "remote/other_S1_L001_I1_001.fastq.gz",
    "remote/sample_S1_L001_R1_001.fastq.gz",
    "remote/sample_S1_L001_R2_001.fastq.gz",
    "remote/other_S1_L001_I2_001.fastq.gz"
  )
  commands <- concat_10X_fastq_commands(inputs, "output")
  expect_length(commands, 4L)
  expect_true(any(grepl("other_S1_L001_I1_001.fastq.gz", commands, fixed = TRUE)))
  expect_true(any(grepl("other_S1_L001_I2_001.fastq.gz", commands, fixed = TRUE)))
  expect_identical(commands, concat_10X_fastq_commands(rev(inputs), "output"))
  expect_identical(concat_10X_fastq_commands(character()), character())
})

test_that("ambiguous or invalid inputs are rejected", {
  input <- "sample_S1_L001_R1_001.fastq.gz"
  expect_error(concat_10X_fastq_commands(NA_character_), "non-NA")
  expect_error(concat_10X_fastq_commands(1), "character vector")
  expect_error(concat_10X_fastq_commands(""), "non-empty")
  expect_error(concat_10X_fastq_commands("sample.fastq.gz"), "naming scheme")
  expect_error(concat_10X_fastq_commands(c(input, file.path("copy", input))),
               "Duplicate")
  expect_error(concat_10X_fastq_commands(input, NA_character_), "save_path")
  expect_error(concat_10X_fastq_commands(file.path("out", "sample", input), "out"),
               "overwrite")
})
