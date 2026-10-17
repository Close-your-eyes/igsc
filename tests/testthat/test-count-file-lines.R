check_text_line_counts <- function(directory) {
  contents <- c("", "one", "one\n", "\n", "one\n\nthree\n", "one\r\ntwo\r\n")
  counts <- c(0, 1, 1, 1, 3, 2)
  for (i in seq_along(contents)) {
    path <- file.path(directory, paste0("text ' $literal; [", i, "] `literal`.txt"))
    writeBin(charToRaw(contents[i]), path)
    compressed_path <- paste0(path, ".gz")
    con <- gzfile(compressed_path, "wb")
    writeBin(charToRaw(contents[i]), con)
    close(con)
    expect_identical(count_file_lines(path), counts[i])
    expect_identical(count_file_lines(compressed_path), counts[i])
  }
}

test_that("text line counts include blanks and unterminated final lines", {
  skip_if_not(all(nzchar(Sys.which(c("rg", "gzip")))))
  directory <- tempfile("text counts ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  check_text_line_counts(directory)
})

test_that("text line counts work with the grep fallback", {
  skip_on_os("windows")
  commands <- Sys.which(c("grep", "gunzip", "gzip"))
  skip_if_not(all(nzchar(commands)))
  directory <- tempfile("text counts fallback ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  bin <- file.path(directory, "bin")
  dir.create(bin)
  linked <- file.symlink(unname(commands), file.path(bin, names(commands)))
  skip_if_not(all(linked))
  original_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = original_path), add = TRUE)
  Sys.setenv(PATH = bin)
  expect_identical(unname(Sys.which("rg")), "")
  check_text_line_counts(directory)
})

test_that("line counting rejects invalid paths", {
  for (path in list(NULL, character(), NA_character_, "", c("a", "b"), 1)) {
    expect_error(count_file_lines(path), "single non-empty file path")
  }
  expect_error(count_file_lines(tempfile()), "existing file")
  expect_error(count_file_lines(tempdir()), "existing file")
})
