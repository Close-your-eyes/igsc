tt <- vroom::vroom_lines("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq0/bamtofastq_S1_L001_R1_001.fastq.gz",
                         n_max = 40)

tt1 <- vroom::vroom_lines("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq0/bamtofastq_S1_L001_R2_001.fastq.gz",
                         n_max = 40)


line1 <- tt[seq(1,length(tt),4)]

concat_10X_fastq("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB",
                 exec = T, method = "binary")

n <- 4
cmd <- sprintf(
  "gzip -dc %s | awk 'NR %% %d == 0'",
  shQuote("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq/bamtofastq_S1_L001_R1_001.fastq.gz"), n
)

con <- pipe(cmd)
out <- tryCatch(
  vroom::vroom_lines(con),
  finally = close(con)
)

get_fasta_seq_bounds()

zz <- vroom::vroom_lines("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq/bamtofastq_S1_L001_R1_001.fastq.gz",
                         n_max = 400)
zz[399]


file_path <- "/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq/bamtofastq_S1_L001_R1_001.fastq.gz"
cmd <- paste0("rg -z --count '^' ", file_path)
ll <- as.numeric(trimws(system(cmd, intern = T)))
ll/4

make_ids <- function(
    n,
    alphabet = c(letters, as.character(0:9))) {

  # stringr::str_pad() ?

  stopifnot(
    is.numeric(n), length(n) == 1L, is.finite(n),
    n >= 0, n == floor(n),
    is.character(alphabet), length(alphabet) >= 2L,
    !anyNA(alphabet), !anyDuplicated(alphabet),
    all(nchar(alphabet) == 1L)
  )
  if (n == 0) return(character(0))

  width <- 1L
  capacity <- length(alphabet)

  while (capacity < n) {
    width <- width + 1L
    capacity <- capacity * length(alphabet)
  }

  grid <- expand.grid(
    rep(list(alphabet), width),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Reverse columns to match your original ID order.
  grid <- grid[seq_len(n), rev(seq_len(width)), drop = FALSE]

  do.call(paste0, grid)
}


ids <- make_ids(1e5, alphabet = as.character(0:9))
tail(ids)
lobstr::obj_size(ids)

ids <- make_ids(ll/4)

#file_path <- "/Users/chris/chunks/bamtofastq_S1_L001_R1_001.part_001.fastq.gz"
files <- list.files("/Users/chris/chunks", full.names = T)
zz <- vroom::vroom_lines("/Users/chris/chunks/bamtofastq_S1_L001_R1_001.part_001.fastq.gz")
zz

# files: your pasted vector of FASTQ filenames
# ids: character vector of new read IDs, without "@"

r1 <- sort(files[grepl("_R1_", basename(files))])
r2 <- sub("_R1_", "_R2_", r1, fixed = TRUE)

stopifnot(
  length(r1) > 0L,
  all(r2 %in% files),
  all(file.exists(c(r1, r2))),
  is.character(ids),
  !anyNA(ids),
  all(nzchar(ids)),
  !any(grepl("[[:space:]]", ids)),
  !any(startsWith(ids, "@"))
)

out_dir <- "/Users/chris/renamed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (mate_files in list(r1, r2)) {
  offset <- 0

  for (file in mate_files) {
    message("Processing: ", basename(file))

    x <- vroom::vroom_lines(file)

    # Assumes standard four-line FASTQ records.
    if (length(x) == 0L || length(x) %% 4L != 0L) {
      stop("Empty file or invalid four-line FASTQ: ", file)
    }

    headers <- seq.int(1L, length(x), by = 4L)
    plus <- headers + 2L
    n_reads <- length(headers)

    stopifnot(
      all(startsWith(x[headers], "@")),
      all(startsWith(x[plus], "+"))
    )

    if (offset + n_reads > length(ids)) {
      stop("Not enough IDs for: ", file)
    }

    idx <- offset + seq_len(n_reads)
    x[headers] <- paste0("@", ids[idx])
    x[plus] <- "+"  # Remove any repeated old header on the "+" line.

    vroom::vroom_write_lines(
      x,
      file.path(out_dir, basename(file))
    )

    offset <- offset + n_reads
    rm(x)
  }

  if (offset != length(ids)) {
    stop("The number of reads does not match the number of IDs.")
  }
}


rr <- vroom::vroom_lines("/Users/chris/fastq_out/Pat4_rep1_blood_0_1_H3HNYBGXB/concatfastq/bamtofastq/bamtofastq_S1_L001_R2_001_renamed_constant.fastq.gz",
                         n_max = 40)
rr
