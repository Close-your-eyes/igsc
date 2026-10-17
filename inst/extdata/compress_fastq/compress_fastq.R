tt <- vroom::vroom_lines("/Volumes/CMS_SSD_2TB/R_scRNAseq/2019_SLE_LN/data/Sequencing_data/hg38_CR7_wo_introns/GEX/Pat4_rep1_blood/outs/mhc_reads_bamboozled.fastq.gz")


line1 <- tt[seq(1,length(tt)/10,4)]
nchar(line1)

namehash <- function(x, mc.cores) {
  hashed <- unlist(
    parallel::mclapply(
      enc2utf8(x),
      digest::digest,
      algo = "xxh3_128",
      serialize = FALSE,
      mc.cores = 4L,
      mc.set.seed = FALSE
    ),
    use.names = FALSE
  )
  return(hashed)
}

namehash(line1)
line1x <- namehash(tt[seq(1,length(tt),4)])
line1stream <- paste(line1x,collapse = "")


line2 <- tt[seq(2,length(tt)/10,4)]

make_ids <- function(n) {
  if (n == 0L) return(character(0))

  alphabet <- c(letters, LETTERS, as.character(0:9))
  i <- seq_len(n) - 1L
  ids <- rep("", n)

  repeat {
    ids <- paste0(alphabet[i %% 62L + 1L], ids)
    i <- i %/% 62L
    if (all(i == 0L)) break
  }

  ids
}

make_id_stream <- function(n) {
  ids <- make_ids(n)
  if (!length(ids)) return("")

  alphabet <- c(letters, LETTERS, as.character(0:9))
  width <- nchar(ids[1L])

  paste0(alphabet[width], paste0(ids, collapse = ""))
}

restore_ids <- function(stream) {
  if (!nzchar(stream)) return(character(0))

  alphabet <- c(letters, LETTERS, as.character(0:9))
  width <- match(substr(stream, 1L, 1L), alphabet)
  payload <- substring(stream, 2L)

  starts <- seq.int(1L, nchar(payload), by = width)
  substring(payload, starts, starts + width - 1L)
}


pack_dna <- function(dna) {
  codes <- match(strsplit(paste0(dna, collapse = ""), "", fixed = TRUE)[[1L]],
                 c("A", "T", "C", "G")) - 1L

  bits <- as.raw(rbind(codes %% 2L, codes %/% 2L))
  bits <- c(bits, raw((-length(bits)) %% 8L))

  packBits(bits, type = "raw")
}

unpack_dna <- function(data) {
  n_reads <- (length(data) * 4L) %/% 98L
  if (n_reads == 0L) return(character(0))

  bits <- matrix(as.integer(rawToBits(data)), nrow = 2L)
  codes <- colSums(bits * c(1L, 2L))

  dna <- paste0(
    c("A", "T", "C", "G")[codes[seq_len(n_reads * 98L)] + 1L],
    collapse = ""
  )

  starts <- seq.int(1L, n_reads * 98L, by = 98L)
  substring(dna, starts, starts + 97L)
}

pack_hashes <- function(hashes) {
  hex <- paste0(hashes, collapse = "")
  if (!nzchar(hex)) return(raw(0))

  starts <- seq.int(1L, nchar(hex), by = 2L)
  as.raw(strtoi(substring(hex, starts, starts + 1L), base = 16L))
}

unpack_hashes <- function(data) {
  if (!length(data)) return(character(0))

  hex <- paste0(format(data), collapse = "")
  starts <- seq.int(1L, nchar(hex), by = 32L)

  substring(hex, starts, starts + 31L)
}

pack_qualities <- function(qualities) {
  values <- as.integer(charToRaw(paste0(qualities, collapse = "")))
  if (!length(values)) return(raw(0))

  dictionary <- sort(unique(values))
  codes <- match(values, dictionary) - 1L
  b <- max(1L, ceiling(log2(length(dictionary))))

  bits <- matrix(raw(length(codes) * b), nrow = b)
  for (j in seq_len(b)) {
    bits[j, ] <- as.raw(bitwAnd(bitwShiftR(codes, j - 1L), 1L))
  }
  bits <- as.raw(bits)
  bits <- c(bits, raw((-length(bits)) %% 8L))

  c(
    as.raw(length(dictionary)),
    as.raw(dictionary),
    packBits(bits, type = "raw")
  )
}

unpack_qualities <- function(stream) {
  if (!length(stream)) return(character(0))

  k <- as.integer(stream[1L])
  dictionary <- stream[seq_len(k) + 1L]
  data <- stream[-seq_len(k + 1L)]
  b <- max(1L, ceiling(log2(k)))

  n_reads <- (length(data) * 8) %/% (98L * b)
  n_values <- n_reads * 98L
  if (!n_values) return(character(0))

  bits <- matrix(
    as.integer(rawToBits(data)[seq_len(n_values * b)]),
    nrow = b
  )
  codes <- colSums(bits * 2^(seq_len(b) - 1L))
  qualities <- rawToChar(dictionary[codes + 1L])

  starts <- seq.int(1L, n_values, by = 98L)
  substring(qualities, starts, starts + 97L)
}

# pack line1
con <- gzfile("line1.bin.gz", open = "wb", compression = 1)
writeBin(make_id_stream(length(tt)/4), con)
close(con)

# pack line2
con <- gzfile("line2.bin.gz", open = "wb", compression = 1)
writeBin(pack_dna(tt[seq(2,length(tt),4)]), con)
close(con)

source("/Users/chris/Documents/Codex/2026-10-17/h/outputs/quality_huffman.R")

# pack line4
con <- gzfile("line4.bin.gz", open = "wb", compression = 1)
writeBin(pack_qualities(tt[seq(4,length(tt),4)]), con)
close(con)



