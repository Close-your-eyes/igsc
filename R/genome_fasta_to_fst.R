#' Title
#'
#' @param genome_file
#' @param subfolder "genome_fst"
#' @param overwrite
#' @param compression
#' @param verbose
#' @param simplify_names
#'
#' @return
#' @export
#'
#' @examples
genome_fasta_to_fst <- function(genome_file,
                                subfolder = "genome_fagz",
                                # fst for random access, but large as a whole data frame read into mem
                                overwrite = F,
                                compression = 100,
                                verbose = T,
                                simplify_names = T) {
  if (missing(genome_file) || length(genome_file) == 0) {
    stop("path to genome_file missing.")
  }
  if (!file.exists(genome_file)) {
    stop("genome_file not found.")
  }
  subfolder <- match.arg(subfolder, c("genome_fst", "genome_fagz"), several.ok = T)

  seq_bounds <- get_fasta_seq_bounds(file_path = genome_file)

  for (y in seq_along(subfolder)) {
    subfolder[y] <- file.path(dirname(genome_file), subfolder[y])
    dir.create(subfolder[y], showWarnings = F)
  }

  for (x in seq_along(seq_bounds$seqname)) {
    if (simplify_names) {
      name <- strsplit(seq_bounds$seqname[x], " ")[[1]][1]
    } else {
      name <- seq_bounds$seqname[x]
    }

    for (y in seq_along(subfolder)) {
      out_path <- file.path(subfolder[y], paste0(name, "_fa.", ifelse(grepl("fst$", subfolder[y]), "fst", "fa.gz")))
      if (verbose && y == 1) {
        message(name, " ", x, "/", length(seq_bounds$seqname), " ≈", format(seq_bounds$bp_approx[x], big.mark = ","), " bp")
      }
      skip <- file.exists(out_path)
      if (!skip && overwrite) {
        skip <- F
      }

      if (!skip) {
        if (y == 1) {
          refseq <- read_fasta(genome_file,
                               start_line = seq_bounds[x, "start_line",drop=T],
                               end_line = seq_bounds[x, "end_line",drop=T])

          # complement seq_bounds file
          seq_bounds[which(seq_bounds$seqname == name), "bp"] <- nchar(refseq)
          seq_bounds[which(seq_bounds$seqname == name), "lead_N"] <- nchar(refseq) - nchar(gsub("^N+", "", refseq))
          seq_bounds[which(seq_bounds$seqname == name), "trail_N"] <- nchar(refseq) - nchar(gsub("N+$", "", refseq))
        }

        if (grepl("fst$", subfolder[y])) {
          refseqdf <- data.frame(x = strsplit(x = refseq, split = "", fixed = T)[[1]])
          names(refseqdf) <- name
          fst::write_fst(refseqdf, path = out_path, compress = compression)
        } else {
          con <- file(description = gsub("\\.gz$", "" , out_path), open = "w")
          writeLines(paste0(">", seq_bounds[x, "seqname",drop=T]), con)
          writeLines(refseq, con)
          close(con)
          system(paste0("gzip -1 ", gsub("\\.gz$", "" , out_path)))
        }
      }
    }
  }
  if ("file_path" %in% names(attributes(seq_bounds))) {
    fst::write_fst(seq_bounds, path = attributes(seq_bounds)[["file_path"]])
  }
  if (tools::file_ext(genome_file) != "gz") {
    system(paste0("gzip -1 ", genome_file)) # not much cpu but large disk saving
  }
}

# binary is as quick as rds but w/o compression

# saveRDS(refseq, "/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/refseq.rds", compress = F)
# refseq <- readRDS("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/refseq.rds")
#
# con <- file("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.bin", "wb")
# writeBin(charToRaw(refseq), con)
# close(con)
#
# con <- file("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.bin", "rb")
# dna_sequence <- rawToChar(readBin(con, raw(), file.info("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.bin")$size))
# close(con)
#
# rle18 <- rle(strsplit(refseq, "")[[1]])
# rest <- paste(rep(rle18$values, rle18$lengths), collapse = "")
# rest == refseq
# ref2 <- fst::read_fst("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/genome_fa_fst/chr18.fst")
# gz_connection <- file(description = "/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.fa", open = "w")
# #vroom::vroom_write_lines(x = refseq, file = gz_connection)
# writeLines(">chr18", gz_connection)
# writeLines(refseq, gz_connection)
# close(gz_connection)
# system(paste0("gzip -3 /Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.fa"))
#
#
# fff <- read_fasta("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.fa.gz")
#
# vroom::vroom_write_lines(refseq, "/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.fa")
# gzfile("/Users/vonskopnik/Documents/2024_igsc_testing/refdata-gex-GRCh38-2020-A/chr18.fa")
# this is to slow/large, random access cannot be handled easily by vroom lines
#writeLines(paste(strsplit(refseq, "")[[1]], sep = "\n"), con)
