#' Read paired reads from RNA sequencing into memory
#'
#' Reads fastq files from Illumina paired sequencing and applied optional filtering
#' by quality metrics.
#'
#' @param fastq_path_r1 path to .fastq or .fastq.gz file of read1
#' @param fastq_path_r2 path to .fastq or .fastq.gz file of read2
#' @param min_len minimal length of reads to return
#' @param filter_ACTG filter for reads with letter ACTG only; remove all reads which contain N for instance
#' @param filter_paired if one of paired reads does not pass a test (min_len, filter_ACTG, quality_filter, max_freq_below_Q30),
#' then remove both reads (TRUE) or only the respective read (FALSE) by setting it to NA
#' @param quality_filter perform quality filtering on reads according to max_freq_below_Q30
#' @param max_freq_below_Q30 maximum frequency of bases per reads that are allowed to be below Q30
#' @param vroom_lines_args arguments to vroom::lines
#' @param mc.cores number of cores to use to parallel computing, speeds up quality filtering;
#' is limited to parallel::detectCores()
#' @param sample return a subset (sample) of paired reads only; can be a value greater
#' than 1 to sample a fixed absolute number (limited to total number of reads) or
#' a value lower than 1 to sample a fraction of reads
#'
#' @return a list with read sequences, their quality strings and statistics
#' @export
#'
#' @examples
#' #' \dontrun{
#' reads_paired <- igsc::read_paired_reads(fastq_path_r1 = r1_path, fastq_path_r2 = r2_path, vroom_lines_args = list(progress = T, skip_empty_rows = T), min_len = 40, mc.cores = 8)
#' }
read_paired_reads <- function(fastq_path_r1,
                              fastq_path_r2,
                              min_len = NULL,
                              filter_ACTG = T,
                              filter_paired = F,
                              quality_filter = T,
                              max_freq_below_Q30 = 0.1,
                              vroom_lines_args = list(progress = F, skip_empty_rows = T),
                              sample = 1,
                              mc.cores = 1) {

  # do not return qual_num - would also allow to make parallel computing easier

  if (any(!grepl("\\.fastq$", fastq_path_r1) || !grepl("\\.fastq$", fastq_path_r2)) && any(!grepl("\\.fastq\\.gz$", fastq_path_r1) || !grepl("\\.fastq\\.gz$", fastq_path_r2))) {
    stop("fastq_path_r1 and fastq_path_r2 have to be paths to fastq or fastq.gz files.")
  }

  if (max_freq_below_Q30 > 1 && quality_filter) {
    message("max_freq_below_Q30 is greater than 1. Please provide a fraction, e.g. 0.1 for 10 %. max_freq_below_Q30 is now set to 1.")
    max_freq_below_Q30 <- 1
  }

  if (mc.cores > 1) {
    mc.cores <- min(mc.cores, parallel::detectCores())
  }

  # filter_ACTG self-explained
  # filter_paired; remove r1 and r2 if one of them is removed due to filter_ACTG or quality_filter; if FALSE r1 or r2 is set NA and the other is kept
  message("reading reads into memory.")
  reads <- lapply(list(r1 = fastq_path_r1, r2 = fastq_path_r2), function(path) {
    if (grepl("\\.fastq\\.gz$", path)) {
      vroom_lines_args <- c(list(file = gzfile(path)), vroom_lines_args)
    } else {
      vroom_lines_args <- c(list(file = path), vroom_lines_args)
    }
    all_lines <- Gmisc::fastDoCall(vroom::vroom_lines, args = vroom_lines_args)
    if (quality_filter) {
      # make numeric qual below
      return(list(seq = all_lines[seq(2, length(all_lines), 4)],
                  qual = all_lines[seq(4, length(all_lines), 4)],
                  qual_num = NULL))
    } else {
      return(list(seq = all_lines[seq(2, length(all_lines), 4)],
                  qual = NULL,
                  qual_num = NULL))
    }
  })
  stat_total_reads <- length(reads[["r1"]][["seq"]])

  message("  ", format(length(reads[["r1"]][["seq"]]), big.mark=","), " paired reads.")
  if (length(unique(lengths(reads))) != 1) {
    stop("Unequal numbers of read 1 and read 2.")
  }


  if (!is.null(min_len)) {
    message("filter for min_len.")
    len_list <- list(r1 = nchar(reads[["r1"]][["seq"]]),
                     r2 = nchar(reads[["r2"]][["seq"]]))
    if (filter_paired) {
      both_meet_min_len <- which(len_list[["r1"]] >= min_len & len_list[["r2"]] >= min_len)
      message("  removed ", format(length(reads[["r1"]][["seq"]])-length(both_meet_min_len), big.mark=","), " paired reads (", round((length(reads[["r1"]][["seq"]])-length(both_meet_min_len))/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      for (x in c("r1", "r2")) {
        reads[[x]][["seq"]] <- reads[[x]][["seq"]][both_meet_min_len]
        if (quality_filter) {
          reads[[x]][["qual"]] <- reads[[x]][["qual"]][both_meet_min_len]
        }
      }
    } else {
      for (x in c("r1", "r2")) {
        reads[[x]][["seq"]][which(len_list[[x]] < min_len)] <- NA
      }
      both_min_len_NA <- which(is.na(reads[["r1"]][["seq"]]) & is.na(reads[["r2"]][["seq"]]))
      message("  removed ", format(length(both_min_len_NA), big.mark=","), " paired reads (", round(length(both_min_len_NA)/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      if (length(both_min_len_NA) > 0) {
        for (x in c("r1", "r2")) {
          reads[[x]][["seq"]] <- reads[[x]][["seq"]][-both_min_len_NA]
          message("  made ", format(length(which(len_list[[x]] < min_len))-length(both_min_len_NA), big.mark=","), " ", x,  " reads NA.")
          if (quality_filter) {
            # also filter which(r1_len < min_len) ??
            reads[[x]][["qual"]] <- reads[[x]][["qual"]][-both_min_len_NA]
          }
        }
      }
    }
  }


  # filter for reads with ATGC only in read seq
  if (filter_ACTG) {
    message("filter for ACTG only.")
    # r1_inds and r2_inds is TRUE for reads that only contain ACTG but N

    if (mc.cores > 1) {
      # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks

      split_sizes <- list(r1 = ceiling(length(reads[["r1"]][["seq"]])/mc.cores),
                          r2 = ceiling(length(reads[["r2"]][["seq"]])/mc.cores))
      inds <- purrr::map(stats::setNames(c("r1", "r2"), c("r1", "r2")), function(x) {
        !unlist(parallel::mclapply(split(reads[[x]][["seq"]], ceiling(seq_along(reads[[x]][["seq"]])/split_sizes[[x]])),
                                   grepl,
                                   pattern = "[^ACTG]",
                                   mc.cores = mc.cores))
      })
    } else {
      inds <- list(r1 = !grepl("[^ACTG]", reads[["r1"]][["seq"]]),
                   r2 = !grepl("[^ACTG]", reads[["r2"]][["seq"]]))
    }
    #inds_diff2 <- setdiff(which(r1_inds), which(r2_inds))
    #inds_diff1 <- setdiff(which(r2_inds), which(r1_inds))

    if (filter_paired) {
      inds_paired <- intersect(which(inds[["r1"]]), which(inds[["r2"]]))
      message("  removed ", format(length(reads[["r1"]][["seq"]])-length(inds_paired), big.mark=","), " paired reads (", round((length(reads[["r1"]][["seq"]])-length(inds_paired))/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      for (x in c("r1", "r2")) {
        reads[[x]][["seq"]] <- reads[[x]][["seq"]][inds_paired]
        if (quality_filter) {
          reads[[x]][["qual"]] <- reads[[x]][["qual"]][inds_paired]
        }
      }
    } else {
      reads[["r1"]][["seq"]][which(!inds[["r1"]])] <- NA
      reads[["r2"]][["seq"]][which(!inds[["r2"]])] <- NA
      both_NA_ACTG <- which(is.na(reads[["r1"]][["seq"]]) & is.na(reads[["r2"]][["seq"]]))
      message("  removed ", format(length(both_NA_ACTG), big.mark=","), " paired reads (", round(length(both_NA_ACTG)/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      if (length(both_NA_ACTG) > 0) {
        for (x in c("r1", "r2")) {
          reads[[x]][["seq"]] <- reads[[x]][["seq"]][-both_NA_ACTG]
          message("  made ", format(length(which(inds[[x]]))-length(both_NA_ACTG), big.mark=","), " ", x, " reads NA.")
          if (quality_filter) {
            reads[[x]][["qual"]] <- reads[[x]][["qual"]][-both_NA_ACTG]
          }
        }
      }
    }
  }

  if (quality_filter) {
    message("filter for quality.")
    # when mc.cores > 1 qual_num is not returned
    if (mc.cores > 1) {
      temp_fun <- function(chunks) {
        qual_num <- methods::as(Biostrings::PhredQuality(x = chunks), "IntegerList")
        qual_b30 <- unlist(lapply(qual_num, function(x) sum(x < 30)))
        b30_freq <- qual_b30/lengths(qual_num)
        return(b30_freq)
      }
      split_sizes <- list(r1 = ceiling(length(reads[["r1"]][["qual"]])/mc.cores),
                          r2 = ceiling(length(reads[["r2"]][["qual"]])/mc.cores))
      b30_freq <- purrr::map(stats::setNames(c("r1", "r2"), c("r1", "r2")), function(x) {
        unlist(parallel::mclapply(split(reads[[x]][["qual"]], ceiling(seq_along(reads[[x]][["qual"]])/split_sizes[[x]])),
                                  temp_fun,
                                  mc.cores = mc.cores))
      })
    } else {
      b30_freq <- list(r1 = "", r2 = "")
      for (x in c("r1", "r2")) {
        reads[[x]][["qual_num"]] <- methods::as(Biostrings::PhredQuality(reads[[x]][["qual"]]), "IntegerList")
        qual_b30 <- unlist(lapply(reads[[x]][["qual_num"]], function(x) sum(x < 30)))
        b30_freq[[x]] <- qual_b30/lengths(reads[[x]][["qual_num"]])
      }
    }

    if (filter_paired) {
      both_qual_filter <- which(b30_freq[["r1"]] <= max_freq_below_Q30 & b30_freq[["r2"]] <= max_freq_below_Q30)
      message("  removed ", format(length(reads[["r1"]][["seq"]])-length(both_qual_filter), big.mark=","), " paired reads (", round((length(reads[["r1"]][["seq"]])-length(both_qual_filter))/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      for (x in c("r1", "r2")) {
        reads[[x]][["seq"]] <- reads[[x]][["seq"]][both_qual_filter]
        reads[[x]][["qual"]] <- reads[[x]][["qual"]][both_qual_filter]
        reads[[x]][["qual_num"]] <- reads[[x]][["qual_num"]][both_qual_filter]
      }
    } else {
      reads[["r1"]][["seq"]][which(b30_freq[["r1"]] > max_freq_below_Q30)] <- NA
      reads[["r2"]][["seq"]][which(b30_freq[["r2"]] > max_freq_below_Q30)] <- NA
      both_qual_filter_NA <- which(is.na(reads[["r1"]][["seq"]]) & is.na(reads[["r2"]][["seq"]]))
      message("  removed ", format(length(both_qual_filter_NA), big.mark=","), " paired reads (", round(length(both_qual_filter_NA)/length(reads[["r1"]][["seq"]])*100, 0), " %).")
      if (length(both_qual_filter_NA) > 0) {
        for (x in c("r1", "r2")) {
          reads[[x]][["seq"]] <- reads[[x]][["seq"]][-both_qual_filter_NA]
          message("  made ", format(length(which(b30_freq[[x]] > max_freq_below_Q30))-length(both_qual_filter_NA), big.mark=","), " ", x, " reads NA.")
          reads[[x]][["qual"]] <- reads[[x]][["qual"]][-both_qual_filter_NA]
          reads[[x]][["qual_num"]] <- reads[[x]][["qual_num"]][-both_qual_filter_NA]
        }
      }
    }
  }

  if (sample != 1) {
    reads <- sample_reads(reads = reads, p = sample)
  }

  message("returning:")
  if (filter_paired) {
    message("  ", format(length(which(!is.na(reads[["r1"]][["seq"]]))), big.mark=","), " paired reads.")
  } else {
    for (x in c("r1", "r2")) {
      message("  ", format(length(which(!is.na(reads[[x]][["seq"]]))), big.mark=","), " reads ", x,".")
      message("  ", format(length(which(is.na(reads[[x]][["seq"]]))), big.mark=","), " reads ", x, " are NA.")
    }
  }


  stat_df <- data.frame(reads_total = stat_total_reads,
                        r1_post_qc_total = length(which(!is.na(reads[["r1"]][["seq"]]))),
                        r2_post_qc_total = length(which(!is.na(reads[["r2"]][["seq"]]))),
                        r1_post_qc_total_na = length(which(is.na(reads[["r1"]][["seq"]]))),
                        r2_post_qc_total_na = length(which(is.na(reads[["r2"]][["seq"]]))))

  return(list(reads = reads,
              stats = stat_df))
}

sample_reads <- function(reads, p = 1) {
  if (p != 1) {
    len <- length(reads_paired[["reads"]][["r1"]][["seq"]])
    if (p > 1) {
      p <- min(p, len)
    } else if (p < 1) {
      p <- round(p*len,0)
    }
    q <- sample(1:len, p)
    for (i in c("r1", "r2")) {
      for (j in c("seq", "qual")) {
        reads_paired[["reads"]][[i]][[j]] <- reads_paired[["reads"]][[i]][[j]][q]
      }
    }
  }
  return(reads)
}


