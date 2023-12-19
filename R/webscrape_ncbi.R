#' Scrape data from NCBI nucleotide pages
#'
#' @param accession accession number to get data from
#'
#' @return
#' @export
#'
#' @examples
webscrape_ncbi <- function(accession) {
  ncbi_text <- rentrez::entrez_fetch("nucleotide", accession, rettype = "text")

  # check ncbi_text ?

  if (!grepl("ORIGIN", ncbi_text, ignore.case = F)) {
    message("ORIGIN not found in text from NCBI.")
    origin <- NULL
  } else {
    origin <- get_origin(x = ncbi_text)
    names(origin) <- accession
  }

  if (!grepl("FEATURES", ncbi_text, ignore.case = F)) {
    message("FEATURES not found in text from NCBI.")
    features <- NULL
  } else {
    if (length(attr(regexec("FEATURES", ncbi_text, fixed = T)[[1]], "match.length")) > 1) {
      message("FEATURES found more than once. Check results.")
    }
    features <- get_features(x = ncbi_text)
  }
  header <- get_header(x = ncbi_text)

  return(list(header = header, features = features, origin = origin))
}

get_header <- function(x) {
  LOCUS <- trimws(strsplit(substr(x, 1, 2000), "DEFINITION")[[1]][1])
  DEFINITION <- trimws(strsplit(strsplit(substr(x, 1, 2000), "DEFINITION")[[1]][2], "\n")[[1]][1])
  ACCESSION <- trimws(strsplit(strsplit(substr(x, 1, 2000), "ACCESSION")[[1]][2], "\n")[[1]][1])
  VERSION <- trimws(strsplit(strsplit(substr(x, 1, 2000), "VERSION")[[1]][2], "\n")[[1]][1])
  return(list(locus = LOCUS,
              definition = DEFINITION,
              accession = ACCESSION,
              version = VERSION))
}

get_origin <- function(x) {

  seq <- strsplit(x, "ORIGIN", fixed = T)[[1]]
  seq <- seq[length(seq)] # get last index in any case
  seq <- gsub(" {1,}", "", gsub("\n", "", seq))
  seq <- toupper(gsub("//", "", gsub("[[:digit:]]{1,}", "", seq)))
  #names(seq) <- "origin"
  return(seq)
}


get_features <- function(x) {

  seq <- strsplit(x, "FEATURES", fixed = T)[[1]][2]
  seq <- strsplit(seq, "ORIGIN", fixed = T)[[1]]
  seq <- seq[-length(seq)] # remove last index
  seqs <- strsplit(seq, "\n")[[1]][-1]
  seqs <- fix_lines2(z = seqs)
  seqs <- trimws(gsub(" {1,}", " ", seqs))
  inds_log <- !grepl("^/", seqs) # redetect feature lines

  inds <- which(inds_log)
  alt_inds <- which(!inds_log)

  meta_data_belonging <- cut(alt_inds, breaks = inds, labels = F, include.lowest = T)
  meta_data_belonging[which(is.na(meta_data_belonging))] <- max(meta_data_belonging, na.rm = T) + 1 # the alt_inds belonging to final feature are returned as NA from cut-function; label them manually

  list_grouped <- setNames(split(seqs[alt_inds], meta_data_belonging), seqs[inds])
  df <- prep_list(x = list_grouped)

  df_distinct <- dplyr::distinct(df, ID, Feature)
  message("Features found:")
  print(table(df_distinct$Feature))

  return(df)
}

fix_translation_lines <- function(z) {
  # treat lines following translation differently - first - add all lines to the starting line of translation
  translation_lines <- which(grepl("^ {1,}/translation=", z))
  translation_ranges <- lapply(setNames(translation_lines, translation_lines), function(y) {
    start_gap <- nchar(z[y]) - nchar(trimws(z[y]))
    while (start_gap == 21) {
      y <- y + 1
      start_gap <- nchar(z[y]) - nchar(trimws(z[y]))
    }
    return(y-1)
  })

  translation_ranges <- seq2(as.numeric(names(translation_ranges)), unlist(unname(translation_ranges)))

  for (i in translation_ranges) {
    z[i[1]] <- paste0(paste(rep(" ", 21), collapse = ""), gsub(" ", "", paste(z[i], collapse = "")))
  }
  translation_ranges <- lapply(translation_ranges, function(x) x[-1])
  z <- z[-unlist(translation_ranges)]
  return(z)
}

fix_lines2 <- function(z, exclude_lines = NULL) {
  # fix line breaks and paste them to the first line, each
  break_lines <- which(!grepl("^ {1,}/", z))
  break_lines <- break_lines[!grepl("^ {5}[[:alpha:]]", z[break_lines])] # find exactly 5 spaces at the beginning


  # give meaningful names to different funs, which find_consecutive_sequences
  runs <- find_consecutive_sequences(x = break_lines)
  lines_to_concat <- lapply(runs, function(x) c(min(x)-1, x))
  # find_consecutive_sequences2(lines_to_concat)

  for (i in lines_to_concat) { #consecutive_sequences
    z[i[1]] <- paste(c(z[i[1]], trimws(z[i[-1]])), collapse = "")
  }
  z <- z[-break_lines]
  return(z)
}

fix_lines <- function(z, exclude_lines = NULL) {
  # fix line breaks and paste them to the first line, each
  lines <- which(!grepl("^/", trimws(z)))
  lines <- setdiff(lines, exclude_lines) # exclude feature lines

  ## maybe make a separate function?
  # Extract the consecutive sequences
  runs <- rle(diff(lines) == 1)
  consecutive_positions <- which(runs$values)
  consecutive_sequences <- Map(function(start, end) lines[start:(end + 1)],
                               consecutive_positions,
                               consecutive_positions + runs$lengths[consecutive_positions] - 1)


  for (i in consecutive_sequences) {
    z[i[1]-1] <- paste(c(z[i[1]-1], z[i]), collapse = " ")
  }
  z <- z[-lines]
  return(z)
}

prep_list <- function(x) {

  x_split <- strsplit(names(x), " ")
  feature <- sapply(x_split, "[", 1)
  range <- sapply(x_split, "[", 2)
  out <- purrr::map_dfr(x, function(y) {
    z <- strsplit(y, "=")
    stack(setNames(gsub("\"", "", unlist(sapply(z, "[", 2))), gsub("^/", "", unlist(sapply(z, "[", 1)))))
  }, .id = "ID")
  out$ind <- as.character(out$ind)
  out$Feature <- rep(feature, lengths(x))
  out$range <- rep(range, lengths(x))


  if (any(grepl("order", out$range))) {
    message("order found in range of features Not sure how to handle that currently. Are the Ns in origin?")
  }
  if (any(grepl(">", out$range)) || any(grepl("<", out$range))) {
    message("'<' or '>' found in range of features This may not handable by other methods.")
  }

  out$range <- gsub("order\\(", "", out$range)
  out$ind <- as.character(out$ind)
  out$complement <- grepl("complement", out$range)
  out$range <- gsub("complement\\(", "", out$range)
  out$range <- gsub("\\)", "", out$range)
  out$range <- gsub("join\\(", "", out$range)
  names(out)[which(names(out) == "ind")] <- "Subfeature"
  names(out)[which(names(out) == "values")] <- "value"
  out <- out[,c(1,4,3,5,6,2)]


  # this puts range as a separate row, not column
  ' out2 <- data.frame(ind = "range", Feature = feature, values = range, ID = names(x))
  out2$ind <- ifelse(grepl("complement", out2$values), "range_revcomp", "range")
  out2$values <- gsub("complement\\(", "", out2$values)
  out2$values <- gsub("\\)", "", out2$values)
  out2$values <- gsub("join\\(", "", out2$values)
  out <- rbind(out, out2)
  names(out)[which(names(out) == "ind")] <- "Subfeature"
  names(out)[which(names(out) == "values")] <- "value"
  out <- out[,c(1,4,3,2)]
  out <- out[order(match(out$ID, names(x))),]
  rownames(out) <- as.character(seq(1,nrow(out)))'

  return(tibble::as_tibble(out))
}

find_consecutive_sequences <- function(x) {
  diffs <- c(1, diff(x))
  starts <- c(1, which(diffs != 1))
  ends <- c(starts[-1] - 1, length(x))

  consecutive_sequences <- Map(function(start, end) x[start:end], starts, ends)

  return(consecutive_sequences)
}

find_consecutive_sequences2 <- function(x) {
  runs <- rle(diff(x) == 1)
  consecutive_positions <- which(runs$values)
  consecutive_sequences <- Map(function(start, end) x[start:(end + 1)],
                               consecutive_positions,
                               consecutive_positions + runs$lengths[consecutive_positions] - 1)

}
