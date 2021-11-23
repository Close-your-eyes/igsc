#' Prepare a data frame with TCR gene segment reference data from IMGT
#'
#' Immunoglobulin (IG) reference data from IMGT do not come in a handy format for processing in R.
#' For T cell receptor (TCR) gene segments, this functions uses data from IMGT (fasta files and one manually prepared xlsx) to
#' create a data frame that can be used subsequently to align TCR sequences from scRNAseq (or other).
#' All necessary files (human or mouse) are included in this package (Oct-2021) but may be downloaded manually from IMGT in case there are major updates.
#' The files included can be retrieved with file.copy(list.files(system.file("extdata", "IMGT_ref/human", package = "igsc")), 'path to your folder') or
#' file.copy(list.files(system.file("extdata", "IMGT_ref/mouse", package = "igsc")), 'path to your folder').
#' These files may also demonstrate the required file names and structures in case you want to provide updated data from IMGT.
#' To skip this function and immediatly obtain its output, ready made data frames are available with
#' imgt_ref <- readRDS(system.file("extdata", "IMGT_ref/human/hs.rds", package = "igsc")) or
#' imgt_ref <- readRDS(system.file("extdata", "IMGT_ref/mouse/mm.rds", package = "igsc")).
#'
#' Sources and how to prepare the data yourself. Data for the xlsx-files are from:
#' http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=human&latin=Homo%20sapiens&group=TRAV,
#' http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=human&latin=Homo%20sapiens&group=TRBV,
#' http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRAV,
#' http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRBV.
#' Fasta-files are made from the data found here:
#' http://www.imgt.org/vquest/refseqh.html.
#' Leader sequences are from "L-PART1+L-PART2" artificially spliced sets, nucleotides (F+ORF+all P).
#' Others are from "L-PART1+V-EXON" artificially spliced sets and Constant gene artificially spliced exons sets.
#' Fasta-formatted sequences from there have to be copied manually and saved as .fasta files in a folder. This folder then becomes the path argument.
#'
#' @param path path to a folder with all necessary files from IMGT; if not provided human or mouse data downloaded roughly Oct-2021 will be used
#' @param organism if no path is provided data will be taken from this package, either human or mouse
#' @param mc use multicore (mclapply from parallel package) for pairwise alignment of TCR segments
#'
#' @return a data frame
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' imgt_df <- imgt_tcr_segment_prep()
#' openxlsx::write.xlsx(imgt_df, "imgt_ref_df.xlsx")
#' saveRDS(imgt_df, "imgt_ref_df.rds")
#' }
imgt_tcr_segment_prep <- function(path, organism = "human", mc = F) {

  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"Biostrings" %in% rownames(utils::installed.packages())) {BiocManager::install("Biostrings")}

  if (missing(path)) {
    if (organism == "human") {
      path <- system.file("extdata", "IMGT_ref/human", package = "igsc")
    }
    if (organism == "mouse") {
      path <- system.file("extdata", "IMGT_ref/mouse", package = "igsc")
    }
  } else {
    if (!is.character(path) || !file.exists(path)) {
      stop("path not found or not a character")
    }
  }
  organism <- match.arg(organism, c("human", "mouse"))
  if (!is.logical(mc)) {
    stop("mc has to be T or F.")
  }

  files <- list.files(path, "\\.fasta", recursive = T, full.names = T)
  files <- files[which(!grepl("leader", basename(files)))]
  ts <- dplyr::bind_rows(lapply(files, function(x) {
    utils::stack(read_fasta(x))
  })) %>%
    dplyr::rename("seq.nt" = values, "meta" = ind)
  ts$Allele <- stringr::str_replace(sapply(stringr::str_split(ts$meta, "\\|"), "[", 2),"/", "")
  ts$Gene <- stringr::str_replace(sapply(stringr::str_split(ts$Allele, "\\*"), "[", 1),"/", "")
  ts$Allele.number <- sapply(stringr::str_split(ts$Allele, "\\*"), "[", 2)
  ts$AccNum <- stringr::str_replace(sapply(stringr::str_split(ts$meta, "\\|"), "[", 1),"/", "")
  ts$Functionality <- stringr::str_replace(sapply(stringr::str_split(ts$meta, "\\|"), "[", 4),"/", "")

  ## all functional V segments start with an ATG
  ts$seq.aa <- unlist(lapply(split(ts, seq(nrow(ts))), function(x) {
    if (!grepl("V", x[,"Gene"]) || grepl("partial.in.5", x[,"meta"]) || x[,"Functionality"] != "F") {
      NA
    } else {
      as.character(Biostrings::translate(Biostrings::DNAStringSet(x[,"seq.nt"])))
    }
  }))

  leader.seq <- utils::stack(read_fasta(list.files(path, "TRV_leader_aa.fasta", recursive = T, full.names = T))) %>%
    dplyr::rename("LEADER" = values, "meta" = ind) %>%
    dplyr::distinct()
  leader.seq$Allele <- stringr::str_replace(sapply(stringr::str_split(leader.seq$meta, "\\|"), "[", 2), "/", "")
  leader.seq$Gene <- stringr::str_replace(sapply(stringr::str_split(leader.seq$Allele, "\\*"), "[", 1), "/", "")
  leader.seq$Allele.number <- sapply(stringr::str_split(leader.seq$Allele, "\\*"), "[", 2)
  leader.seq <- dplyr::select(leader.seq, -meta)

  imgt_cdr_fr <-
    dplyr::bind_rows(lapply(list.files(path, "\\.xlsx", recursive = T, full.names = T), function(x) {openxlsx::read.xlsx(x)})) %>%
    dplyr::mutate(FR1 = stringr::str_replace_all(FR1, "\\.", "")) %>%
    dplyr::mutate(FR2 = stringr::str_replace_all(FR2, "\\.", "")) %>%
    dplyr::mutate(FR3 = stringr::str_replace_all(FR3, "\\.", "")) %>%
    dplyr::mutate(CDR1 = stringr::str_replace_all(CDR1, "\\.", "")) %>%
    dplyr::mutate(CDR2 = stringr::str_replace_all(CDR2, "\\.", "")) %>%
    dplyr::mutate(Allele = stringr::str_replace_all(Allele, "/", "")) %>%
    dplyr::select(-c(Species, Allele, AccNum, Functionality, Domain.label))

  ts <-
    ts %>%
    dplyr::left_join(leader.seq, by = c("Allele", "Gene", "Allele.number")) %>%
    dplyr::left_join(imgt_cdr_fr, by = "Gene")


  al_fun <- function(x) {
    if (is.na(x[,FR]) || x[,"Functionality"] != "F") {
      NA
    } else {
      Biostrings::pairwiseAlignment(subject = x[,"seq.aa"], pattern = x[,FR], type = "local")
    }
  }

  for (FR in c("LEADER", "FR1", "CDR1", "FR2", "CDR2", "FR3")) {
    if (mc) {
      out <- parallel::mclapply(split(ts, seq(nrow(ts))), al_fun, mc.cores = parallel::detectCores())
    } else {
      out <- pbapply::pblapply(split(ts, seq(nrow(ts))), al_fun)
    }
    ts[,paste0(FR, ".start.aa")] <- sapply(out, function(x) ifelse(is.na(x), NA, x@subject@range@start))
    ts[,paste0(FR, ".end.aa")] <- sapply(out, function(x) ifelse(is.na(x), NA, x@subject@range@start + x@subject@range@width - 1))
    ## formula only works since the first ATG is at position 1
    ts[,paste0(FR, ".start.nt")] <- (ts[,paste0(FR, ".start.aa")] - 1)*3 + 1
    ts[,paste0(FR, ".end.nt")] <- ts[,paste0(FR, ".end.aa")]*3
  }
  ts[,"CDR3.start.aa"] <- ts[,"FR3.end.aa"] + 1
  ts[,"CDR3.start.nt"] <- ts[,"FR3.end.nt"] + 1

  return(ts)
}


read_fasta <- function(file, legacy.mode = T, seqonly = F) {
  # taken from the protr package - needs one modification
  lines <- readLines(file)
  if (legacy.mode) {
    comments <- grep("^;", lines)
    if (length(comments) > 0) {
      lines <- lines[-comments]
    }
  }
  ind <- which(substr(lines, 1L, 1L) == ">")
  nseq <- length(ind)
  if (nseq == 0)
    stop("no line starting with a > character found")
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]],
                                                       collapse = ""))
  if (seqonly)
    return(sequences)
  nomseq <- lapply(seq_len(nseq), function(i) {
    # changed here
    #firstword <- stringr::str_split(lines[ind[i]], " ")[[1]][1]
    #substr(firstword, 2, nchar(firstword))
    substr(trimws(lines[ind[i]]), 2, nchar(trimws(lines[ind[i]])))
  })
  names(sequences) <- nomseq
  sequences
}



