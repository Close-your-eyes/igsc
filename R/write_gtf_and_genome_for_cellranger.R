#' Prepare files for creating a reference genome with mkref
#'
#' The function is specifically designed to receive the return values from
#' igsc::webscrape_ncbi. One list entry of data (see the fun argument) would be one returned
#' list from that function. Moreover it is designed to allow passing its results to cellrangers
#' mkref in order to create a reference genome for mapping of single cell RNA seq data.
#' See documentations here: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr,
#' https://www.ensembl.org/info/website/upload/gff.html, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#multi,
#' https://kb.10xgenomics.com/hc/en-us/articles/115003327112-How-can-we-add-genes-to-a-reference-package-for-Cell-Ranger-
#'
#'
#' @param data a named list; each list element has to contain a character (string) named "origin"
#' which becomes one entry in the genome fasta file; the respective name of the list element
#' become the header for the fasta entry; classically this would be one chromosome in a hg38 ref;
#' another element of each list entry is a data frame named "features" which contains the features
#' of origin; the data frame has to have these columns: "range", "complement", "Feature", "value";
#' checkout the return of igsc::webscrape_ncbi to understand the columns; entries in the value
#' column are converted by make.names to replace special symbols by dots
#' @param genome_file_name filename of the genome file; genome.fa by default
#' @param gtf_file_name filename of the gtf file; genes.gtf by default
#' @param save_path path where to save genome_file and gtf_file
#' @param gtf_header header lines for gtf files, each vector element becomes one line
#' @param features_to_become_exon features from the Feature column in features data frames
#' which become an entry called 'exon' in the gtf file; cellranger only considers entry which
#' have 'exon' in the respective column; so even when we do not write exon boundaries into the
#' gtf file but those of, e.g. a whole protein coding sequence, this must be named exon
#' @param other_features_to_write other features from Feature column in features data frames
#' to write into the gtf file
#' @param append append genome_file and gtf_file?; only one of append or overwrite can be TRUE
#' @param overwrite overwrite existing genome_file and gtf_file?; only one of append or overwrite can be TRUE
#' @param genome_file_linewidth max number of character per line in genome_file; 60 by default
#'
#' @return genome_file_name and gtf_file_name written to disk at save_path; their content returned as character vectors
#' @export
#'
#' @examples
#' \dontrun{
#' # example accession from viral genomes
#' # these are small genomes hence the complete genome sequence is found under these accession
#' # for larger genomes, e.g. from bacteria, checkout NCBIs genome repository: https://www.ncbi.nlm.nih.gov/datasets/genome/
#' viral_genome_accessions <- c("NC_007605.1","NC_006273.2","NC_001806.2","NC_001798.2","NC_001348.1")
#' ncbi_data_list <- lapply(stats::setNames(viral_genome_accessions, viral_genome_accessions), function(x) igsc::webscrape_ncbi(accession = x))
#' # select a subset of features only
#' # other modification to features data are possible to restrict what ends up in the gtf file
#' ncbi_data_list <- purrr::map(ncbi_data_list, function(x) {
#' x[["features"]] <- dplyr::filter(x[["features"]], Subfeature == "locus_tag")
#' return(x)
#' })
#' # this write files with viral sequences only
#' write_gtf_and_genome_for_cellranger(data = ncbi_data_list, save_path = file.path(wd, "HHV_ref_genome"))
#'
#' # to append hg38: download the reference genome from 10X genomics then run
#' # write_gtf_and_genome_for_cellranger, select genome.fa and genes.gtf from this ref genome
#' # as genome_file_name and gtf_file_name and set append = TRUE.
#' # it was noted that more UMIs are assigned to ref sequences from human herpesvirus
#' # when a pure HHV ref is used compared to when hg38 is appended
#' }
write_gtf_and_genome_for_cellranger <- function(data,
                                                genome_file_name = "genome.fa",
                                                gtf_file_name = "genes.gtf",
                                                save_path = NULL,
                                                gtf_header = c("##description: made with write_gtf_and_genome_for_cellranger function from https://github.com/Close-your-eyes/igsc", "##provider: CMS", "##conctact: vonskopnik@pm.me", "##format: gtf", paste0("##date: ", Sys.Date())),
                                                features_to_become_exon = c("CDS"),
                                                other_features_to_write = NULL,
                                                append = F,
                                                overwrite = F,
                                                genome_file_linewidth = 60) {


  if (!is.list(data) || is.null(names(data))) {
    stop("data has to be a named list.")
  }
  if (anyDuplicated(names(data)) != 0) {
    stop("names of data have to be unique!")
  }
  if (any(temp <- !purrr::map_lgl(data, function(x) all(c("features", "origin") %in% names(x))))) {
    stop("every list entry of data has to contain 'features' and 'origin' at least. check: ", paste(names(which(temp)), collapse = ","))
  }
  if (any(temp <- !purrr::map_lgl(data, function(x) is.data.frame(x[["features"]])))) {
    stop("every features entry has to be a data frame. check: ", paste(names(which(temp)), collapse = ","))
  }
  if (any(temp <- !purrr::map_lgl(data, function(x) all(c("range", "complement", "Feature", "value") %in% names(x[["features"]]))))) {
    stop("every features data frame has to have columns named: 'range', 'complement', 'Feature', 'value'. check: ", paste(names(which(temp)), collapse = ","))
  }
  if (any(temp <- !purrr::map_lgl(data, function(x) all(c(length(x[["origin"]]) == 1, is.character(x[["origin"]])))))) {
    stop("every origin entry has to be a character of length 1. check: ", paste(names(which(temp)), collapse = ","))
  }
  if (!is.character(features_to_become_exon)) {
    stop("features_to_become_exon has to be a character vector.")
  }
  if (!is.null(other_features_to_write) && !is.character(other_features_to_write)) {
    stop("other_features_to_write has to be a character vector.")
  }
  if (!is.null(other_features_to_write)) {
    if (any(other_features_to_write %in% features_to_become_exon)) {
      rm <- other_features_to_write[which(other_features_to_write %in% features_to_become_exon)]
      stop(paste(rm, collapse = ", "), " from other_features_to_write were also in features_to_become_exon. Please remove.")
    }
  }
  if (is.null(save_path)) {
    stop("Please provide a save_path.")
  }
  if (!dir.exists(save_path) && file.exists(save_path)) {
    stop("save_path seems to be a file, not a path.")
  }
  if (append && overwrite) {
    stop("append (= attach new lines to existing file) and overwrite (= replacing an existing file), both set to TRUE does not make sense. select one only, please.")
  }

  if (any(!grepl("^##", gtf_header))) {
    stop("each entry in gtf_header has to start with two hashtags. please check.")
  }

  genome_file_path <- file.path(save_path, genome_file_name)
  gtf_file_path <- file.path(save_path, gtf_file_name)
  if (file.exists(gtf_file_path) && !append && !overwrite) {
    stop(gtf_file_path, " exists. It must be appended (append = T) or overwritten (overwrite = T) or another file needs to be written.")
  }
  if (file.exists(genome_file_path) && !append && !overwrite) {
    stop(genome_file_path, " exists. It must be appended (append = T) or overwritten (overwrite = T) or another file needs to be written.")
  }

  # remove < > from range column of features data frames
  rm_temp <- F
  for (i in purrr::map_dfr(data, `[[`, "features")[["range"]]) {
    if (grepl("<|>", i)) {
      message("< and/or > detected in range column of at least one features data frame. Will remove those. The indicate uncertainty of feature boundaries.")
      rm_temp <- T
      break
    }
  }
  if (rm_temp) {
    for (i in names(data)) {
      data[[i]][["features"]][["range"]] <- gsub("<|>", "", data[[i]][["features"]][["range"]])
    }
  }
  # check for numbers, dots and comma only in range column
  for (i in purrr::map_dfr(data, `[[`, "features")[["range"]]) {
    if (!grepl("^[0-9.,]+$", i)) {
      stop("At least one range contains other symbols than digits and dots and commas. This should not be. Please check.")
    }
  }
  # check if complement column is logical
  if (any(temp <- !purrr::map_lgl(data, function(x) all(is.logical(x[["features"]][["complement"]]))))) {
    stop("complement columns of features data frames have to be logical (TRUE or FALSE only). check: ", paste(names(which(temp)), collapse = ","))
  }

  dir.create(save_path, showWarnings = F, recursive = T)
  feat_select <- unique(c(features_to_become_exon, other_features_to_write))

  if (append) {
    linewidths <- table(c(nchar(vroom::vroom_lines(genome_file_path, skip = 1, n_max = 10))))
    linewidth_given <- as.numeric(names(which.max(linewidths)))
    if (linewidth_given != genome_file_linewidth) {
      message("Due to appending, genome_file_linewidth is changed to: ", linewidth_given, ".")
      genome_file_linewidth <- linewidth_given
    }
  }

  lines_to_genome_and_gtf <- purrr::map(stats::setNames(names(data), names(data)), function(x) {

    # allow for other column names, and check above
    features <-
      data[[x]][["features"]] %>%
      dplyr::filter(Feature %in% feat_select)

    original_values <- features$value
    features$value <- gsub("\\.", "", make.names(features$value))

    if (anyDuplicated(features$value)) {
      message(length(which(duplicated(features$value))), " duplicated feature name(s) was/were made altered with trailing numbers to achieve uniqueness.")
      features$value <- make.unique(features$value)
    }

    features$Feature[which(features$Feature %in% features_to_become_exon)] <- "exon"

    # notify when 0 rows remained
    if (nrow(features) > 0) {
      lines_to_gtf <- sapply(1:nrow(features), function(y) {
        paste(x, # Chromosome
              "unknown", # unused
              features[y,"Feature",drop=T], # feature; convert entry in features to
              as.character(min(as.numeric(unlist(strsplit(strsplit(features[y,"range",drop=T], ",")[[1]], "\\.\\."))))), # range start
              as.character(max(as.numeric(unlist(strsplit(strsplit(features[y,"range",drop=T], ",")[[1]], "\\.\\."))))), # range end,
              ".", # unused
              ifelse(features[y,"complement",drop=T], "-", "+"), # when complement, then to minus strand
              ".", # unused
              paste0("gene_id \"", gsub("\\.", "", make.names(features[y,"value",drop=T])), "\"; ", "transcript_id \"", gsub("\\.", "", make.names(features[y,"value",drop=T])), "\";"), # attributes
              sep = "\t")
      })

      return(list(genome = c(paste0(">", x), splitString(data[[x]][["origin"]], n = genome_file_linewidth)),
                  gtf = lines_to_gtf))
    } else {
      message("zero rows in features data frame remained for ", x, ". check 'features_to_become_exon' and 'other_features_to_write' !?")
      return(NULL)
    }
  })


  if (overwrite && file.exists(genome_file_path)) {
    file.remove(genome_file_path)
  }
  if (overwrite && file.exists(gtf_file_path)) {
    file.remove(gtf_file_path)
  }

  # write files
  vroom::vroom_write_lines(unlist(sapply(lines_to_genome_and_gtf, "[", "genome"), use.names = F),
                           genome_file_path,
                           append = append)

  ## gtf header?
  if (append) {
    message("Since gtf file is being appended, the header lines provided are ignored.")
    gtf_header <- NULL
  }
  vroom::vroom_write_lines(c(gtf_header,
                             unlist(sapply(lines_to_genome_and_gtf, "[", "gtf"), use.names = F)),
                           gtf_file_path,
                           append = append)

  message(genome_file_path)
  message(gtf_file_path)

  return(lines_to_genome_and_gtf)
}


splitString <- function(inputString, n) {
  # Split the input string into individual characters
  chars <- strsplit(inputString, '')[[1]]

  # Calculate the number of substrings
  numSubstrings <- ceiling(length(chars) / n)

  # Create a matrix to store substrings
  substrings <- matrix('', nrow = numSubstrings, ncol = n)

  # Populate the matrix with substrings
  for (i in seq_along(chars)) {
    substrings[(i - 1) %/% n + 1, i %% n + 1] <- chars[i]
  }

  # Convert matrix to a list of strings
  result <- apply(substrings, 1, paste, collapse = '')

  # Remove empty strings
  result <- result[result != '']

  return(result)
}
