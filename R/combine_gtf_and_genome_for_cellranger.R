#' Make a custom/hybrid reference for read alignment
#'
#' Genome fasta files and gtf files are concatenated respectively and some
#' checks are run. Pairs of genome_files and gtf_files expected.
#'
#' @param genome_files vector of paths to genome fasta files (.fa, .fna, .fasta)
#' @param gtf_files vector of paths to matching gtf files (e.g. genes.gtf)
#' @param gtf_exons_only remove all lines but exons from gtf files; cellranger
#' only considers lines with Feature = exon; can be a vector to handle gtf files
#' differently
#' @param gtf_features_to_exon list of vectors of length(gtf_files); do not provide a
#' list of length 1 because recycling logic; features to become exons
#' for Prokaryotic or minimal annotations
#' @param fasta_name_two_parts make names in genome file in 10X Genomics style:
#' chr1 1 or chrM MT instead of chr1 and chrM only
#' @param save_path where to save
#' @param save_names save names
#' @param overwrite overwrite existing target files?
#' @param gtf_header header for gtf
#' @param gtf_rm_exon remove original exon entries from gtf, e.g. when CDS are
#' turned into exon
#'
#' @returns
#' @export
#'
#' @examples
#' # hg38 and 5 herpes virus genomes
#' root <- "/Volumes/CMS_SSD_2TB/reference_genomes/"
#' combine_gtf_and_genome_for_cellranger(genome_files = c(paste0(root, "release_114/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
#'                                                        list.files(paste0(root, "GRCh38_2020_A_HHV1to5_appended/viral_ref_genomes"),
#'                                                                   pattern = "fna$", recursive = T, full.names = T)),
#'                                       gtf_files = c(paste0(root, "release_114/Homo_sapiens.GRCh38.114.gtf.gz"),
#'                                                     list.files(paste0(root, "GRCh38_2020_A_HHV1to5_appended/viral_ref_genomes"),
#'                                                                pattern = "gtf$", recursive = T, full.names = T)),
#'                                       save_path = paste0(root, "GRCh38_2020_A_HHV1to5_appended"),
#'                                       save_names = c("genome2.fa", "gtf2.gtf"),
#'                                       gtf_features_to_exon = c(list(c("")), rep(list(c("CDS")), 5)),
#'                                       gtf_rm_exon = c(F, rep(T, 5)))
combine_gtf_and_genome_for_cellranger <- function(genome_files,
                                                  gtf_files,
                                                  gtf_exons_only = T,
                                                  gtf_features_to_exon = list(c("")),
                                                  gtf_rm_exon  = F,
                                                  fasta_name_two_parts = T,
                                                  save_path = dirname(genome_files[1]),
                                                  save_names = c("genome1.fa", "genes1.gtf"),
                                                  overwrite = F,
                                                  gtf_header = c(
                                                    "##description: made with combine_gtf_and_genome_for_cellranger function from https://github.com/Close-your-eyes/igsc",
                                                    paste0("##source files: ", paste(basename(gtf_files), collapse = ", ")),
                                                    "##provider: CMS",
                                                    "##conctact: vonskopnik@pm.me",
                                                    "##format: gtf",
                                                    paste0("##date: ", Sys.Date())
                                                  )
) {

  if (length(genome_files) != length(gtf_files)) {
    stop("genome_files and gtf_files must have same lengths.")
  }

  if (any(!grepl("^##", gtf_header))) {
    stop("each entry in gtf_header has to start with two hashtags. please check.")
  }

  if (!overwrite && (file.exists(file.path(save_path, save_names[1]))) || (file.exists(file.path(save_path, save_names[2])))) {
    stop("at least one target file exists on disk. rename or delete or set overwrite = T or change save_path/save_names.")
  }

  gtf_exons_only <- brathering::recycle(gtf_exons_only, gtf_files, one_or_equal = T)
  #gtf_features_to_exon <- brathering::recycle(gtf_features_to_exon[[1]], gtf_files, one_or_equal = T)) # cannot predict argument format
  gtf_rm_exon <- brathering::recycle(gtf_rm_exon, gtf_files, one_or_equal = T)

  names(genome_files) <- basename(genome_files)
  genome_gtf_list <- purrr::pmap(
    list(fasta_path = genome_files,
         gtf_path = gtf_files,
         gtf_exons_only = gtf_exons_only,
         gtf_features_to_exon = gtf_features_to_exon,
         gtf_rm_exon = gtf_rm_exon),
    process_files,
    fasta_name_two_parts = fasta_name_two_parts
  )
  genome_gtf_list <- purrr::map2(
    genome_files,
    gtf_files,
    process_files,
    gtf_exons_only = gtf_exons_only,
    fasta_name_two_parts = fasta_name_two_parts,
    gtf_features_to_exon = gtf_features_to_exon
  )

  # check duplicates
  fasta_names <- purrr::map(genome_gtf_list, ~names(.x[["genome"]]))
  if (anyDuplicated(unlist(fasta_names))) {
    genome_gtf_fasta_names <<- fasta_names
    stop("duplicate fasta names across files. check global variable genome_gtf_fasta_names.")
  }

  # https://kb.10xgenomics.com/s/article/115003327112-How-can-we-add-genes-to-a-reference-package-for-Cell-Ranger
  # A semicolon-delimited list of key-value pairs of the form key "value".
  # The attribute keys transcript_id and gene_id are required; gene_name is optional, but
  # if present will be preferentially displayed in reports.
  gtfs <- purrr::map(genome_gtf_list, ~process_gtf_attribute_col(.x[["gtf"]],
                                                                 attr_keep = c("gene_id", "transcript_id", "gene_name"),
                                                                 attr_rename = c("gene" = "gene_name"),
                                                                 attr_as = "kv",
                                                                 rm_index = T))
  # replace gtfs in list
  genome_gtf_list <- purrr::map2(genome_gtf_list, gtfs, function(x, y) {
    x[["gtf"]] <- y[["gtf"]]
    return(x)
  })

  gene_ids <- purrr::map(gtfs, ~unique(.x[["attr"]][["gene_id"]]))
  if (anyDuplicated(unlist(gene_ids))) {
    message("duplicate gene_id across gtfs. this could be a problem for aligner.")
  }

  # write combined files
  write_fasta(seqs = purrr::reduce(sapply(genome_gtf_list, "[", 1), c),
              file = file.path(save_path, save_names[1]))

  write_gtf(gtf_df_list = genome_gtf_list,
            header = gtf_header,
            file = file.path(save_path, save_names[2]))

}

process_files <- function(fasta_path,
                          gtf_path,
                          gtf_exons_only = T,
                          fasta_name_two_parts = T,
                          gtf_features_to_exon = c(""),
                          gtf_rm_exon = F) {

  message(basename(fasta_path), " + ", basename(gtf_path))
  gtf1 <- igsc::read_gtf(gtf_path, process_attr_col = F)[["gtf"]]
gtf_nest <- gtf1 |>
  dplyr::distinct(feature)
  if (gtf_rm_exon) {
    gtf1 <- dplyr::filter(gtf1, feature != "exon")
  }
  gtf1$feature[which(gtf1$feature %in% gtf_features_to_exon)] <- "exon"

  if (gtf_exons_only) {
    before <- nrow(gtf1)
    gtf1 <- dplyr::filter(gtf1, feature == "exon")
    if (nrow(gtf1) == 0) {
      stop("no rows found with exon.")
    } else {
      message("all gtf rows: ", formatC(before, big.mark = ","), ". exons only: ", formatC(nrow(gtf1), big.mark = ","), ".")
    }
  }

  genome1 <- igsc::read_fasta(fasta_path)

  namedf <- tidyr::drop_na(match_gtf_fasta_names(fasta = genome1, gtf = gtf1, fasta_name_two_parts = fasta_name_two_parts))
  gtf_map <- stats::setNames(namedf$gtf, namedf$gtf_orig)
  fa_map <- stats::setNames(namedf$fa, namedf$fa_orig)

  gtf1 <- dplyr::filter(gtf1, seqname %in% names(gtf_map))
  gtf1$seqname <- gtf_map[gtf1$seqname]

  genome1 <- genome1[which(names(genome1) %in% names(fa_map))]
  names(genome1) <- fa_map[names(genome1)]

  return(list(genome = genome1, gtf = gtf1, namedf = namedf))
}

match_gtf_fasta_names <- function(fasta,
                                  gtf,
                                  fasta_name_two_parts = T) {
  ## oriented on references as made by 10X Genomics: fasta_name_two_parts

  names_fa <- names(fasta)
  if (anyDuplicated(names_fa)) {
    stop("duplicate names in genome fasta file.")
  }
  names_gtf <- unique(gtf$seqname)


  names_fa_short <- sapply(strsplit(names_fa, " "), "[", 1)
  names(names_fa_short) <- gsub("chr", "", names_fa) # when genome fasta comes from 10X and contains chr already, strip it to make the function below work for all the same
  names_gtf_short <- sapply(strsplit(names_gtf, " "), "[", 1)
  fa_orig_short_map <- stats::setNames(names_fa, names_fa_short)
  gtf_orig_short_map <- stats::setNames(names_gtf, names_gtf_short)
  ## match names_fa_short and names_gtf_short here?

  names_fa_gtf_10X_adj <- list(fa = names_fa_short, gtf = names_gtf_short)
  for (i in names(names_fa_gtf_10X_adj)) {
    # get grouped indices
    chroms_num_inds <- c(which(!is.na(suppressWarnings(as.numeric(names_fa_gtf_10X_adj[[i]])))), which(names_fa_gtf_10X_adj[[i]] %in% c("X", "Y")))
    chroms_MT_inds <- which(names_fa_gtf_10X_adj[[i]] %in% c("MT"))
    chroms_other_inds <- setdiff(seq_along(names_fa_gtf_10X_adj[[i]]), c(chroms_num_inds, chroms_MT_inds))
    # rename fasta names
    names_fa_gtf_10X_adj[[i]][chroms_num_inds] <- paste0("chr", names_fa_gtf_10X_adj[[i]][chroms_num_inds], " ", names_fa_gtf_10X_adj[[i]][chroms_num_inds])
    names_fa_gtf_10X_adj[[i]][chroms_MT_inds] <- "chrM MT"
    names_fa_gtf_10X_adj[[i]][chroms_other_inds] <- paste0(names_fa_gtf_10X_adj[[i]][chroms_other_inds], " ", names_fa_gtf_10X_adj[[i]][chroms_other_inds])
  }
  fa_short_10X_map <- stats::setNames(names_fa_short, names_fa_gtf_10X_adj[["fa"]])
  gtf_short_10X_map <- stats::setNames(names_gtf_short, names_fa_gtf_10X_adj[["gtf"]])

  names_fa_gtf_10X_adj_df <- dplyr::full_join(data.frame(fa = names_fa_gtf_10X_adj[["fa"]], name = names_fa_gtf_10X_adj[["fa"]]),
                                              data.frame(gtf = names_fa_gtf_10X_adj[["gtf"]], name = names_fa_gtf_10X_adj[["gtf"]]),
                                              by = "name")
  names_fa_gtf_10X_adj_df$fa_short <- fa_short_10X_map[names_fa_gtf_10X_adj_df$fa]
  names_fa_gtf_10X_adj_df$fa_orig <- fa_orig_short_map[names_fa_gtf_10X_adj_df$fa_short]
  names_fa_gtf_10X_adj_df$gtf_short <- gtf_short_10X_map[names_fa_gtf_10X_adj_df$gtf]
  names_fa_gtf_10X_adj_df$gtf_orig <- gtf_orig_short_map[names_fa_gtf_10X_adj_df$gtf_short]
  # for gtf: only first half
  names_fa_gtf_10X_adj_df$gtf <- sapply(strsplit(names_fa_gtf_10X_adj_df$gtf, " "), "[", 1)

  if (anyDuplicated(names_fa_gtf_10X_adj_df$name)) {
    names_fa_gtf_10X_adj_df <<- names_fa_gtf_10X_adj_df
    stop("processing creating duplicate names. global variable created: names_fa_gtf_10X_adj_df. check it.")
  }

  ## what if no matches between fasta and gtf? fix it ?

  if (!anyNA(names_fa_gtf_10X_adj_df$fa)) {
    message("all names from fasta appear in gtf.")
  } else {
    miss <- names_fa_gtf_10X_adj_df$name[which(is.na(names_fa_gtf_10X_adj_df$fa))]
    message("names from fasta which do not appear in gtf (", length(miss), "): ", paste(miss, collapse = ", "))
    message("these will be removed from gtf.")
  }
  if (!anyNA(names_fa_gtf_10X_adj_df$gtf)) {
    message("all names from gtf appear in fasta.")
  } else {
    miss <- names_fa_gtf_10X_adj_df$name[which(is.na(names_fa_gtf_10X_adj_df$gtf))]
    message("names from gtf which do not appear in fasta (", length(miss), "): ", paste(miss, collapse = ", "))
    message("these will be removed from fasta.")
  }
  return(names_fa_gtf_10X_adj_df)
}
