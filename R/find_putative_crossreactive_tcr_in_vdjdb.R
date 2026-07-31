#' Find Putative Cross-Reactive TCRs in VDJdb
#'
#' Identifies groups of similar human T-cell receptor CDR3 sequences in VDJdb
#' that are associated with multiple epitope species. Similarity is measured
#' using Levenshtein distance within each TCR-chain/MHC-class group.
#'
#' Before calculating distances, the function standardizes selected epitope
#' gene and species labels, retains records with an information score greater
#' than one, and removes duplicate annotations.
#'
#' @param vdjdb_df path to tsv or tsv.ghz file or
#'   a data frame containing VDJdb records. It must contain the
#'   columns `Species`, `Gene`, `CDR3`, `MHC.class`, `Info`,
#'   `Epitope.gene`, and `Epitope.species`. Column names are normalized with
#'   [make.names()] before processing.
#'
#' @param max_lvdist A non-negative integer specifying the maximum allowed
#'   Levenshtein distance between CDR3 sequences. Self-matches, with distance
#'   zero, are included during matching and subsequently excluded when they
#'   are the only match. Defaults to `2`.
#'
#' @param nthread A positive integer specifying the number of threads passed to
#'   [stringdist::stringdistmatrix()]. Defaults to `20`.
#' @param substitutionMatrix passed to DECIPHER::DistanceMatrix
#'
#' @return A named list with one element for each observed combination of TCR
#'   chain and MHC class. Each element is a data frame with one row per focal
#'   CDR3 sequence and the following columns:
#'
#'   \describe{
#'     \item{CDR3_1}{The focal CDR3 amino-acid sequence.}
#'     \item{data}{A nested data frame containing matching CDR3 sequences,
#'       their Levenshtein distances, and their VDJdb annotations.}
#'     \item{n_rows}{Number of matched annotation records.}
#'     \item{n_species}{Number of distinct epitope species among the matches.}
#'     \item{avg_lv}{Mean Levenshtein distance across matched records.}
#'     \item{max_lv}{Maximum Levenshtein distance across matched records.}
#'     \item{species}{Comma-separated, alphabetically sorted epitope-species
#'       labels.}
#'     \item{lv0}{Character indicator equal to `"1"` when more than one
#'       distance-zero annotation is present and `"0"` otherwise.}
#'   }
#'
#' @details
#' Records are restricted to `Species == "HomoSapiens"` and `Info > 1`.
#' Selected inconsistencies in `Epitope.gene` and `Epitope.species` are
#' normalized before analysis.
#'
#' The distance matrix is calculated separately for every `Gene`/`MHC.class`
#' combination. Consequently, CDR3 sequences from different chains or MHC
#' classes are not compared.
#'
#' A result represents sequence similarity and annotation diversity, not
#' experimentally established TCR cross-reactivity. Additional biological
#' validation is required.
#'
#' @importFrom rlang %||%
#'
#' @seealso
#' [stringdist::stringdistmatrix()], [tidyr::nest()]
#'
#' @examples
#' \dontrun{
#' vdjdb <- vroom::vroom("vdjdb.tsv.gz", .name_repair = make.names)
#'
#' candidates <- find_putative_crossreactive_tcr_in_vdjdb(
#'   vdjdb_df = vdjdb,
#'   max_lvdist = 2,
#'   nthread = 4
#' )
#'
#' names(candidates)
#' candidates[["TRB_MHCI"]]
#' }
#'
#' @export
find_putative_crossreactive_tcr_in_vdjdb <- function(vdjdb_df,
                                                     max_lvdist = 2,
                                                     nthread = 20,
                                                     substitutionMatrix = "BLOSUM62") {
  vdjdb_df <- "/Volumes/CMS_SSD_2TB/VDJdb/vdjdb_20260725.tsv.gz"
  vdjdbhs <- read_vdjdb_fixed(vdjdb_df) |>
    dplyr::filter(Info>1) |>
    dplyr::select(chain, CDR3, MHC.class, dplyr::starts_with("Epi")) |>
    dplyr::mutate(chain_mhc = paste0(chain, "_", MHC.class)) |>
    dplyr::distinct()

  # names(vdjdbhs)
  # check epitope genes
  # eg <- tibble::enframe(table(vdjdbhs$Epitope.gene)) |> dplyr::arrange(name)

  vdjdbhs <- split(vdjdbhs, vdjdbhs$chain_mhc)

  out <- purrr::map(vdjdbhs, function(x) {

    dist <- stringdist::stringdistmatrix(a = unique(x$CDR3),
                                         b = unique(x$CDR3),
                                         method = "lv",
                                         nthread = nthread,
                                         useNames = T)
    dist_df <- brathering::mat_to_df_long(dist,
                                          rownames_to = "CDR3_1",
                                          colnames_to = "CDR3_2",
                                          values_to = "lv",
                                          row_col_type = "char") |>
      dplyr::filter(dplyr::between(lv, 0,max_lvdist)) # 0 includes self-matches

    # tt <- vdjdbhs$TRB_MHCI |> dplyr::add_count(CDR3)
    # tt2 <- vdjdbhs$TRB_MHCI |> dplyr::distinct(CDR3, Epitope.species) |> dplyr::add_count(CDR3)

    dist_blos <- DECIPHER::DistanceMatrix(Biostrings::AAStringSet(unique(x$CDR3)),
                                          substitutionMatrix = substitutionMatrix)
    rownames(dist_blos) <- unique(x$CDR3)
    colnames(dist_blos) <- unique(x$CDR3)
    dist_blos_df <- brathering::mat_to_df_long(dist_blos,
                                               rownames_to = "CDR3_1",
                                               colnames_to = "CDR3_2",
                                               values_to = substitutionMatrix,
                                               row_col_type = "char")

    # tt <- vdjdbhs$TRB_MHCI |> dplyr::add_count(CDR3)
    # tt2 <- vdjdbhs$TRB_MHCI |> dplyr::distinct(CDR3, Epitope.species) |> dplyr::add_count(CDR3)

    dist_df_nest <- dist_df |>
      dplyr::left_join(dist_blos_df) |>
      dplyr::left_join(x, by = c("CDR3_2" = "CDR3")) |>
      tidyr::nest(data = -CDR3_1)

    dist_df_nest$n_rows <- purrr::map_dbl(dist_df_nest$data, ~nrow(.x) %||% 0)
    dist_df_nest$n_species <- purrr::map_dbl(dist_df_nest$data, ~length(unique(.x$Epitope.species)) %||% 0)
    dist_df_nest$avg_lv <- purrr::map_dbl(dist_df_nest$data, ~mean(.x$lv) %||% 0)
    dist_df_nest$max_lv <- purrr::map_dbl(dist_df_nest$data, ~max(.x$lv) %||% 0)
    dist_df_nest$species <- purrr::map_chr(dist_df_nest$data, ~paste(sort(unique(.x$Epitope.species)), collapse = ",") %||% "0")
    # one CDR3 with lv=0 is the self-hit
    dist_df_nest$lv0 <- purrr::map_chr(dist_df_nest$data, ~as.character(as.numeric(sum(.x$lv == 0)>1)) %||% 0)

    dist_df_nest <- dist_df_nest |>
      dplyr::filter(!(n_rows == 1 & avg_lv == 0))

    return(dist_df_nest)
  })

  return(out)
}



read_vdjdb_fixed <- function(path) {

  if (is.character(path)) {
    df <- vroom::vroom(path, .name_repair = make.names)
  } else {
    names(df) <- make.names(df)
  }

   df |>
    dplyr::filter(Species == "HomoSapiens") |>
    dplyr::rename("chain" = Gene) |>
    ## suggested fixes
    dplyr::mutate(Epitope.gene = ifelse(grepl("gag", Epitope.gene, ignore.case = T) , "GAG", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(grepl("synthetic", Epitope.gene, ignore.case = T) , "synthetic", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(grepl("nef", Epitope.gene, ignore.case = T) , "NEF", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(grepl("yeiH", Epitope.gene, ignore.case = T) , "yeiH", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(grepl("gp160", Epitope.gene, ignore.case = T) , "gp160", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(grepl("Mbp", Epitope.gene, ignore.case = T) , "Mbp", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(Epitope.gene == "IE-1", "IE1", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(Epitope.gene == "MART1", "MLANA", Epitope.gene)) |>
    dplyr::mutate(Epitope.gene = ifelse(Epitope.gene == "INSDRIP", "INS-DRiP", Epitope.gene)) |>
    dplyr::mutate(Epitope.species = ifelse(grepl("synthetic", Epitope.species, ignore.case = T) , "synthetic", Epitope.species)) |>
    dplyr::mutate(Epitope.species = ifelse(grepl("HIV", Epitope.species) , "HIV", Epitope.species))
}
