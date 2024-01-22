#' Read  and preprocess the contig_annotation file from 10X Genomics' cellranger output and return a concatenated data frame
#'
#' This function has been reduced from read_cellranger_outs. Most important change is that CDR3 sequences
#' from equal clonotype_ids are not shared across the clonotype. In principle this function will only read the
#' filtered_contig_annotations.csv and filtered_contig.fasta and joins them.
#'
#' @param vdj_outs_path named vector of paths to the outs-folder (or any folders containing the necessary file, specified above);
#' names will be added as "sample"-column to the output data frame
#'
#' @return data frame (cl_long)
#' @export
#'
#' @examples
read_cellranger_outs2 <- function(vdj_outs_path) {


  if (is.null(names(vdj_outs_path))) {
    stop("Please provide names with vdj_outs_path, indicating the sample name.")
  }

  ff <- c("filtered_contig_annotations.csv", "filtered_contig.fasta", "consensus_annotations.csv", "consensus.fasta")

  for (i in vdj_outs_path) {
    if (any(!ff %in% list.files(i))) {
      stop(paste(ff[which(!ff %in% list.files(i))], collapse = ", "), " not found in ", i.)
    }
  }

  filt_cont_ann <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig_annotations\\.csv$", recursive = T, full.names = T)
  filt_cont_fast <- sapply(vdj_outs_path, list.files, pattern = "filtered_contig\\.fasta$", recursive = T, full.names = T)

  cons_ann <- sapply(vdj_outs_path, list.files, pattern = "consensus_annotations\\.csv$", recursive = T, full.names = T)
  cons_fast <- sapply(vdj_outs_path, list.files, pattern = "consensus\\.fasta$", recursive = T, full.names = T)

  contig_annotations <-
    purrr::map_dfr(filt_cont_ann, vroom::vroom, delim = ",", show_col_types = F, col_types = vroom::cols(productive = vroom::col_character()), .id = "sample") %>%
    dplyr::rename("contig_id" = contig_id, "clonotype_id" = raw_clonotype_id, "consensus_id" = raw_consensus_id) %>%
    dplyr::mutate(high_confidence = tolower(high_confidence), full_length = tolower(full_length), productive = tolower(productive), is_cell = tolower(is_cell)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(consensus_id = gsub(paste0(clonotype_id, "_"), "", consensus_id)) %>%
    dplyr::mutate(consensus_id = gsub("_", "", consensus_id)) %>%
    dplyr::mutate(consensus_id = ifelse(consensus_id == "None", "None", paste(clonotype_id, consensus_id, sep = "_"))) %>%
    dplyr::ungroup()

  contig_fasta <-
    purrr::map_dfr(filt_cont_fast, function(x) utils::stack(read_fasta(x)), .id = "sample") %>%
    dplyr::rename("contig_seq" = values, "contig_id" = ind)

  consensus_annotations <-
    purrr::map_dfr(cons_ann, vroom::vroom, delim = ",", show_col_types = F, col_types = vroom::cols(productive = vroom::col_character()), .id = "sample") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(consensus_id = gsub(paste0(clonotype_id, "_"), "", consensus_id)) %>%
    dplyr::mutate(consensus_id = gsub("_", "", consensus_id)) %>%
    dplyr::mutate(consensus_id = paste(clonotype_id, consensus_id, sep = "_")) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, clonotype_id, consensus_id, chain, cdr3, cdr3_nt) %>%
    dplyr::rename("cdr3_consensus" = cdr3, "cdr3_nt_consensus" = cdr3_nt)

  consensus_fasta <-
    purrr::map_dfr(filt_cont_fast, function(x) utils::stack(read_fasta(x)), .id = "sample") %>%
    dplyr::rename("consensus_seq" = values, "contig_id" = ind)


  cl_long <-
    contig_annotations %>%
    dplyr::left_join(contig_fasta, by = c("contig_id", "sample")) %>%
    dplyr::left_join(consensus_fasta, by = c("contig_id", "sample")) %>%
    dplyr::left_join(consensus_annotations, by = c("sample", "chain", "clonotype_id", "consensus_id")) %>%
    dplyr::mutate(barcode = stringr::str_replace(barcode, "-1$", ""))

  return(cl_long)
}




.add_clonotype_names <- function(cl_long,
                                strict_TCR_biology = F,
                                names_to_avoid = NULL) {

  ## this function is not made for general usage yet
  ## several columns are not being checked for presence: V_IMGT, J_IMGT, sample, patient, etc.

  cl_wide_names <-
    cl_long %>%
    ## use sample here, not patient, to avoid grouping by same barcode from different samples!!!
    # or make a column "ID" of sample and barcode -- this could be checked for uniqueness within the function first
    # use collapse.order.fun to have the same order for every transcriptome in case of multi annotations
    dplyr::select(chain, cdr3, sample, barcode) %>%
    tidyr::pivot_wider(names_from = chain, values_from = c(cdr3), values_fn = collapse.order.fun) %>%
    tidyr::nest(data = c(sample, barcode))

  if (strict_TCR_biology) {
    cl_wide__TCR_biology_excluded <-
      cl_wide_names %>%
      dplyr::filter((stringr::str_count(TRB, ",") > 0 | is.na(TRB)) | (stringr::str_count(TRA, ",") > 1 | is.na(TRA))) %>%
      tidyr::unnest(data) %>%
      dplyr::mutate(TRA = strsplit(TRA, ","), TRB = strsplit(TRB, ",")) %>%
      # keep_empty = T is not strictly needed when there are NAs and not empty values; below drop_na is used to remove NAs
      tidyr::unnest(TRA, keep_empty = T) %>%
      tidyr::unnest(TRB, keep_empty = T) %>%
      tidyr::pivot_longer(cols = c(TRA, TRB), names_to = "chain", values_to = "cdr3") %>%
      tidyr::drop_na() %>%
      dplyr::distinct() %>%
      dplyr::mutate(cl_name = NA)
    cl_wide_names <-
      cl_wide_names %>%
      dplyr::filter(stringr::str_count(TRB, ",") == 0 | is.na(TRB)) %>%
      dplyr::filter(stringr::str_count(TRA, ",") <= 1 | is.na(TRA))
  }
  message("Picking random names.")
  cl_wide_names <-
    cl_wide_names %>%
    dplyr::mutate(cl_name = igsc::pick_randomNames(n = nrow(.), max_iter = 10000, names_to_avoid = names_to_avoid)) %>% #
    tidyr::unnest(data)
  #any(cl_wide_names$cl_name %in% names_to_avoid)

  cl_long_join <-
    cl_wide_names %>%
    dplyr::mutate(TRA = strsplit(TRA, ","), TRB = strsplit(TRB, ",")) %>%
    # keep_empty = T is probably strictly not needed when there are NAs and not empty values; below drop_na is used to remove NAs
    tidyr::unnest(TRA, keep_empty = T) %>%
    tidyr::unnest(TRB, keep_empty = T) %>%
    tidyr::pivot_longer(cols = c(TRA, TRB), names_to = "chain", values_to = "cdr3") %>%
    tidyr::drop_na() %>%
    dplyr::distinct() %>%
    dplyr::left_join(cl_long, by = dplyr::join_by(sample, barcode, chain, cdr3)) # to get clonotype_id

  cl_wide <-
    cl_long_join %>%
    dplyr::select(cl_name, chain, barcode, sample, cdr3, clonotype_id) %>%
    tidyr::pivot_wider(names_from = chain, values_from = c(cdr3, clonotype_id),
                       names_sep = "_",
                       values_fn = collapse.unique.fun) %>%
    tidyr::separate(sample, into = c("patient", "replicate", "body_fluid"), remove = F)
  #any(cl_wide$cl_name %in% names_to_avoid)

  ## this is not beautiful but it worked somehow;
  ## toggling strict_TCR_biology between T and F gave comparable results for clonotypes which are not actually affected by 'forbidden' multi-annotations
  message("Collapsing clonotypes.")
  cl_wide_start <- cl_wide

  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide2 <- igsc:::.collapse.clonotypes(cl_wide, split_cols = c("cdr3_TRB", "sample"))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide3 <- igsc:::.collapse.clonotypes(cl_wide2, split_cols = c("cdr3_TRA", "sample"))
  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide4 <- igsc:::.collapse.clonotypes(cl_wide3, split_cols = c("cdr3_TRB", "sample"))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide5 <- igsc:::.collapse.clonotypes(cl_wide4, split_cols = c("cdr3_TRA", "sample"))
  message("---------------------------- cdr3_TRA+cdr3_TRB and patient ----------------------------")
  cl_wide6 <- igsc:::.collapse.clonotypes(cl_wide5, split_cols = c("cdr3_TRA", "cdr3_TRB", "patient"))
  message("---------------------------- cdr3_TRB and sample ----------------------------")
  cl_wide7 <- igsc:::.collapse.clonotypes(cl_wide6, split_cols = c("cdr3_TRB", "sample"))
  message("---------------------------- cdr3_TRA and sample ----------------------------")
  cl_wide8 <- igsc:::.collapse.clonotypes(cl_wide7, split_cols = c("cdr3_TRA", "sample"))
  message("---------------------------- cdr3_TRA+cdr3_TRB and patient ----------------------------")
  cl_wide9 <- igsc:::.collapse.clonotypes(cl_wide8, split_cols = c("cdr3_TRA", "cdr3_TRB", "patient"))
  # run this at very last, once
  message("---------------------------- cdr3_TRA+cdr3_TRB and sample ----------------------------")
  cl_wide9 <- igsc:::.collapse.clonotypes(cl_wide9, split_cols = c("cdr3_TRA", "cdr3_TRB", "sample"))

  # only allow same clonotype name across patients, when TRA and TRB are annotated
  shared_cl_name <-
    cl_wide9 %>%
    dplyr::group_by(cl_name) %>%
    dplyr::summarise(n_pat = dplyr::n_distinct(patient)) %>%
    dplyr::filter(n_pat > 1) %>%
    dplyr::pull(cl_name)

  cl_wide9_split <- split(cl_wide9, cl_wide9$cl_name %in% shared_cl_name)
  cl_wide9_split_TRUE <- split(cl_wide9_split[["TRUE"]], cl_wide9_split[["TRUE"]]$cl_name)
  for (i in 1:length(cl_wide9_split_TRUE)) {
    if (anyNA(cl_wide9_split_TRUE[[i]]$cdr3_TRA) || anyNA(cl_wide9_split_TRUE[[i]]$cdr3_TRB)) {
      cl_wide9_split_TRUE[[i]] <- split(cl_wide9_split_TRUE[[i]], cl_wide9_split_TRUE[[i]]$patient)
    }
  }
  cl_wide9_split_TRUE <- purrr::list_flatten(cl_wide9_split_TRUE)

  new_names <- igsc::pick_randomNames(n = length(cl_wide9_split_TRUE), max_iter = 10000, names_to_avoid = c(cl_wide9$cl_name, names_to_avoid))
  for (i in 1:length(cl_wide9_split_TRUE)) {
    cl_wide9_split_TRUE[[i]]$cl_name <- new_names[i]
  }
  cl_wide9_split_TRUE <- dplyr::bind_rows(cl_wide9_split_TRUE)
  cl_wide9 <- dplyr::bind_rows(cl_wide9_split_TRUE, cl_wide9_split[["FALSE"]])


  # run and notify of result
  # get.clonotype.levels.per.ref
  # get.cdr3.levels.per.clonotype
  # get.clname.levels.per.clonotypeid

  # prep cl_wide for joining
  cl_long9 <-
    cl_wide9 %>%
    dplyr::select(-c(clonotype_id_TRB, clonotype_id_TRA)) %>%
    dplyr::mutate(cdr3_TRA = strsplit(cdr3_TRA, ","), cdr3_TRB = strsplit(cdr3_TRB, ",")) %>%
    ## very important to set keep_empty = T; otherwise empty rows (e.g. no cdr3_TRA) are dropped
    tidyr::unnest(cdr3_TRA, keep_empty = T) %>%
    tidyr::unnest(cdr3_TRB, keep_empty = T) %>%
    tidyr::pivot_longer(cols = c(cdr3_TRA, cdr3_TRB), names_to = "chain", values_to = "cdr3") %>%
    tidyr::drop_na() %>%
    dplyr::distinct() %>%
    dplyr::select(-c(patient, replicate, body_fluid)) %>%
    dplyr::mutate(chain = gsub("cdr3_", "", chain))

  if (strict_TCR_biology) {
    # run get.cdr3.levels.per.clonotype(cl_wide9)
    # then split those clonotypes which have multiple TRB
    multi_TRB <-
      igsc:::.get.cdr3.levels.per.clonotype(cl_wide9) %>%
      dplyr::filter(cdr3_TRB > 1) %>%
      dplyr::pull(cl_name)
    if (length(multi_TRB) > 0) {
      cl_wide9_1 <-
        cl_wide9 %>%
        dplyr::filter(!cl_name %in% multi_TRB) #| is.na(cl_name)
      cl_wide9_2 <-
        cl_wide9 %>%
        dplyr::filter(cl_name %in% multi_TRB)
      split_val <- purrr::reduce(lapply(unique(c("cl_name", "cdr3_TRB")), function(x) cl_wide9_2[,x,drop=T]), paste, sep = "_")
      cl_wide9_2_split <- split(cl_wide9_2, split_val)
      new_names <- pick_randomNames(n = length(cl_wide9_2_split), names_to_avoid = unique(cl_wide9_1[,"cl_name",drop=T], names_to_avoid), max_iter = 10000)
      for (i in seq_along(cl_wide9_2_split)) {
        cl_wide9_2_split[[i]]$cl_name <- new_names[i]
      }
      cl_wide9_2 <- dplyr::bind_rows(cl_wide9_2_split)
      cl_wide9 <- dplyr::bind_rows(cl_wide9_1, cl_wide9_2)
    }
    cl_wide9 <- dplyr::bind_rows(cl_wide9, cl_wide__TCR_biology_excluded)
  }

  cl_long9_join <-
    cl_long %>%
    dplyr::left_join(cl_long9, by = dplyr::join_by(sample, barcode, chain, cdr3))

  message("Making cl_wide from cl_long.")
  cl_wide9_join <-
    cl_long9_join %>%
    dplyr::select(sample, barcode, chain, cl_name,         cdr3, cdr3_nt, contig_id, v_gene, d_gene, j_gene, c_gene, reads, umis, clonotype_id, consensus_id, contig_seq, consensus_seq, V_imgt, J_imgt) %>%
    tidyr::pivot_wider(names_from = chain, values_from = c(cdr3, cdr3_nt, contig_id, v_gene, d_gene, j_gene, c_gene, reads, umis, clonotype_id, consensus_id, contig_seq, consensus_seq, V_imgt, J_imgt), values_fn = collapse.unique.fun) %>% #, values_fn = collapse.unique.fun
    tidyr::separate(sample, into = c("patient", "replicate", "body_fluid"), remove = F)

  message("Ordering columns of cl_wide.")
  # use alphabetical order from TRA and TRB to order other columns
  TRA_order <- sapply(strsplit(cl_wide9_join$cdr3_TRA, ","), order, simplify = F)
  for (colname in grep("_TRA$", names(cl_wide9_join), value = T)) {
    cl_wide9_join[,colname] <- purrr::map2_chr(.x = cl_wide9_join[,colname,drop=T], .y = TRA_order, function(x,y) {
      if (grepl(",", x)) {
        paste(strsplit(x, ",")[[1]][y], collapse = ",")
      } else {
        x
      }
    })
  }
  TRB_order <- sapply(strsplit(cl_wide9_join$cdr3_TRB, ","), order, simplify = F)
  for (colname in grep("_TRB$", names(cl_wide9_join), value = T)) {
    cl_wide9_join[,colname] <- purrr::map2_chr(.x = cl_wide9_join[,colname,drop=T], .y = TRB_order, function(x,y) {
      if (grepl(",", x)) {
        paste(strsplit(x, ",")[[1]][y], collapse = ",")
      } else {
        x
      }
    })
  }
  for (colname in grep("_TRA$|_TRB$", names(cl_wide9_join), value = T)) {
    cl_wide9_join[which(cl_wide9_join[,colname,drop=T] == "NA"),colname] <- NA
  }

  if (strict_TCR_biology) {
    if (length(intersect(cl_wide9_join %>% dplyr::filter(is.na(cl_name)) %>% dplyr::distinct(barcode) %>% dplyr::pull(barcode),
                         cl_long9_join %>% dplyr::filter(is.na(cl_name)) %>% dplyr::distinct(barcode) %>% dplyr::pull(barcode))) !=
        length(cl_wide9_join %>% dplyr::filter(is.na(cl_name)) %>% dplyr::distinct(barcode) %>% dplyr::pull(barcode))) {
      warning("Barcodes associated with is.na(cl_name) are not equal between cl_wide and cl_long.")
    }
  }


  #cl_nest and cl_wide_join should have same nrow
  return(list(cl_wide_start = cl_wide_start,
              cl_long = cl_long9_join,
              cl_wide = cl_wide9_join,
              cl_nest = cl_long9_join %>%
                tidyr::nest(data = c(chain, cdr3, cdr3_nt, contig_id, v_gene, d_gene, j_gene, c_gene,
                                     reads, umis, clonotype_id, consensus_id, contig_seq, consensus_seq, V_imgt, J_imgt, length))))
}



collapse.fun <- function(x) paste(x, collapse = ",")
collapse.unique.fun <- function(x) {
  if (length(unique(x)) == 1) {
    as.character(unique(x))
  } else {
    paste(x, collapse = ",")
  }
}
collapse.order.fun <- function(x) paste(sort(x), collapse = ",")
collapse.unique.order.fun <- function(x) {
  if (length(unique(x)) == 1) {
    as.character(unique(x))
  } else {
    paste(sort(x), collapse = ",")
  }
}

.check.clonotype.id.levels <- function(cl_wide,
                                      id_cols = c("clonotype_id_TRA", "clonotype_id_TRB"),
                                      group_cols = c("sample")) {

  # maybe this function is useless

  if (length(id_cols) != 2) {
    stop("id_cols has to have length 2.")
  }

  ## test if any clonotype_id between TRA and TRB is different
  cl_wide_sub <-
    cl_wide %>%
    tidyr::drop_na(dplyr::all_of(id_cols)) %>%
    dplyr::mutate(!!id_cols[1] := strsplit(!!rlang::sym(id_cols[1]), ","),
                  !!id_cols[2] := strsplit(!!rlang::sym(id_cols[2]), ",")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!id_cols[1] := paste(unique(!!rlang::sym(id_cols[1])), collapse = ","),
                  !!id_cols[2] := paste(unique(!!rlang::sym(id_cols[2])), collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!rlang::syms(c(id_cols[1], group_cols))) %>%
    dplyr::summarise(n = dplyr::n_distinct(!!rlang::sym(id_cols[1])), .groups = "drop")

  return(cl_wide_sub)
}

.collapse.clonotypes <- function(cl_wide,
                                clonotype_col = "cl_name",
                                cdr3_cols = c("cdr3_TRB", "cdr3_TRA"),
                                split_cols = NULL) {
  # e.g. c("clonotype_id_TRB", "cdr3_TRB", "sample")
  # e.g. c("clonotype_id_TRA", "cdr3_TRA", "sample")
  # e.g. c("cdr3_TRA", "cdr3_TRB", "patient")
  #cdr3_col_to_unnest = NULL) {

  if (is.null(split_cols)) {
    stop("Provide a value for split_cols.")
  }
  '  if (!is.null(cdr3_col_to_unnest)) {
    cdr3_col_to_unnest <- match.arg(cdr3_col_to_unnest, c("cdr3_TRB", "cdr3_TRA"))
    if (!cdr3_col_to_unnest %in% split_cols) {
      stop("Shouldnt cdr3_col_to_unnest be in split_cols as well?")
    }
  }'

  if (length(cdr3_cols) != 2) {
    stop("cdr3_cols should be of length 2.")
  }

  #split_cols <- match.arg(split_cols, c("clonotype_id_TRB", "clonotype_id_TRA"))

  ## this collapses cl_name to one based on same values in split_cols
  ## so within groups of split_cols, the same cl_name is assigned
  ## in principle this uses the annotation from 10X (clonotype_id) to group same clonotypes together

  ## this can only work within sample, as clonotype_id may be equal for completely different clonotypes across samples
  ## so when cl_wide contains multiple samples, these have to be split by split_cols
  ## or this could be done outside of this function

  ## run this function separately on cdr3_TRB and cdr3_TRA
  ## this probably has benefits as NAs in one of these columns cannot be used for grouping

  ## maybe instead of the above; unnesting cdr3_TRB and cdr3_TRA is needed to --> extra fun
  ## alternative: just join all additional columns from cl_long in the very end
  ## catch all collapsible clonotypes (think of multi annotations)

  ## some clonotypes collapse to one clonotype_id, already by 10X, when there is one multi-annotated clone
  ## which then links two somewhat different groups of clonotypes with different CDR3s

  n_barcode_start <- length(unique(cl_wide$barcode))
  barcode_count_start <- sort(table(cl_wide$barcode))
  n_unique_cl_start <- length(unique(cl_wide$cl_name))

  '  if (length(cdr3_col) > 1 && unnest_cdr3_col) {
    message("cdr3_col > 1; setting unnest_cdr3_col to FALSE.")
    unnest_cdr3_col <- F
  }

  if (!is.null(cdr3_col_to_unnest)) {
    cl_wide <-
      cl_wide %>%
      dplyr::mutate(!!cdr3_col_to_unnest := strsplit(!!rlang::sym(cdr3_col_to_unnest), ",")) %>%
      tidyr::unnest(!!rlang::sym(cdr3_col_to_unnest))
  }'

  split_val <- purrr::reduce(lapply(split_cols, function(x) cl_wide[,x,drop=T]), paste, sep = "_")
  cl_wide_split <- split(cl_wide, split_val)
  cl_wide_split_NA <- cl_wide_split[which(grepl("NA_NA", names(cl_wide_split)))]
  cl_wide_split_nonNA <- cl_wide_split[which(!grepl("NA_NA", names(cl_wide_split)))] ##only helper
  # one clonotype annotation from 10X per cl_name
  cl_wide_split_1 <- cl_wide_split_nonNA[sapply(cl_wide_split_nonNA, function(x) length(unique(x[,clonotype_col,drop=T])) == 1)]
  # multiple clonotype annotation from 10X per cl_name
  cl_wide_split_multi <- cl_wide_split_nonNA[sapply(cl_wide_split_nonNA, function(x) length(unique(x[,clonotype_col,drop=T])) > 1)]

  other_cdr3_col <- setdiff(cdr3_cols, split_cols)
  if (length(cl_wide_split_multi) > 0) {
    if (length(other_cdr3_col) > 0) {
      split_val <- sapply(cl_wide_split_multi, function (x) {
        #!!!
        # this is a critical step which decides which clonotypes are collapsed
        # how to decide that algorithmically?
        #!!!
        # the order of if clauses below matters!

        other_chains <- unique(x[,other_cdr3_col,drop=T])
        if (length(other_chains) == 1) {
          return(T)
        }
        if (length(which(is.na(other_chains))) == 1 && length(other_chains) == 2) {
          # one NA and other only one level
          # can be collapsed
          return(T)
        }
        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && !any(grepl(",", other_chains))) {
          # more than one other chain but no multi annotation
          return(F)
        }
        if (length(which(is.na(other_chains))) == 0 && length(other_chains) >= 2) {
          # one NA and other only one level
          # can be collapsed
          return(F)
        }
        if (length(other_chains) == 2 && length(which(grepl(",", other_chains))) == 1) {
          return(T)
        }
        other_chains_split <- unlist(strsplit(other_chains, ","))
        if (length(which(is.na(other_chains))) == 1 && length(other_chains[which(!is.na(other_chains))]) == 2 && length(which(grepl(",", other_chains))) == 1 && length(unique(other_chains_split[which(!is.na(other_chains_split))])) == 2) { #
          return(T)
        }
        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && any(grepl(",", other_chains)) && length(unique(other_chains_split[which(!is.na(other_chains_split))])) > 2) { #
          # multiple levels with multi annotation
          return(F)
        }

        other_chains_split2 <- strsplit(other_chains[which(!is.na(other_chains))], ",")
        if (length(which(is.na(other_chains))) == 1 && length(other_chains[which(!is.na(other_chains))]) == 3 && length(which(grepl(",", other_chains))) == 1 && length(unique(other_chains_split[which(!is.na(other_chains_split))])) == 2 && all(unlist(lapply(other_chains_split2[which(lengths(other_chains_split2) == 1)], function(y) y %in% other_chains_split2[[which(lengths(other_chains_split2) == 2)]])))) {
          return(T)
        }

        if (length(unique(other_chains[which(!is.na(other_chains))])) > 1 && !any(grepl(",", other_chains))) {
          return(F)
        }


        stop("New case: check conditions on collapsing clonotypes.")

        '      max <- max(lengths(strsplit(x[,other_cdr3_col,drop=T], ",")))
      unique_other <- unique(unlist(strsplit(x[,other_cdr3_col,drop=T], ",")))
      unique_other <- unique_other[which(!is.na(unique_other))]
      length(unique_other) <= max'
      })
    } else {
      ## when cdr3_cols is cdr3_TRA and cdr3_TRB
      split_val <- rep(T, length(cl_wide_split_multi))
    }

    cl_wide_split_multi_1 <- cl_wide_split_multi[which(split_val)]
    cl_wide_split_multi_2 <- cl_wide_split_multi[which(!split_val)]
  } else {
    cl_wide_split_multi_1 <- NULL
    cl_wide_split_multi_2 <- NULL
  }


  if (length(cl_wide_split_multi_1) > 0) {
    print_table(table(sapply(cl_wide_split_multi_1, function(x) length(unique(x$cl_name)))),
                c("n unique cl_name", "n groups"),
                paste0(length(cl_wide_split_multi_1), " groups with collapsible clonotypes."))


    cl_wide_split_multi_1 <- lapply(cl_wide_split_multi_1, function(x) {

      x[,clonotype_col] <- sort(unique(x[,clonotype_col,drop=T]))[1] # alphabetic first
      '
      # make this an outside option?
      if (!is.null(cdr3_col_to_unnest)) {
        x[,clonotype_col] <- names(sort(table(x[,clonotype_col,drop=T]), decreasing = T))[1] # most frequent first
      } else {
        x[,clonotype_col] <- sort(unique(x[,clonotype_col,drop=T]))[1] # alphabetic first
      }'
      return(x)
    })
  }


  cl_wide <- Reduce(dplyr::bind_rows, list(dplyr::bind_rows(cl_wide_split_NA),
                                           dplyr::bind_rows(cl_wide_split_1),
                                           dplyr::bind_rows(cl_wide_split_multi_1),
                                           dplyr::bind_rows(cl_wide_split_multi_2)))

  n_unique_cl_end <- length(unique(cl_wide$cl_name))
  message("Number of unique clonotypes at before and after: ", n_unique_cl_start, ", ", n_unique_cl_end, ". (", n_unique_cl_end-n_unique_cl_start, ", ", round(((n_unique_cl_end-n_unique_cl_start)/n_unique_cl_start)*100, 2), " %)")

  n_barcode_end <- length(unique(cl_wide$barcode))
  barcode_count_end <- sort(table(cl_wide$barcode))

  # run check functions and notify how many different clonotypes per xx still exist
  check1 <- igsc:::.get.cdr3.levels.per.clonotype(cl_wide = cl_wide)
  print_table(table(check1[which(check1$total > 1), "n_cdr3_TRA_cdr3_TRB",drop=T]),
              c("n comb of cdr3_TRA+cdr3_TRB", "n observation"),
              caption = "CDR3 levels per cl_name")


  check2 <- igsc:::.get.clname.levels.per.clonotypeid(cl_wide)
  print_table(table(dplyr::distinct(check2, cl_names, n)$n),
              c("n unique cl_name per clonotype id", "n observation"),
              caption = "cl_name levels per clonotype_id")

  print("\n")
  print("\n")
  print("\n")

  if (n_barcode_end != n_barcode_start) {
    warning("Number of unique barcodes has changed.")
  }
  if (!identical(barcode_count_start, barcode_count_end)) {
    warning("Count of barcodes has changed.")
  }

  return(cl_wide)
}

.check.clonotype.changes <- function(cl_wide_before,
                                    cl_wide_after,
                                    clonotype_col = "cl_name") {

  combined_cl <-
    stack(table(cl_wide_before[,clonotype_col,drop=T])) %>%
    dplyr::rename("prev" = values) %>%
    dplyr::left_join(stack(table(cl_wide_after[,clonotype_col,drop=T])) %>%
                       dplyr::rename("new" = values), by = "ind") %>%
    dplyr::mutate(new = ifelse(is.na(new), 0, new)) %>% # now missing cl_name (assigned completely to another cl_name) get a 0
    dplyr::mutate(diff = new-prev)

  return(combined_cl)
}

.compare.cl.wide.df <- function(cl_wide1,
                               cl_wide2,
                               ref_cols = c("clonotype_id_TRB", "clonotype_id_TRA", "patient"),
                               clonotype_col = "cl_name") {

  cl_wide1_summ <- igsc:::.get.clonotype.levels.per.ref(cl_wide = cl_wide1,
                                                ref_cols = ref_cols,
                                                clonotype_col = clonotype_col)
  names(cl_wide1_summ)[which(names(cl_wide1_summ) %in% c("cl_names_str", "cl_names"))] <- paste0(names(cl_wide1_summ)[which(names(cl_wide1_summ) %in% c("cl_names_str", "cl_names"))], "1")

  cl_wide2_summ <- igsc:::.get.clonotype.levels.per.ref(cl_wide = cl_wide2,
                                                ref_cols = ref_cols,
                                                clonotype_col = clonotype_col)
  names(cl_wide2_summ)[which(names(cl_wide2_summ) %in% c("cl_names_str", "cl_names"))] <- paste0(names(cl_wide2_summ)[which(names(cl_wide2_summ) %in% c("cl_names_str", "cl_names"))], "2")

  cl_wide_summ <-
    cl_wide1_summ %>%
    dplyr::left_join(cl_wide2_summ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cl1_len = length(cl_names1), cl2_len = length(cl_names2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cl_len_diff = cl1_len - cl2_len)

  return(cl_wide_summ)
}

.get.clonotype.levels.per.ref <- function(cl_wide,
                                         ref_cols = c("clonotype_id_TRB", "clonotype_id_TRA", "patient"),
                                         clonotype_col = "cl_name") {
  #cl_names_count = stack(table(cl_name)))
  cl_wide %>%
    dplyr::group_by(!!!rlang::syms(ref_cols)) %>%
    dplyr::summarise(cl_names_str = paste(unique(cl_name), collapse = ","),
                     cl_names = list(unique(cl_name)), .groups = "drop")
}

.get.cdr3.levels.per.clonotype <- function(cl_wide,
                                          clonotype_col = "cl_name",
                                          cdr3_cols = c("cdr3_TRA", "cdr3_TRB")) {

  # cdr3_cols must be length 2 currently
  ## 3 different comb of cdr3_TRA and cdr3_TRB is possible (NA in TRA and TRB)
  ## 4 is also possible in case of multi annotations

  cl_wide %>%
    tidyr::drop_na(!!rlang::sym(clonotype_col)) %>%
    dplyr::group_by(!!!rlang::syms(clonotype_col)) %>%
    dplyr::summarise(!!paste0("n_", paste(cdr3_cols, collapse = "_")) := dplyr::n_distinct(!!!rlang::syms(cdr3_cols)), # , na.rm = T # triple bang needed to work in n_distinct; like in group_by
                     !!cdr3_cols[1] := dplyr::n_distinct(!!!rlang::syms(cdr3_cols[1]), na.rm = T),
                     !!cdr3_cols[2] := dplyr::n_distinct(!!!rlang::syms(cdr3_cols[2]), na.rm = T),
                     total = dplyr::n())
}

.get.clname.levels.per.clonotypeid <- function(cl_wide,
                                              clonotype_col = "cl_name",
                                              clonotype_id_col = c("clonotype_id_TRB", "clonotype_id_TRA"),
                                              grouping_cols = c("sample")) {

  cl_wide %>%
    tidyr::drop_na(!!rlang::sym(clonotype_col)) %>%
    dplyr::group_by(!!!rlang::syms(c(clonotype_id_col, grouping_cols))) %>%
    dplyr::summarise(cl_names = paste(unique(cl_name), collapse = ","), n = dplyr::n_distinct(!!!rlang::syms(clonotype_col)), .groups = "drop")

}


.compare.cl.wide.by.barcode <- function(cl_wide1, cl_wide2) {
  cl_wide1 %>%
    dplyr::select(barcode, sample, cl_name) %>%
    dplyr::rename("cl_name1" = cl_name) %>%
    dplyr::left_join(cl_wide2 %>%
                       dplyr::select(barcode, sample, cl_name) %>%
                       dplyr::rename("cl_name2" = cl_name))
}



.get.number.of.cell.and.unique.cdr3.per.cl <- function(cl_long) {

  cl_long %>%
    dplyr::group_by(patient, cl_name, chain, cdr3) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(patient, cl_name, chain) %>%
    dplyr::mutate(n_unique_cdr3 = dplyr::n_distinct(cdr3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count_and_ncdr3_greate_1 = count > 1 & n_unique_cdr3 > 1)
}

print_table <- function(stack_table, colnames, caption) {
  temp_tab <- data.frame(stack_table)
  colnames(temp_tab) <- colnames
  print(knitr::kable(temp_tab, format = "simple", caption = caption))
}





