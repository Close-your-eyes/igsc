wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(ggplot2)

im_path <- file.path(wd, "R.export", gsub("-", "", Sys.Date()))
dir.create(im_path, recursive = T, showWarnings = F)

save_path <- file.path(wd, "alignment_results2")
dir.create(save_path, recursive = T, showWarnings = F)

# source the script with functions used herein
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "aln_funs.R"))


### rewrite by using functions from igsc

## ---- get all fastq files (68, 2 from each sample with read1 and read2) -----
fastq_files <- list.files(file.path(wd, "Sequencing_data", "Analysis", "data", "fastq"), pattern = "\\.fastq\\.gz$", full.names = T)
fastq_files <- fastq_files[which(!grepl("undetermined", fastq_files, ignore.case = T))]
fastq_files_df <- data.frame(path = fastq_files, filename = basename(fastq_files),
                             sample = sapply(strsplit(basename(fastq_files), "_"), "[", 1),
                             read_number = sapply(strsplit(basename(fastq_files), "_"), "[", 4))

## some ref seq of pathogens were downloaded manually
local_pathogen_ref <- list.files(file.path(wd, "pathogen_ref_sequences_origin"), full.names = T)
local_pathogen_ref <- local_pathogen_ref[which(!grepl("\\.xlsx|\\.fasta", local_pathogen_ref))]
names(local_pathogen_ref) <- sapply(strsplit(basename(local_pathogen_ref), "__"), "[", 2)
cds_list_local <- purrr::map(local_pathogen_ref, get_cds_from_local_pathogen_ref, .progress = T)

# for training:
# sample_name <- "UrineCMV16"
# reads_paired <- get_paired_reads(fastq_gz_path_r1 = fastq_files_df %>% dplyr::filter(sample == sample_name) %>% dplyr::filter(read_number == "R1") %>% dplyr::pull(path), fastq_gz_path_r2 = fastq_files_df %>% dplyr::filter(sample == sample_name) %>% dplyr::filter(read_number == "R2") %>% dplyr::pull(path), min_len = 40, mc.cores = 2)
# saveRDS(reads_paired, file = file.path(wd, "UrineCMV16_reads.RDS"), compress = F)
# reads_paired <- readRDS(file.path(wd, "UrineCMV16_reads.RDS"))

## ---- get reads from fastq and align to refs ------
urine_sample_names <- setNames(unique(fastq_files_df$sample), unique(fastq_files_df$sample))
for (sample_name in urine_sample_names) {
  print(sample_name)
  # calculate read tables within get_paired_reads? then pass the tables to align_reads and have the loop over ref_seq_list outside?
  temp_path <- file.path(wd, paste0(sample_name, "_reads.RDS"))
  if (!file.exists(temp_path)) {
    reads_paired <- igsc::read_paired_reads(fastq_path_r1 = fastq_files_df %>% dplyr::filter(sample == sample_name) %>% dplyr::filter(read_number == "R1") %>% dplyr::pull(path),
                                            fastq_path_r2 = fastq_files_df %>% dplyr::filter(sample == sample_name) %>% dplyr::filter(read_number == "R2") %>% dplyr::pull(path),
                                            min_len = 40,
                                            mc.cores = 8)
    # saveRDS(reads_paired, file = temp_path, compress = T) maybe zap
  } else {
    message("read RDS")
    reads_paired <- readRDS(temp_path)
  }
  
  ## check igsc fun: read_indices return is NULL?!
  r1_table <- table(reads_paired[["reads"]][["r1"]][["seq"]], useNA = "no")
  r2_table <- table(reads_paired[["reads"]][["r2"]][["seq"]], useNA = "no")
  read_align_list <- igsc::align_reads(r1 = reads_paired[["reads"]][["r1"]][["seq"]],
                                       r2 = reads_paired[["reads"]][["r2"]][["seq"]],
                                       r1_table = r1_table,
                                       r2_table = r2_table,
                                       ref_seq_list = cds_list_local,
                                       mc.cores = 8)
  read_align_list[["ref_seq_list"]] <- cds_list_local
  message("Writing rds file.")
  dir.create(file.path(save_path, sample_name), recursive = T, showWarnings = F)
  saveRDS(read_align_list, file = file.path(save_path, sample_name, "read_align_list.RDS"), compress = T)
}

## ---- read read_align_list.RDS first, then save plots -----
aln_rds_files <- list.files(file.path(wd, "alignment_results2"), pattern = "read_align_list.RDS", recursive = T, full.names = T)
names(aln_rds_files) <- basename(dirname(aln_rds_files))
aln_res <- purrr::map(aln_rds_files, readRDS, .progress = T)


purrr::map(aln_res, function(x) {
  # make this multicore in igsc?
  plots_aligned_reads <-
   igsc:::plot_aligned_reads(match_df_list = purrr::discard(x[["read_matches"]][["match_df_list"]][["r1_r2"]], ~is.null(.x)),
                              ref_seq_list = x[["ref_seq_list"]])
  plots_aligned_reads <-
    purrr::discard(plots_aligned_reads, ~is.null(.x))
  
  message("saving plots to disk.")
  for (i in names(plots_aligned_reads)) {
    dir.create(file.path(wd, "alignment_results", sample_name, gsub(" ", "_", i)), recursive = T, showWarnings = F)
    for (j in names(plots_aligned_reads[[i]])) {
      ggsave(plot = plots_aligned_reads[[i]][[j]][["plot"]],
             path = file.path(wd, "alignment_results", sample_name, gsub(" ", "_", i)),
             filename = paste0("ref_seq_ind_", j, ".png"),
             device = "png",
             dpi = 200,
             width = 6,
             height = min(nlevels(plots_aligned_reads[[i]][[j]][["data"]][["pattern.group"]])/3+1, 49.9))
    }
  }
  
})

## ----- inspect alignment results -------

## read all alignment results
aln_rds_files <- list.files(file.path(wd, "alignment_results"), pattern = "read_align_list.RDS", recursive = T, full.names = T)
names(aln_rds_files) <- basename(dirname(aln_rds_files))
aln_res <- purrr::map(aln_rds_files, readRDS, .progress = T)
#zap::zap_write()

# aln_res[["UrineCMV1"]][["read_matches"]][["match_df_list"]][["r1_r2"]]

# pull out r1_r2 matches
matches <- purrr::map(aln_res, purrr::pluck, "read_matches", "match_df_list", "r1_r2")
# rm NULLs
matches <- purrr::map(matches, purrr::discard, ~is.null(.x))
# bind results across pathogens into one data frame for each sample
matches <- purrr::map(matches, dplyr::bind_rows, .id = "pathogen")
# summaries across read_inds to see whether some reads aligned to multiple pathogens
matches_summary_1 <- parallel::mclapply(matches, function(x) {
    dplyr::summarise(x,
                     n_total_matches = dplyr::n(),
                     n_pathogen = dplyr::n_distinct(pathogen),
                     pathogens = paste(unique(pathogen), collapse = ","), .by = read_ind)
}, mc.cores = 10)


## read results from  Illuminas tool:
illumina_res <-
  openxlsx::read.xlsx(file.path(wd, "Results", "Illumnina_Explify_results", "CMV_urine_1_34_711532164_RPIP_summary.xlsx")) %>%
  dplyr::select(Accession, Detected.Microorganisms) %>%
  dplyr::mutate(Detected.Microorganisms = strsplit(Detected.Microorganisms, "; ")) %>%
  tidyr::unnest(Detected.Microorganisms) %>%
  dplyr::rename("sample" = Accession, "pathogen" = Detected.Microorganisms) %>%
  dplyr::mutate(detected_by_Illumina_app = T) %>%
  dplyr::mutate(pathogen = trimws(gsub(" \\(.{1,}\\)", "", pathogen)))

# prepare for some summarizing plots

# total matched reads
matched_reads_per_sample <- stack(purrr::map(matches_summary_1, nrow))
ggplot(matched_reads_per_sample, aes(x = ind, y = values)) +
  geom_col() +
  ggprism::theme_prism(base_size = 12, base_fontface = NULL, base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  labs(x = "sample", y = "total matched reads")

# number or pathogens per read
matches_summary_1_1 <- dplyr::bind_rows(matches_summary_1, .id = "sample")
ggplot(matches_summary_1_1, aes(x = n_pathogen)) +
  geom_bar() +
  ggprism::theme_prism(base_size = 12, base_fontface = NULL, base_line_size = 0.5) +
  #theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  #labs(x = "sample", y = "total matched reads") +
  ggh4x::facet_wrap2(vars(sample), axes = "all", scales = "free_y")

# which pathogens, count of reads per sample
matches_2 <-
  dplyr::bind_rows(matches, .id = "sample") %>%
  dplyr::group_by(sample, pathogen) %>%
  dplyr::tally(name = "single or paired read matches (n)") %>%
  dplyr::left_join(illumina_res) %>%
  dplyr::mutate(detected_by_Illumina_app = ifelse(is.na(detected_by_Illumina_app), F, detected_by_Illumina_app)) %>%
  dplyr::ungroup()

plotx <- ggplot(matches_2, aes(x = sample, y = `single or paired read matches (n)`)) +
  geom_col(aes(fill = detected_by_Illumina_app)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), legend.title = element_text()) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 5), breaks = c(1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)) +
  #scale_y_log10() +
  facet_wrap(vars(pathogen), scales = "free")
ggsave(plot = plotx, filename = "n_read_per_pathogen.pdf", width = 20, height = 12, device = cairo_pdf, path = im_path)

## paired reads only
matches_2_paired_only <-
  dplyr::bind_rows(matches, .id = "sample") %>%
  tidyr::drop_na() %>%
  dplyr::group_by(sample, pathogen) %>%
  dplyr::tally(name = "single or paired read matches (n)") %>%
  dplyr::left_join(illumina_res) %>%
  dplyr::mutate(detected_by_Illumina_app = ifelse(is.na(detected_by_Illumina_app), F, detected_by_Illumina_app)) %>%
  dplyr::ungroup()

plotx <- ggplot(matches_2_paired_only, aes(x = sample, y = `single or paired read matches (n)`)) +
  geom_col(aes(fill = detected_by_Illumina_app)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), legend.title = element_text()) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 5), breaks = c(1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)) +
  #scale_y_log10() +
  facet_wrap(vars(pathogen), scales = "free")
ggsave(plot = plotx, filename = "n_paired_read_per_pathogen.pdf", width = 20, height = 12, device = cairo_pdf, path = im_path)



# what separates those detected or not detected by Illumina app?
# e.g. Dialister pneumosintes
for (xxx in c("Dialister pneumosintes", "Rothia mucilaginosa")) {
  matches_3 <-
    dplyr::bind_rows(matches, .id = "sample") %>%
    dplyr::filter(pathogen == xxx) %>%
    dplyr::group_by(pathogen, sample, ref_seq_ind) %>%
    dplyr::tally(name = "reads_per_ref_seq_ind") %>%
    dplyr::left_join(illumina_res) %>%
    dplyr::mutate(detected_by_Illumina_app = ifelse(is.na(detected_by_Illumina_app), F, detected_by_Illumina_app)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ref_seq_ind = as.character(ref_seq_ind))
  
  ploty <- ggplot(matches_3, aes(x = ref_seq_ind, y = reads_per_ref_seq_ind)) +
    geom_col(aes(fill = detected_by_Illumina_app)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    facet_grid(cols = vars(sample), rows = vars(pathogen), scales = "free")
  ggsave(plot = ploty, filename = paste0("ref_seq_inds_", xxx, ".pdf"), width = length(unique(matches_3$sample))*2, height = 3, device = cairo_pdf, path = im_path)
}



