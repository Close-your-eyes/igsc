wd <- dirname(rstudioapi::getActiveDocumentContext()$path)

library(igsc)

## ---- find fasta and gtf files on disk -----
genomes <- list.files("/Volumes/CMS_SSD_2TB/reference_genomes/viral_ref_genomes/NCBI_original_data",
                      full.names = T,
                      recursive = T,
                      pattern = "\\.fna")
genomes <- grep("GCF", genomes, value = T)
names(genomes) <- basename(dirname(dirname(genomes)))

genes <- list.files("/Volumes/CMS_SSD_2TB/reference_genomes/viral_ref_genomes/NCBI_original_data",
                    full.names = T,
                    recursive = T,
                    pattern = "\\.gtf")
genes <- grep("GCF", genes, value = T)
names(genes) <- basename(dirname(dirname(genes)))

viral_genomes_lens <- purrr::map_int(genomes, ~nchar(igsc::read_fasta(.x)))

new_root_path <- "/Volumes/CMS_SSD_2TB/reference_genomes/viral_ref_genomes/20260319_modified_data"

## ---- prepare GTF files ------
virus_names <- purrr::set_names(names(genes))
attr_rename <- list(c("gene" = "gene_name", "gene_biotype" = "gene_type"),
                    c("product" = "gene_name", "gene_biotype" = "gene_type"),
                    c("gene" = "gene_name", "gene_biotype" = "gene_type"),
                    c("gene" = "gene_name", "gene_biotype" = "gene_type"),
                    c("gene" = "gene_name", "gene_biotype" = "gene_type"))
gene_name_replace = c("glycoprotein" = "",
                      "helicase-primase primase subunit" = "",
                      "helicase-primase subunit" = "",
                      "single-stranded DNA binding protein" = "ssDNAbp",
                      "hypothetical protein" = "HYPPROT")

viral_gtf <- purrr::pmap(
  list(x = virus_names,
       y = attr_rename,
       z = viral_genomes_lens), function(x,y,z) {
         
         gtf <- igsc::read_gtf(file_path = genes[[x]],
                               process_attr_col = F)[["gtf"]]
         
         igsc::process_gtf_attribute_col(gtf = gtf,
                                         attr_rename = y,
                                         genome_length = z,
                                         attr_keep = c("gene_id", "transcript_id", "gene_name", "gene_type"),
                                         attr_as = "kv",
                                         gene_name_prefix = paste0(x, "_"),
                                         gene_name_force = "gene_id",
                                         gene_name_replace = gene_name_replace,
                                         rm_index = T,
                                         rm_exon = T,
                                         exons_only = T,
                                         aggregate_overlapping_exon_ranges = T,
                                         aggregate_exons = T,
                                         check_for_rotation = T,
                                         features_to_exon = c("CDS"))[["gtf"]]
       })


## rotate EBV genome
ebvgenome <- rotate_genome_string(igsc::read_fasta(genomes[["EBV"]]), cut = 144792)
igsc::write_fasta(ebvgenome, file = file.path("/Volumes/CMS_SSD_2TB/reference_genomes/viral_ref_genomes/20260319_modified_data/EBV",
                                              basename(genomes[["EBV"]])))

## copy other genome fasta
for (i in names(genomes)[which(names(genomes) != "EBV")]) {
  file.copy(genomes[i], file.path(new_root_path, i, basename(genomes[i])))
}

## write modified gtf to disk
viral_gtf_new_path <- file.path(new_root_path, names(genes), paste0(basename(dirname(genes)), ".gtf"))
purrr::map2(viral_gtf, viral_gtf_new_path, function(x,y) igsc::write_gtf(gtf_df = x, file = y))


## ---- reduce human reference --------
# top30 frequent:
top30 <- c("TMSB4X", "RPS18", "RPS2", "MALAT1", "RPL41", "RPS6", "RPL13A",
           "EEF1A1", "RPLP1", "RPS19", "RPL13", "RPS14", "RPS12", "RPL26",
           "RPS3", "RPL10", "MT-CO1", "RPS15A", "RPS23", "RPS27", "RPL28",
           "RPL11", "RPS4X", "RPL32", "RPL37A", "RPL23A", "RPL3", "RPL27A",
           "RPL19", "RPL15")
hugtf <- read_gtf(file_path = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz",
                  features = "exon",
                  gene_names = top30,
                  process_attr_col = T)[["gtf"]]

## shorten ref_genome to top30 most frequent genes
genome <- igsc::read_fasta("/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2024-A/fasta/genome.fa.gz")
names(genome) <- sapply(strsplit(names(genome), " "), "[", 1)
out <- make_short_refseq_from_gtf(hugtf, genome, overhang_whole = c(20,20))
out$gtf_df <- make_kv_attr_col(out$gtf_df)
write_gtf(out$gtf_df, file = file.path(new_root_path, "human_genome_top30", "genes_hs_top30.gtf"))
write_fasta(out$refseqs, file = file.path(new_root_path, "human_genome_top30", "genome_hs_top30.fna"))



## ---- combine gtf and genomes ------
## viral genomes and top30 human
genomes_new_path <- list.files(new_root_path,
                               pattern = "\\.fna",
                               full.names = T,
                               recursive = T)
gtf_new_path <- list.files(new_root_path,
                           pattern = "\\.gtf",
                           full.names = T,
                           recursive = T)
combine_gtf_and_genome_for_cellranger(genome_files = genomes_new_path,
                                      gtf_files = gtf_new_path,
                                      save_path = dirname(new_root_path),
                                      save_names = c("genome_hhv15_hs30.fa", "genes_hhv15_hs30.gtf"),
                                      overwrite = T)





## ---- unrelated extra -------
## test output from get_seq_from_refseq
genome <- igsc::read_fasta("/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2024-A/fasta/genome.fa.gz",
                           seqname = "chr6")
hugtf <- read_gtf(file_path = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz",
                  features = "exon",
                  gene_names = top30,
                  process_attr_col = T)[["gtf"]]

names(genome)
gtf_rps18 <- hugtf |> dplyr::filter(gene_name == "RPS18")
# whole range
con1 <- get_seq_from_refseq(gtf_df = gtf_rps18, refseq = genome["chr6 6"],
                            what = "whole")
# single exons
con2 <- get_seq_from_refseq(gtf_df = gtf_rps18, refseq = genome["chr6 6"],
                            what = "single")
# exons grouped
con3 <- get_seq_from_refseq(gtf_df = gtf_rps18, refseq = genome["chr6 6"],
                            what = "single", names = gtf_rps18$exon_number, group_by = gtf_rps18$transcript_name)

# plot wo groups
aln <- pwalign_multi(subject = con1, patterns = con2, order_patterns = T)
aln$plot

# plot with y axis groups
aln2.1 <- pwalign_multi(subject = con1, patterns = con3,
                        compare_seq_df_wide_args = list(
                          change_nonref = T,
                          change_ref = F,
                          rm_pure_NA_non_ref = F))
aln2.1$plot

# remove gaps for plotting only
aln2.2 <- pwalign_multi(subject = con1,
                        patterns = con3,
                        order_patterns = T,
                        compare_seq_df_wide_args = list(
                          change_nonref = T,
                          rm_pure_NA_non_ref = T))
aln2.2$plot



