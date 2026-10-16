
library(igsc)

# fast settings for hla-typing:
# unique_exon2_3: remove alleles with duplicate sequences here --> fewer combinations in pairwise matching
# rev_comp_minus: no increase in speed but needed to align reads from BAM files which have seq from (-)-strand

hla <- hla_df_from_xml(
  "/Volumes/CMS_SSD_2TB/hla.xml.gz",
  lapply_fun = parallel::mclapply,
  mc.cores = 3,
  unique_exon2_3 = T,
  rm_intron_exon_bounds = T,
  rev_comp_minus = T)

# fast settings for hla-typing:
# provide reference_gtf to filter actual hla-gene regions only
# mclapply + mc.cores: run parallel read matching
# maxmis set to 0: much less matching reads but way quicker for read alignment!
# make_reads_distinct: no major increase as only unique are aligned once. but reduces redundancy.
# pairwise_method to heuristic: reduce tested pairwise combinations massively. groups similar alleles.
# allele_diff reduced to 2 or so: less candidates for pairwise testing.

bam <- "/Volumes/CMS_SSD_2TB/R_scRNAseq/2019_SLE_LN/data/Sequencing_data/hg38_CR7_wo_introns/GEX/Pat4_rep1_blood/outs/possorted_genome_bam.bam"
hlares <- hla_typing_from_bam(
  bam = bam,
  hla_ref = hla,
  reference_gtf = "/Volumes/CMS_SSD_2TB/reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz",
  genes = c("A", "B", "C"),
  lapply_fun = parallel::mclapply,
  make_reads_distinct = T,
  pairwise_method = "heuristic",
  mc.cores = 4,
  maxmis = 0,
  allele_diff = 2)



# write reads for arcasHLA
range <- GenomicRanges::GRanges(seqnames = "chr6",
                                ranges = "28000000..34000000")
reads <- get_bam_reads(bam = bam,
                       granges = range,
                       scores = F,
                       decode_bam_flags = F)
write_fastq(reads)
