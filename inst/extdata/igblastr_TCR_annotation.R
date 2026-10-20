wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# dir.create(im_path <- file.path(wd, "R_images"), showWarnings = F)
options(im_path = im_path)

BiocManager::install("igblastr")

igblastr::install_igblast()
igblastr::install_IMGT_germline_db(release = "202639-1")
igblastr::list_germline_dbs(preinstalled.only=F)
igblastr::use_germline_db("IMGT-202639-1.Homo_sapiens.IGH+IGK+IGL")

db_name <- igblastr::install_IMGT_germline_db(
  release = "202639-1",
  organism = "Homo sapiens",
  tcr.db = TRUE
)
igblastr::show_intdata_disagreements("IMGT-202631-1.Homo_sapiens.TRA+TRB+TRG+TRD")


igblastr::use_germline_db("IMGT-202639-1.Homo_sapiens.TRA+TRB+TRG+TRD")
ann <- igblastr::igblastn(
  query = file.path(wd, "tcr_test.fasta"),
  ig_seqtype = "TCR",
  domain_system = "imgt",
  c_region_db = NULL,
  outfmt = "AIRR",
  out = file.path(wd, "tcr_annotations.tsv")
)

ann
as.character(igblastr::extract_sequence_region(ann, region_type = "cdr3", as.aa = T))
# igblastr uses the IMGT/AIRR definition of CDR3, which excludes the conserved starting cysteine (C) and ending phenylalanine/tryptophan (F/W). Those residues belong to FR3 and FR4, respectively. The AIRR junction includes both flanking residues.
# https://docs.airr-community.org/en/stable/datarep/rearrangements.html#junction-versus-cdr3

# To get the sequence including the starting C, use:
ann$junction_aa



# ---- Mapping between igblast output and consensus_annotations.csv from 10X -------

  # |reference    |annotations  |notes                                                  |
  # |:------------|:------------|:------------------------------------------------------|
  # |consensus_id |sequence_id  |Normalize identifier: _consensus1 becomes _consensus_1 |
  # |clonotype_id |sequence_id  |Extract clonotype prefix from sequence_id              |
  # |length       |sequence     |Count nucleotides in sequence                          |
  # |chain        |locus        |Direct match                                           |
  # |v_gene       |v_call       |IgBLAST includes allele suffixes                       |
  # |d_gene       |d_call       |Assignments may differ or include multiple candidates  |
  # |j_gene       |j_call       |IgBLAST includes allele suffixes                       |
  # |c_gene       |NA           |Not present in the supplied IgBLAST TSV                |
  # |full_length  |complete_vdj |Related flag; criteria are not necessarily identical   |
  # |productive   |productive   |Reference uses true/false; IgBLAST uses T/F            |
  # |reads        |NA           |Read-count metadata absent from the IgBLAST TSV        |
  # |umis         |NA           |UMI-count metadata absent from the IgBLAST TSV         |
  #
  #
  # |reference_aa |annotations_aa |reference_nt |annotations_nt |notes                                                                                                                             |
  # |:------------|:--------------|:------------|:--------------|:---------------------------------------------------------------------------------------------------------------------------------|
  # |fwr1         |fwr1_aa        |fwr1_nt      |fwr1           |Direct match                                                                                                                      |
  # |cdr1         |cdr1_aa        |cdr1_nt      |cdr1           |Direct match                                                                                                                      |
  # |fwr2         |fwr2_aa        |fwr2_nt      |fwr2           |Direct match                                                                                                                      |
  # |cdr2         |cdr2_aa        |cdr2_nt      |cdr2           |Direct match                                                                                                                      |
  # |fwr3         |fwr3_aa        |fwr3_nt      |fwr3           |Remove final C from amino acids and final three nucleotides                                                                       |
  # |cdr3         |junction_aa    |cdr3_nt      |junction       |Reference CDR3 matches IgBLAST junction, including flanking C and F                                                               |
  # |fwr4         |fwr4_aa        |fwr4_nt      |fwr4           |Remove initial F from amino acids. For nucleotides, remove initial three nucleotides and append the next nucleotide from sequence |
  #
  #
  # |reference  |annotations    |conversion_to_reference |
  # |:----------|:--------------|:-----------------------|
  # |fwr1_start |fwr1_start     |Subtract one            |
  # |fwr1_end   |fwr1_end       |Unchanged               |
  # |cdr1_start |cdr1_start     |Subtract one            |
  # |cdr1_end   |cdr1_end       |Unchanged               |
  # |fwr2_start |fwr2_start     |Subtract one            |
  # |fwr2_end   |fwr2_end       |Unchanged               |
  # |cdr2_start |cdr2_start     |Subtract one            |
  # |cdr2_end   |cdr2_end       |Unchanged               |
  # |fwr3_start |fwr3_start     |Subtract one            |
  # |fwr3_end   |fwr3_end       |Subtract three          |
  # |cdr3_start |cdr3_start     |Subtract four           |
  # |cdr3_end   |cdr3_end       |Add three               |
  # |fwr4_start |fwr4_start     |Add two                 |
  # |fwr4_end   |j_sequence_end |Unchanged               |
  #
  #
  # |reference   |closest_annotations |notes                                                           |
  # |:-----------|:-------------------|:---------------------------------------------------------------|
  # |v_start     |v_sequence_start    |Different alignment boundaries; no simple coordinate conversion |
  # |v_end       |v_sequence_end      |Alignment endpoints may differ                                  |
  # |v_end_ref   |v_germline_end      |Reference sequences and offsets differ                          |
  # |j_start     |j_sequence_start    |Alignment starts may differ beyond coordinate convention        |
  # |j_start_ref |j_germline_start    |Reference alignment starts may differ                           |
  # |j_end       |j_sequence_end      |Numerically identical in all seven supplied rows                |


#---- conversion ----
ann <- tcr_annotations

# Remove allele suffixes while preserving distinct gene candidates.
gene_names <- function(x) {
  vapply(strsplit(as.character(x), ",", fixed = TRUE), function(z) {
    paste(unique(sub("\\*.*$", "", z)), collapse = ",")
  }, character(1))
}

ann_ref_format <- data.frame(
  clonotype_id = sub("_consensus_.*$", "", ann$sequence_id),
  consensus_id = sub(
    "_consensus_", "_consensus", ann$sequence_id, fixed = TRUE
  ),
  length = nchar(ann$sequence),
  chain = ann$locus,
  v_gene = gene_names(ann$v_call),
  d_gene = gene_names(ann$d_call),
  j_gene = gene_names(ann$j_call),
  c_gene = NA_character_,
  full_length = NA,  # complete_vdj is not an exact equivalent
  productive = tolower(as.character(as.logical(ann$productive))),
  reads = NA_integer_,
  umis = NA_integer_,

  # No equivalent leader-inclusive V start or reference offset.
  v_start = NA_integer_,
  v_end = ann$v_sequence_end,
  v_end_ref = NA_integer_,

  # Convert IgBLAST J alignment coordinates to reference numbering.
  j_start = ann$j_sequence_start - 1L,
  j_start_ref = ann$j_germline_start - 1L,
  j_end = ann$j_sequence_end,
  stringsAsFactors = FALSE
)

# Transfer region sequences and convert coordinate numbering.
regions <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")

for (region in regions) {
  ann_ref_format[[region]] <- ann[[paste0(region, "_aa")]]
  ann_ref_format[[paste0(region, "_nt")]] <- ann[[region]]

  ann_ref_format[[paste0(region, "_start")]] <-
    ann[[paste0(region, "_start")]] - 1L

  ann_ref_format[[paste0(region, "_end")]] <-
    ann[[paste0(region, "_end")]]
}

# FR3: move the terminal cysteine into the reference CDR3.
ann_ref_format$fwr3 <- substr(
  ann$fwr3_aa, 1L, nchar(ann$fwr3_aa) - 1L
)
ann_ref_format$fwr3_nt <- substr(
  ann$fwr3, 1L, nchar(ann$fwr3) - 3L
)
ann_ref_format$fwr3_end <- ann$fwr3_end - 3L

# Reference CDR3 includes the conserved flanking residues.
ann_ref_format$cdr3 <- ann$junction_aa
ann_ref_format$cdr3_nt <- ann$junction
ann_ref_format$cdr3_start <- ann$cdr3_start - 4L
ann_ref_format$cdr3_end <- ann$cdr3_end + 3L

# FR4: exclude the initial phenylalanine.
ann_ref_format$fwr4 <- substring(ann$fwr4_aa, 2L)
ann_ref_format$fwr4_start <- ann$fwr4_start + 2L

# Reference FR4 extends through the J alignment endpoint.
ann_ref_format$fwr4_end <- ann$j_sequence_end
ann_ref_format$fwr4_nt <- mapply(
  function(sequence, start, end) {
    if (anyNA(c(sequence, start, end))) return(NA_character_)
    substr(sequence, start + 1L, end)
  },
  sequence = ann$sequence,
  start = ann_ref_format$fwr4_start,
  end = ann_ref_format$fwr4_end,
  USE.NAMES = FALSE
)

tcr_annotations_ref <- vroom::vroom(file.path(wd, "tcr_annotations_ref.csv"))
# Match the reference column order.
ann_ref_format <- ann_ref_format[
  , names(tcr_annotations_ref), drop = FALSE
]


waldo::compare(ann_ref_format, tcr_annotations_ref, max_diffs = Inf)


# 10X algorithm:
# https://github.com/10XGenomics/cellranger/blob/main/lib/rust/vdj_ann/src/annotate.rs
