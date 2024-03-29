% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attach_imgt_allele.R
\name{attach_imgt_alleles}
\alias{attach_imgt_alleles}
\title{Find best matching IMGT reference alleles for TCR V and J gene segments}
\usage{
attach_imgt_alleles(
  seq,
  imgt_ref,
  pick.by = "alignment",
  lapply_fun = lapply,
  seq_col = "consensus_seq_cr",
  V_col = "V_cr",
  J_col = "J_cr",
  ...
)
}
\arguments{
\item{seq}{data frame with clonotype data in long format, e.g. cl_long preferentially from igsc::read_cellranger_outs or
a named vector of consensus sequences, named by V and J segments, separated by "\emph{\strong{"; e.g. stats::setNames(object = cl_long$consensus_seq_cr, nm = paste0(cl_long$V_cr, '}}', cl_long$J_cr))}

\item{imgt_ref}{IMGT reference data frame created with igsc::imgt_tcr_segment_prep or a named vector of sequences and respective allele names;
e.g. stats::setNames(imgt_ref$seq.nt, imgt_ref$Allele)}

\item{pick.by}{match disambiguate alleles from IMGT by best match (alignment) or just randomly (random);
random is intended to speed up the process for testing or when exact IMGT alleles are not necessary}

\item{lapply_fun}{function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested}

\item{...}{additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen}
}
\value{
a data frame of unique seq entries in column
}
\description{
Entries from seq are made unique and are matched against reference alleles in imgt_ref.
If names are provided to seq (or respective V and J columns in case of a data frame) matching will be
quicker as possibles alleles are initially narrowed down by string matching.
}
\examples{
\dontrun{
imgt_ref <- readRDS(system.file("extdata", "IMGT_ref/human/hs.rds", package = "igsc"))
ata <- attach_imgt_alleles(seq = cl_long, imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)
cl_long <-
dplyr::left_join(cl_long, ata, by = c("consensus_seq_cr" = "seq")) \%>\%
tidyr::separate(VJ_IMGT, into = c("V_imgt", "J_imgt"), sep = "___")

## or pass vector, named or unnamed
ata <- attach_imgt_alleles(seq = cl_long$consensus_seq_cr, imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)

ata <- attach_imgt_alleles(seq = stats::setNames(object = cl_long$consensus_seq_cr, nm = paste0(cl_long$V_cr, '___', cl_long$J_cr)),
imgt_ref = imgt_ref, pick.by = "random", lapply_fun = lapply)
}
}
