% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_cellranger_outs2.R
\name{read_cellranger_outs2}
\alias{read_cellranger_outs2}
\title{Read  and preprocess the contig_annotation file from 10X Genomics' cellranger output and return a concatenated data frame}
\usage{
read_cellranger_outs2(vdj_outs_path)
}
\arguments{
\item{vdj_outs_path}{named vector of paths to the outs-folder (or any folders containing the necessary file, specified above);
names will be added as "sample"-column to the output data frame}
}
\value{
data frame (cl_long)
}
\description{
This function has been reduced from read_cellranger_outs. Most important change is that CDR3 sequences
from equal clonotype_ids are not shared across the clonotype. In principle this function will only read the
filtered_contig_annotations.csv and filtered_contig.fasta and joins them.
}
