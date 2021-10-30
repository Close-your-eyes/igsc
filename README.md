
<!-- README.md is generated from README.Rmd. Please edit that file -->

# igsc

<!-- badges: start -->
<!-- badges: end -->

## Installation

Please install packages from bioconductor manually.

``` r
install.packages("BiocManager")
BiocManager::install("DECIPHER")
BiocManager::install("Biostrings")
```

Please install igsc from GitHub. This requires devtools.

``` r
install.packages("devtools")
devtools::install_github("Close-your-eyes/igsc")
```

## What it is

This small package provides a standardized workflow to extract the V(D)J
sequences from CellRangers output (10X Genomics). Only TCRs are handled
at the moment. To gain confidence reference sequences from the
[IMGT](http://www.imgt.org) database are included in an R-friendly
format (data frame) which may be aligned to one’s sequencing data. Based
on this alignment the V(D)J sequence of TRA and TRB are retrieved and
are made available for further processing.  
The data from IMGT have been downloaded in Oct-2021. Using the files
included in this package as template one can also download and process
updated versions from IMGT.

## Vignettes (tutorials)

[prepare TCR
sequences](https://close-your-eyes.github.io/igsc/articles/prepare_tcr_seq_data.html)
