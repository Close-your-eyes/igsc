## code to prepare `acp` dataset goes here
acp <- c("A" = "lightpink",
         "T" = "moccasin",
         "U" = "moccasin", # other col?
         "C" = "darkseagreen1",
         "G" = "lightblue1",
         "match" = "grey80",
         "." = "grey80",
         "mismatch" = "tomato2",
         "x" = "tomato2",
         "gap" = "plum1",
         "-" = "plum1",
         "insertion" = "grey30",
         "ambiguous" = "cornsilk2")

## credit to ggmsa for collecting the color scales NT and AA
## to give credit color scales are integrated with dependency

if (!requireNamespace("BiocManager", quietly = T)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ggmsa", quietly = T)) {
  BiocManager::install("ggmsa")
}
if (!requireNamespace("Biostrings", quietly = T)) {
  BiocManager::install("Biostrings")
}
if (!requireNamespace("dplyr", quietly = T)) {
  utils::install.packages("dplyr")
}
if (!requireNamespace("usethis", quietly = T)) {
  utils::install.packages("usethis")
}
if (!requireNamespace("devtools", quietly = T)) {
  utils::install.packages("devtools")
}
if (!requireNamespace("HDMD", quietly = T)) {
  devtools::install_github("cran/HDMD")
}

scheme_NT <- merge(data.frame(ggmsa:::scheme_NT, rowname = rownames(ggmsa:::scheme_NT)),
                   stats::setNames(utils::stack(acp), c("biostrings", "rowname")),
                   by = "rowname",
                   all = T)
scheme_NT[1,2:5] <- NA # have the gap replace by our own defined color
for (i in 1:nrow(scheme_NT)) {
  for (j in 1:ncol(scheme_NT)) {
    if (is.na(scheme_NT[i,j])) {
      scheme_NT[i,j] <- scheme_NT[i,"biostrings"]
    }
  }
}
rownames(scheme_NT) <- scheme_NT$rowname
scheme_NT <- scheme_NT[,-which(colnames(scheme_NT) == "rowname")]
for (i in unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N"))[!unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET, "N")) %in% rownames(scheme_NT)]) {
  scheme_NT <- rbind(scheme_NT, "grey90")
  rownames(scheme_NT)[nrow(scheme_NT)] <- i
}

scheme_AA <- rbind(ggmsa:::scheme_AA, "*" = do.call(grDevices::rgb, c(maxColorValue = 255, as.list(grDevices::col2rgb("grey70")[,1]))))
for (i in unique(c(Biostrings::AA_ALPHABET, "N"))[!unique(c(Biostrings::AA_ALPHABET, "N")) %in% rownames(scheme_AA)]) {
  scheme_AA <- rbind(scheme_AA, "grey90")
  rownames(scheme_AA)[nrow(scheme_AA)] <- i
}


scheme_AA <- as.matrix(scheme_AA)
scheme_NT <- as.matrix(scheme_NT)

## aa properties
aa_polar <- c("S","T","N","Q", "C")
aa_pos_charge <- c("K","R","H")
stop_codon <- c("*")
aa_non_polar <- c("A","V","L","I","M","G","P")
aa_neg_charge <- c("D","E")
aromatic_non_polar <- c("W","F","Y")
aa_main_prop <- c(stats::setNames(rep("polar", length(aa_polar)), aa_polar),
                  stats::setNames(rep("non polar", length(aa_non_polar)), aa_non_polar),
                  stats::setNames(rep("basic pos charge", length(aa_pos_charge)), aa_pos_charge),
                  stats::setNames(rep("acidic neg charge", length(aa_neg_charge)), aa_neg_charge),
                  stats::setNames(rep("stop", length(stop_codon)), stop_codon),
                  stats::setNames(rep("non polar aromtic", length(aromatic_non_polar)), aromatic_non_polar))

## from Peptides::aaComp
aa_list <- list(tiny = c("A", "C", "G", "S", "T"),
                small = c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V"), # tiny included
                aliphatic = c("A", "I", "L", "V"),
                aromatic = c("F", "H", "W", "Y"),
                nonpolar = c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"),
                polar = c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"),
                charged = c("B", "D", "E", "H", "K", "R", "Z"), ## all that are basic or acidic
                basic = c("H", "K", "R"),
                acidic = c("B", "D", "E", "Z"))

# https://www.cryst.bbk.ac.uk/education/AminoAcid/the_twenty.html
#Sometimes it is not possible two differentiate two closely related amino acids, therefore we have the special cases:
#asparagine/aspartic acid - asx - B
#glutamine/glutamic acid - glx - Z

aa_df_long <- stack(aa_list)
aa_df_long$ind <- as.character(aa_df_long$ind)
names(aa_df_long) <- c("aa", "property")
aa_df_nest <- dplyr::summarise(dplyr::group_by(aa_df_long, aa), properties = list(property), .groups = "drop")
aa_df_nest2 <-  dplyr::summarise(dplyr::group_by(aa_df_nest, properties), aa = list(aa))
aa_info <- list(aa_main_prop = aa_main_prop,
                aa_list = aa_list,
                aa_df_long = aa_df_long,
                aa_df_nest = aa_df_nest,
                aa_df_nest2 = aa_df_nest2,
                aa_names = utils::stack(Biostrings::AMINO_ACID_CODE[Biostrings::AA_STANDARD]))

# chem_col <- stack(igsc:::scheme_AA[,"Chemistry_AA"]) %>% dplyr::group_by(values) %>% dplyr::summarise(aa = paste(ind, collapse = ","))

usethis::use_data(scheme_AA, scheme_NT, acp, aa_info,
                  overwrite = T, internal = T)



