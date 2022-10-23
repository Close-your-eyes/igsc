## code to prepare `acp` dataset goes here
acp <- c("A" = "lightpink",
         "T" = "moccasin",
         "U" = "moccasin", # other col?
         "C" = "darkseagreen1",
         "G" = "lightblue1",
         "-" = "white",
         "match" = "grey70",
         "mismatch" = "tomato2",
         "gap" = "mediumpurple1",
         "insertion" = "black",
         "ambiguous" = "goldenrod2")

## credit to ggmsa for collecting the color scales NT and AA
## to give credit color scales are integrated with dependency
scheme_NT <- merge(data.frame(ggmsa:::scheme_NT, rowname = rownames(ggmsa:::scheme_NT)),
                   stats::setNames(utils::stack(acp), c("biostrings", "rowname")),
                   by = "rowname",
                   all = T)
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
usethis::use_data(scheme_AA, scheme_NT, acp, overwrite = T, internal = T)



