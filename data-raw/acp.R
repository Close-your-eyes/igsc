## code to prepare `acp` dataset goes here
acp <- c("A" = "lightpink",
         "T" = "moccasin",
         "C" = "darkseagreen1",
         "G" = "lightblue1",
         "-" = "white",
         "match" = "grey70",
         "mismatch" = "tomato2",
         "gap" = "mediumpurple1",
         "insertion" = "black",
         "ambiguous" = "goldenrod2")
usethis::use_data(acp, overwrite = T, internal = T)
