###################### Quant ###################
library(data.table)
library(tidyverse)


files <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/20120203_OSW_NewLib/", pattern = "peptidesOut", full.names = T)
files

peptides <- lapply(files, fread)

filename <- basename(files)
filename <- c("pepFDR", "pephbFDR", "protFDR")

names(peptides) <- filename

head(peptides[[1]])

pep.wide <- lapply(filename, function(x){
  df <- peptides[[x]]
  df$concentration <- NULL
  df$run_id <- str_extract(df$run_id, "S[0-9]+")
  wide <- dcast(df, peptide_sequence~run_id, drop = FALSE, fill = NaN, sum, value.var = "peptide_intensity")
  wide
})
head(pep.wide[[1]])


write.table(pep.wide[[1]], file = "D:/projects/pca_urine_spectral_lib/data/openswath/20120203_OSW_NewLib/20210826_pepFDR_peptidesWide.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(pep.wide[[2]], file = "D:/projects/pca_urine_spectral_lib/data/openswath/20120203_OSW_NewLib/20210826_pephbFDR_peptidesWide.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(pep.wide[[3]], file = "D:/projects/pca_urine_spectral_lib/data/openswath/20120203_OSW_NewLib/20210826_protFDR_peptidesWide.tsv",
            sep = "\t", row.names = F, quote = F)

