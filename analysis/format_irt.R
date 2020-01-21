#!/usr/env R
# This is an R script for reformatting commercial Biognosys iRT peptides
# into OpenMS and open-readable formats

library(data.table)

data <- "/project/6002011/annieha/pca_urine_spectral_lib/data/irt"

irt_file <- "/project/6002011/annieha/pca_urine_spectral_lib/data/irt/irt_transitions.tsv"

irt <- fread(irt_file)

nrt <- list.files(data, pattern="nrt.tsv", full.names=T)
nrt <- fread(nrt)

targeted_col <- c("PrecursorMz","LibraryIntensity","PrecursorCharge","FragmentType","ProductCharge","FragmentSeriesNumber","PeptideSequence", "transition_group_id","ProteinName")

#> str(irt)
#Classes ‘data.table’ and 'data.frame':  33 obs. of  12 variables:
# $ Q1 monoisotopic                              : num  487 487 487 645 645 ...
# $ Q1 average                                   : num  488 488 488 645 645 ...
# $ Q3                                           : num  860 503 803 800 604 ...
# $ relative intensity (approximate, TSQ-Vantage): num  100 53.1 50.2 100 17 ...
# $ rank                                         : int  1 2 3 1 2 3 1 2 3 1 ...
# $ precursor charge                             : int  2 2 2 2 2 2 2 2 2 2 ...
# $ fragment type                                : chr  "y" "y" "y" "y" ...
# $ fragment charge                              : int  1 1 1 1 1 1 1 1 1 1 ...
# $ fragment number                              : int  8 4 7 8 6 10 8 9 6 8 ...
# $ nominal sequence                             : chr  "LGGNEQVTR" "LGGNEQVTR"
#"LGGNEQVTR" "GAGSSEPVTGLDAK" ...
# $ sequence id                                  : chr  "iRT Kit_a" "iRT Kit_a"
#"iRT Kit_a" "iRT Kit_b" ...
# $ transition id                                : chr  "LGGNEQVTR.y8.1+"
#"LGGNEQVTR.y4.1+" "LGGNEQVTR.y7.1+" "GAGSSEPVTGLDAK.y8.1+" ...


irt <- subset(irt, select = c("Q1 monoisotopic", 
  "relative intensity (approximate, TSQ-Vantage)", 
  "precursor charge", "fragment type", "fragment charge", 
  "fragment number", "nominal sequence", "sequence id", "transition id"))

setnames(irt, c("Q1 monoisotopic", "relative intensity (approximate, TSQ-Vantage)",
  "precursor charge", "fragment type", "fragment charge", "fragment number", 
  "nominal sequence", "sequence id", "transition id"), targeted_col)

setnames(nrt, colnames(nrt), c("PeptideSequence", "PrecursorMz",
  "NormalizedRetentionTime"))
nrt$PrecursorMz <- NULL
irt <- merge(irt, nrt, by = c("PeptideSequence"), all.x=TRUE)
irt$ModifiedPeptideSequence <- irt$PeptideSequence

fragmentlabel <- lapply(1:nrow(irt), function(x){
   f <- paste0(irt$FragmentType[x], irt$FragmentSeriesNumber[x])
   f})
irt$FragmentSeriesNumber <- unlist(fragmentlabel)

irt$ProductMz <- irt$PrecursorMz / irt$ProductCharge

write.table(irt, file=paste0(data, "/irt_assaylib.tsv"), sep = "\t", quote=F,
row.names=F)


