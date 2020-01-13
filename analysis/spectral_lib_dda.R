# Script for analyzing and looking at properties of the prostate secretome

library(data.table)
library(ggplot2)
library(stringr)

dda <- list.files("/Volumes/Files/tklab/dda_prostate_library", full.names = T)

DDAstats <- lapply(dda, function(x){
  msms <- list.files(x, pattern="msms.txt", full.names = T)
  if(msms != ""){
    res <- fread(msms)
    res}
})
# not enough RAM

############# Load subset of data #######################################
pca_urine_exo <- list.files(dda[4], pattern = "msms.txt", full.names = T)
pca_urine_exo <- fread(pca_urine_exo, verbose=TRUE)
pca_urine <- list.files(dda[6], pattern = "msms.txt", full.names = T)
pca_urine <- fread(pca_urine)
