# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
library(VennDiagram)

#################### Functions for Use ###################

# function for uniform modification tags
reformat_mods <- function(col){
  modifiedSequence <- col
  modifiedSequence <- gsub("_", "", modifiedSequence)
  modifiedSequence <- gsub("(ox)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("(ph)","Phospho", modifiedSequence)
  modifiedSequence <- gsub("C","C(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("(ac)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:4\\)","(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:1\\)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:35\\)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:21\\)","Phospho", modifiedSequence)
  
  column <- modifiedSequence
  return(column)
}

# filtering of data

AllPeptides_filter <- function(Peptidesdf){
  #filters out reverse and potential contaminents
  df <- Peptidesdf
  df <- df[Reverse != "+" | `Potential contaminant` != "+"]
  return(df)
}


peptideCountOverObservations <- function(peptidesdf){
  
}









