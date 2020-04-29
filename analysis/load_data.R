#!/bin/env R

# This is a script for anaalyzing the differences between runs within 1
# experimment. Identifying the spans between irt peptides and the other common
# peptides
library(viridisLite)
library(data.table)
library(ggplot2)
library(tidyr)
library(VennDiagram)

source("/home/annieha/source/tk_orbitrap_dia/analysis/Utilities.R")

basefolder <- "/project/6002011/annieha/pca_urine_spectral_lib"

# test out the data first

pca_exo <- list.files(paste0(basefolder, "/data/Analyzed_DDA/pca_urine_exo"),
  full.names=T)

if(!exist(pca_exo.msms){
  pca_exo.ev <- fread(pca_exo[1])
  pca_exo.msms <- fread(pca_exo[2])
  pca_exo.pep <- fread(pca_exo[3])
  pca_exo.protein <- fread(pca_exo[4])
}

### Identify iRT peptides and the SUC2 peptides

irt <- list.files(paste0(basefolder, "/data/irt"), pattern="assay", 
  full.names=T)

irt <- fread(irt)

pca_exo.submsms <- subset(pca_exo.msms, select=c("Raw file", "Sequence",
	"Modified sequence", "Proteins", "Charge", "Retention time",
	"Precursor Intensity", "m/z"))

setnames(pca_exo.submsms, colnames(pca_exo.submsms), c("raw",
  "PeptideSequence", "ModifiedPeptideSequence", 
  "Proteins", "PrecursorCharge", "RetentionTime", 
  "Intensity", "PrecursorMz"))

pca_exo.submsms$ModifiedPeptideSequence <- reformat_mods(pca_exo.submsms$ModifiedPeptideSequence)

irt <- subset(irt, select=c("ModifiedPeptideSequence", "NormalizedRetentionTime"))
irt <- irt[!duplicated(ModifiedPeptideSequence)]
pca_exo.msms_irt <- merge(pca_exo.submsms, irt, by="ModifiedPeptideSequence")

raw_files <- pca_exo.submsms$raw %>% unique()

########## Calcute the rsq values between the iRT and the exp #########

setnames(pca_exo.msms_irt, "NormalizedRetentionTime", "iRT")

pca_exo.msms_irt_filtered <- lapply(raw_files, function(x){
  sub <- pca_exo.msms_irt[raw == x]
  sub <- sub[order(-Intensity)]
  sub <- sub[!duplicated(PeptideSequence)]
  sub
})

pca_exo.msms_irt_filtered <- do.call('rbind', pca_exo.msms_irt_filtered)
pca_exo.msms_irt <- pca_exo.msms_irt_filtered %>% subset(select = c("ModifiedPeptideSequence", 
  "raw", "RetentionTime"))

pca_exo.msms_irt <- reshape(pca_exo.msms_irt,idvar = "ModifiedPeptideSequence", timevar="raw",
  direction="wide")

colnames(pca_exo.msms_irt) <- c("Peptide", raw_files)

sd_na <- lapply(1:11, function(x){
  s <- sd(pca_exo.msms_irt[x, 2:156], na.rm = T)
  na <- sum(is.na(pca_exo.msms_irt[x, 2:156]))
  c(s,na)
})

sd_na <- do.call('rbind', sd_na) %>% as.data.frame()

pca_exo.msms_irt$sd <-sd_na$V1
pca_exo.msms_irt$na <- sd_na$V2


############## Look at the common peptides ###########################

pca_exo.msms_sub <- lapply(raw_files, function(x){
  df <- pca_exo.msms[`Raw file` == x]
  df <- df[order(PEP)]
  df <- subset(df, select=c("Raw file", "Proteins", "Gene Names", "Modified sequence",
	 "m/z", "Charge", "Retention time", "PEP", "Precursor Intensity", "Sequence", 
	"Masses", "Intensities"))
  setnames(df, colnames(df), c("raw", "proteins", "genes", "ModifiedPeptideSequence", 
	"PrecursorMz", "Charge", "RetentionTime", "pep", "Intensity", "PeptideSequence", 
	"masses", "intensities"))
  df <- df[!duplicated(ModifiedPeptideSequence, Charge)]
  df
})

pep_id_names <- colnames(pca_exo.pep)
pep_id_names <- pep_id_names[grep("Identification*", pep_id_names)]
pca_exo.pep_id <- subset(pca_exo.pep, select=c("Sequence", pep_id_names))

pca_pep_count <- lapply(pca_exo.msms_sub, function(x){
  df <- x
  count <- df$PeptideSequence %>% unique() %>% length()
  c(count, unique(df$raw))
})

pep_count_run <- do.call('rbind', pca_pep_count)
pep_count_run <- pep_count_run %>% as.data.frame(stringsAsFactors=FALSE)
colnames(pep_count_run) <- c("count", "run") 
pep_count_run$count <- as.numeric(pep_count_run$count)
 

pca_pep <- lapply(pca_exo.msms_sub, function(x){
  df <- x
  sequence <- df$PeptideSequence %>% unique()
  sequence
})

names(pca_pep) <- raw_files

common <- Reduce(intersect, pca_pep)
q <- quantile(pep_count_run$count)
# quantile


# Upper 50% quantile: 7750 peptides
##### Further filter by upper 50% quantile:
run_id <- pep_count_run[count >= q[3]]$run # 78 runs

names(pca_exo.msms_sub) <- raw_files

pca_exo.msms_upper <- pca_exo.msms_sub[run_id]

pca_upper_pep <- pca_pep[run_id]
common <- Reduce(intersect, pca_upper_pep)


upper_pepdf <- lapply(pca_exo.msms_upper, function(x){
  run <- x
  run <- run[PeptideSequence %in% common]
  run <- run[order(pep)]
  run <- run[order(-Intensity)]
  run$ModifiedPeptideSequence <- reformat_mods(run$ModifiedPeptideSequence)
  run <- run[!duplicated(ModifiedPeptideSequence, Charge)]
  run <- subset(run, select=c("raw", "proteins", "genes", "ModifiedPeptideSequence", "PrecursorMz", "Charge", "RetentionTime", "PeptideSequence", "Intensity"))
  run
})

statsdf <- lapply(1:length(run_id), function(x){
  run <- upper_pepdf[[x]]
  df <- subset(run, select=c( "ModifiedPeptideSequence","RetentionTime"))
  setnames(df, "RetentionTime", run_id[x])
  df
})

pep_stats <- Reduce(function(x, y) merge(x, y, by = "ModifiedPeptideSequence", all=TRUE), statsdf)

sd_na <- lapply(1:nrow(pep_stats), function(x){
  s <- sd(pep_stats[x, 2:ncol(pep_stats)], na.rm = T)
  na <- sum(is.na(pep_stats[x, 2:ncol(pep_stats)]))
  c(s,na)
})

sd_na <- do.call('rbind', sd_na) %>% as.data.frame()

pep_sd <- data.frame(ModifiedPeptideSequence = pep_stats$ModifiedPeptideSequence, SD = sd_na$V1, na = sd_na$V2)



