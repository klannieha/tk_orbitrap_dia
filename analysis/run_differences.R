#!/bin/env R

# This is a script for anaalyzing the differences between runs within 1
# experimment. Identifying the spans between irt peptides and the other common
# peptides

library(data.table)
library(ggplot2)
library(tidyr)
library(VennDiagram)

source("/home/annieha/source/tk_orbitrap_dia/analysis/Utilities.R")

basefolder <- "/project/6002011/annieha/pca_urine_spectral_lib"

# test out the data first

pca_exo <- list.files(paste0(basefolder, "/data/Analyzed_DDA/pca_urine_exo"),
  full.names=T)

pca_exo.ev <- fread(pca_exo[1])
pca_exo.msms <- fread(pca_exo[2])
pca_exo.pep <- fread(pca_exo[3])
pca_exo.protein <- fread(pca_exo[4])


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

pdf(paste0(basefolder,"/results/pca_urine_exo_irt_all.pdf"),onefile = TRUE)

for(r in raw_files){
  msms <- pca_exo.msms_irt[raw == r]
#  msms <- msms[order(PEP)]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  p <- ggplot(msms, aes(x = NormalizedRetentionTime, y = RetentionTime)) + geom_point(size=1, alpha=0.6)
	geom_abline(slope=1)+ theme_bw()
  print(p)
}
dev.off()


png(paste0(basefolder,"/results/pca_urine_exo_irt_all.png"))
p <- ggplot(pca_exo.msms_irt, aes(x = NormalizedRetentionTime,
  y = RetentionTime, color = raw)) + geom_point(size=1, alpha=0.6, show.legend=F) + 
  geom_line(show.legend=F) +theme_bw()
print(p)
dev.off()

########## Calcute the rsq values between the iRT and the exp #########

rsq <- unlist(lapply(raw_files, function(x){
  msms <- pca_exo.msms_irt[raw == x]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  ln <- lm(NormalizedRetentionTime ~ RetentionTime, msms)
  r <- summary(ln)$r.squared
  r
}))

rsq_tb <- data.frame(r_squared = rsq, run = 1:155)

p <- ggplot(rsq_tb, aes(x = run, y = r_squared)) + geom_point() + theme_bw()
+ ylim(0, 1)
ggsave(file=paste0(basefolder, "/results/pca_urine_exo_irt_rsq.png"), plot =print(p))


lnmod <- lapply(raw_files, function(x){
  msms <- pca_exo.msms_irt[raw == x]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  ln <- lm(iRT ~ RetentionTime, msms)
  ln
})

# Look at the distribution of RT diff

#pdf(paste0(basefolder, "/results/pca_urine_exo_irt_resid.pdf"), onefile=TRUE)

p <- ggplot(pca_exo.msms_irt, aes(x = PeptideSequence, y = RetentionTime)) +
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 30))
ggsave("pca_urine_exo_irt_box.png", plot = print(p))


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

pca_pep_count <- lapply(raw_files, function(x){
  df <- pca_exo.msms_sub[raw == x]
  count <- df$PeptideSequence %>% unique() %>% length()
  c(count, x)
})



pca_pep <- lapply(raw_files, function(x){
  df <- pca_exo.msms_sub[raw == x]
  sequence <- df$PeptideSequence %>% unique()
  sequence
})

common <- Reduce(intersect, pca_pep)



