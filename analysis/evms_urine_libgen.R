############### EVMS Urine data for library ###############
library(data.table)
library(ggplot2)
library(ggExtra)
library(tidyr)
library(stringr)

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")
basefolder <- "D:/projects/pca_urine_spectral_lib/"

# Load data

pca_urine <- paste0(basefolder, "data/pca_dda/Virginia_eps_urines")

pca_urine.allPep <- list.files(pca_urine, pattern = "peptides.txt", full.names = T) %>% fread()

spectral_lib <- list.files(paste0(basefolder, "data/library"), pattern = "assaylib.tsv", full.names = T)
spectral_lib <- fread(spectral_lib)

# Missing msms.txt and evidence.txt

################ Basic analysis #################

str(pca_urine.allPep)

urine.peptides <- pca_urine.allPep$Sequence %>% unique()
urine.proteins <- pca_urine.allPep$`Gene names` %>% unique()
spectrallib.peptides <- spectral_lib$PeptideSequence %>% unique()
spectrallib.proteins <- spectral_lib$UniProtID %>% unique()
common_peptides <- intersect(urine.peptides, spectrallib.peptides)
# 23967
common_proteins <- intersect(urine.proteins, spectrallib.proteins)
# 2469

AllPep <- union(urine.peptides, spectrallib.peptides)

AllProt <- union(urine.proteins, spectrallib.proteins)

Combine_peptide <- data.frame(PeptideSequence = AllPep, stringsAsFactors = F)
Combine_peptide$DirectEPS <- !(AllPep %in% spectrallib.peptides)
Combine_peptide$Urine <- !(AllPep %in% urine.peptides)
Combine_peptide$Both <- AllPep %in% common_peptides
identification <- lapply(1:nrow(Combine_peptide), function(x){
  r <- Combine_peptide[x,]
  if(r$DirectEPS){
    id <- "directEPS_only"
  }else if (r$Urine){
    id <- "urine_only"
  }else if(r$Both){
    id <- "both"
  }
  id
})
Combine_peptide$identification <- unlist(identification)

cohort_counts <- table(Combine_peptide$identification) %>% as.data.table()
colnames(cohort_counts) <- c("Identification", "Peptides")

Combine_protein <- data.frame(Protein = AllProt, stringsAsFactors = F)
Combine_protein$DirectEPS <- !(AllProt %in% spectrallib.proteins)
Combine_protein$Urine <- !(AllProt %in% urine.proteins)
Combine_protein$Both <- AllProt %in% common_proteins
identification <- lapply(1:nrow(Combine_protein), function(x){
  r <- Combine_protein[x,]
  if(r$DirectEPS){
    id <- "directEPS_only"
  }else if (r$Urine){
    id <- "urine_only"
  }else if(r$Both){
    id <- "both"
  }
  id
})


Combine_protein$identification <- unlist(identification)

cohort_counts$Protein <- table(Combine_protein$identification)
cohort_counts <- melt(cohort_counts, id.vars = c('Identification'))

ggplot(cohort_counts, aes(x = variable, y = value, fill = Identification)) + 
  geom_bar(stat="identity", position = "stack") + theme_bw() + 
  ylab("Count") + xlab("") + geom_text(aes(label = value), position = position_stack(vjust = 0.5), vjust = -0.25)


pca_urine.pepScores <- subset(pca_urine.allPep, select=c("Sequence", "Score", "PEP"))

eps_lib <- ms_unique
eps_lib <- subset(eps_lib, select = c("sequence","score", "pep"))
eps_lib$cohort <- "DirectEPS"

setnames(eps_lib, c("sequence", "score", "pep"), c("Sequence", "EPS_score", "EPS_PEP"))
eps_lib <- eps_lib[order(PEP)]
eps_lib <- eps_lib[!duplicated(Sequence)]
pca_urine.pepScores <- pca_urine.pepScores[order(PEP)]
pca_urine.pepScores <- pca_urine.pepScores[!duplicated(Sequence)]

pca_urine.pepScores$cohort <- "EVMS_urine"
pepScores <- rbind(eps_lib, pca_urine.pepScores)

pepScores <- merge(eps_lib, pca_urine.pepScores, by="Sequence")

ggplot(pepScores, aes(x = Score.x, y = Score.y)) + geom_point(size=0.5) + theme_bw()
ggplot(pepScores, aes(x = PEP.x, y = PEP.y)) + geom_point(size=0.5)

