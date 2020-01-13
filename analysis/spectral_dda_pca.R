# !/usr/env/ R

library(ggplot2)
library(data.table)
library(stringr)
library(VennDiagram)
library(reshape2)
library(matrixStats)

############### Utility functions ############
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

############### Load data #####################


irt <- list.files('../dia_pipline/library/irt', pattern="iRTassays.tsv", full.names=T)
pca <- list.files('../../dda_prostate_library', pattern="pca_", full.names=T)

#Allmsms <- lapply(pca, function(x){
#  sample <- basename(x)
#  msms <- list.files(x, pattern="msms.txt", full.names=T)
#  if(msms != ""){
#    msms <- fread(msms)
#    msms[, sample_name := sample]
#    msms}
#})
irt <- fread(irt, fill=TRUE)
# Try out with one sample first

pca_cell_line_exo <- fread(list.files(pca[1], pattern="msms.txt", full.name=T))
pca_cell_line_exo_pep <- fread(list.files(pca[1], pattern="peptides.txt", full.name=T))
head(pca_cell_line_exo)
pca_cell_line_exo <- pca_cell_line_exo[!grep("*_Blank_*", `Raw file`)]

raw <- unique(pca_cell_line_exo$`Raw file`)
length(raw)
#17 raw files
pca_pep_col <- colnames(pca_cell_line_exo_pep)
pca_pep_col <- pca_pep_col[grep("^Experiment(.*)", pca_pep_col)]
pca_pep_col <- pca_pep_col[!grepl("Blank", pca_pep_col)]
pep_exp <- subset(pca_cell_line_exo_pep, select=c(pca_pep_col))
pep_check <- rowSums(is.na(pep_exp))
sum(pep_check == 0)
head(pep_exp)

msms <- lapply(1:length(raw), function(x){
  res <- pca_cell_line_exo[`Raw file` == raw[x]]
  res <- res[PEP <= 0.01]
  res <- res[order(PEP)]
})

msms_peptides <- lapply(msms, function(x){
  res <- subset(x, select=c("Modified sequence", "Raw file", "Charge", "Precursor Intensity", "Retention time"))
  colnames(res)
  pep <- unique(x$Sequence)
  pep
})

length(Reduce(intersect, msms_peptides))

msms_subset <- lapply(msms, function(x){
  n
})

################ Look at the iRTs in the run ################
irt <- subset(irt, select=c("Tr_recalibrated", "FullUniModPeptideName", "PrecursorCharge"))
colnames(irt) <- c("irt", "sequence", "charge")
irt <- irt[!duplicated(sequence,charge)]
irt <- irt[!is.na(sequence)]
irt$sequence <- reformat_mods(irt$sequence)
pca_cell_line_exo <- fread(list.files(pca[1], pattern="msms.txt", full.names=T))

pca_cell_line_exo <- pca_cell_line_exo[PEP <= 0.01]
pca_cell_line_exo <- pca_cell_line_exo[order(PEP, `Raw file`)]

ms <- subset(pca_cell_line_exo, select=c("Raw file", "Retention time", "Modified sequence", "Charge"))

raw <- unique(ms$raw)
colnames(ms) <- c("raw", "rt", "sequence", "charge")
#ms <- ms[!duplicated(raw,sequence, charge)]
ms$sequence <- reformat_mods(ms$sequence)
ms <- lapply(raw, function(x){
  s <- ms[raw == x]
  s <- s[!duplicated(sequence,charge)]
  s
})

#ms$sequence <- reformat_mods(ms$sequence)

msms_irt <- merge(ms, irt, by=c("sequence", "charge"))

ggplot(msms_irt, aes(x = rt, y = irt, color=raw)) + geom_point(size=1, alpha=.8) + theme_light(base_size=16)

# Check the overlaps first



#### Load one more just to see ###########################
pca_urine <- fread("../../dda_prostate_library/pca_urine/peptides.txt")
pca_urine_col <- colnames(pca_urine)
pca_urine_protein <- fread("../../dda_prostate_library/pca_urine/proteinGroups.txt")
pca_urine_protein <- pca_urine_protein[Peptides >= 2]

pca_urine_msms <- fread("../../dda_prostate_library/urine/msms.txt") # different runs
# Filter the ones that are contaminants, reverse
pca_urine <- pca_urine[Reverse != "+" | `Potential contaminant` != "+"]
sum(duplicated(pca_urine$Sequence))

pca_urine_count <- subset(pca_urine, select=pca_urine_col[grep("Identification", pca_urine_col)])
pca_urine_count[pca_urine_count == "By MS/MS"] <- 1
pca_urine_count[pca_urine_count == "" ] <- 0
pca_urine_count[pca_urine_count == "By matching"] <- 0
pca_urine_count <- as.data.frame(lapply(pca_urine_count, function(x) x <- as.numeric(x)))
pca_urine_count$count <- rowSums(pca_urine_count)
pca_urine_count <- pca_urine_count[pca_urine_count$count != 0,]
pca_counts <- as.data.frame(table(pca_urine_count$count))

setnames(pca_counts, "Var1", "Observations")
ggplot(pca_counts,aes(x = Observations, y = Freq)) + theme_classic() + geom_bar(stat="identity", position = position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

test <- pca_urine_count[1:10]
test[test != ""] <- 1
test[test == ""] <- 0
test <- as.data.frame(lapply(test, function(x) x <- as.numeric(x)))


test <- replace


pca_urine_sub <- subset(pca_urine_msms, 
                        select=c("Raw file", "Sequence", "Modified sequence", 
                                 "Proteins", "Charge", "m/z", "Mass", "PEP", "Retention time",
                                 "Precursor Intensity", "Masses", "Intensities", "Number of Matches",
                                 "Gene Names", "Evidence ID"))

setnames(pca_urine_sub,
         c("Raw file", "Sequence", "Modified sequence", 
            "m/z", "Mass", "Retention time",
            "Precursor Intensity", "Number of Matches",
            "Gene Names", "Evidence ID"),
         c("raw", "sequence", "ModifiedSequence", "mz", "mass", "RetentionTime",
           "PrecursorIntensity", "NumMatches", "Gene", "EvidenceID"))

sum(is.na(pca_urine_sub$PrecursorIntensity))

runs <- unique(pca_urine_sub$raw)
# 302

pca_urine_sub <- pca_urine_sub[grep("_postdre", raw)]
pca_urine_sub$ModifiedSequence <- reformat_mods(pca_urine_sub$ModifiedSequence)

urine_msms <- lapply(runs[grep("_postdre", runs)], function(x){
  res <- pca_urine_sub[raw == x]
  res <- res[order(PEP)]
  res <- res[!duplicated(ModifiedSequence, Charge)]
  res
})

urine_peptides_intersect <- lapply(urine_msms, function(x){
  peptides <- x$ModifiedSequence
  peptides
})

urine_peptides_intersect <- Reduce(intersect, urine_peptides_intersect)

urine_intersect <- lapply(urine_msms, function(x){
  res <- x
  res <- res[ModifiedSequence %in% urine_peptides_intersect]
  res <- subset(res, select=c("raw", "ModifiedSequence", "RetentionTime"))
  name <- paste0(unique(res$raw), "_RT")
  setnames(res,"RetentionTime", name)
  res <- subset(res, select=c("ModifiedSequence", name))
  res
})
names(urine_intersect) <- unlist(lapply(urine_msms, function (x) unique(x$raw)))

urine_intersect_rt_long <- lapply(1:length(urine_intersect),function(x){
  res <- urine_intersect[[x]]
  run <- names(urine_intersect)[x]
  colnames(res)[2] <- "RT"
  res$run <- run
  res <- spread(res, ModifiedSequence,RT)
  res
})

urine_intersect_rt_long <- do.call('rbind', urine_intersect_rt_long)
rownames(urine_intersect_rt_long) <- urine_intersect_rt_long$run
#colnames(urine_intersect_rt_long) <- urine_intersect_rt_long[1,]
urine_intersect_rt_long$run <- names(urine_intersect)
pca_data <- urine_intersect_rt_long
pca_data$run <- NULL
#urine_intersect_rt_long <- data.frame(sapply(urine_intersect_rt_long,
#                                             function(x) as.numeric(as.character(x))), stringsAsFactors = FALSE)
urine_intersect_rt_long.pca <- prcomp(pca_data)
library(ggbiplot)
ggbiplot(urine_intersect_rt_long.pca, varname.size = 0, var.axes = 0, obs.scale = 1, 
         groups=urine_intersect_rt_long$run) + theme(legend.position = NULL) +
  theme_classic(base_size=16) +
  scale_color_viridis_d(option="D") +
  guides(color=FALSE)

head(urine_intersect_rt_long.pca)


library(reshape)

urine_intersect_rt <- Reduce(function(x,y) merge(x,y, by=c("ModifiedSequence")), urine_intersect)
lnmod <- lm(`170825_postdreredo_4_3_5_2_RT` ~ `170825_postdreredo_4_3_5_1_RT`, urine_intersect_rt)

urine_intersect_runs <- colnames(urine_intersect_rt)[-1]
picks <- sample(urine_intersect_runs, 2)

ggplot(urine_intersect_rt, 
       aes(x = `170908_postdreredo_3_3_16_1_RT`, y = `170825_postdreredo_4_3_4_2_RT`)) +
  geom_point(size=1, alpha=.8) + theme_classic(base_size = 16) + 
  geom_abline(slope=1, color="red") +
  geom_smooth(linetype=0)

library(devtools)



############## Look at Protein count over runs #######################

pca_prot_col <- colnames(pca_urine_protein)
pca_prot_col[-grep("postDRE", pca_prot_col)]
pca_urine_protein <- pca_urine_protein[-grep("REV_", pca_urine_protein$`Protein IDs`)]
pca_urine_protein <- pca_urine_protein[`Q-value` <= 0.01]
pca_urine_protein <- pca_urine_protein[Reverse != "+"]
pca_urine_protein <- pca_urine_protein[`Potential contaminant` != "+"]

pca_urine_protein_count <- subset(pca_urine_protein, select = pca_prot_col[grep("iBAQ postDRE", pca_prot_col)])

# Remove Nans
#pca_urine_protein_count <- na.omit(pca_urine_protein_count)
library(dplyr)
pca_urine_protein_count %>% as.matrix() %>%  rowMedians() %>% head()

pca_urine_protein_count$median <- pca_urine_protein_count %>% as.matrix() %>% rowMedians
head(pca_urine_protein_count)
rowSums(pca_urine_protein_count == 0) %>% head()

pca_urine_protein_count$median <- unlist(lapply(1:nrow(pca_urine_protein_count), function(x){
  int <- as.numeric(as.vector(pca_urine_protein_count[x,]))
  int <- sort(int)
  int <- int[int != 0]
  med <- median(int)
  med
  }))
pca_urine_protein_count$observations <- rowSums(pca_urine_protein_count[,1:151] != 0)
pca_urine_protein_count$MedianIntensity <- pca_urine_protein_count$median/1e10
pca_protein_counts <- subset(pca_urine_protein_count, select=c("observations", "median", "MedianIntensity"))
library(scales)
library(viridis)
pca_protein_counts$logMedian <- log10(pca_protein_counts$median)
ggplot(pca_protein_counts, aes(x = observations, y = median)) + 
  geom_point(aes(color=median), alpha = .6) + scale_y_log10() + 
  scale_color_viridis(option="C", discrete = F) + theme_classic(base_size = 16)
