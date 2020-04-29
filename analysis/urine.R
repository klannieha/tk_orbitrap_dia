############### urine library #####################

library(data.table)
library(ggplot2)
library(ggExtra)
library(tidyr)
library(stringr)

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")
basefolder <- "D:/projects/pca_urine_spectral_lib/"

# load data

pca_postDREurine.msms <- list.files(paste0(basefolder, "data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og"),
                                    pattern = "msms.txt", full.names = T) %>% fread()

pca_postDREurine.msms$`ModifiedPeptideSequence` <- reformat_mods(pca_postDREurine.msms$`Modified sequence`)
pca_postDREurine.ev <- list.files(paste0(basefolder, "data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og"),
                                  pattern = "evidence.txt", full.names = T) %>% fread()

irt <- list.files(paste0(basefolder, "/data/irt"), pattern="irt_assay", full.names = T) %>% fread()
irt <- subset(irt, select = c("ModifiedPeptideSequence", "NormalizedRetentionTime")) 

irt <- irt[!duplicated(irt)]
irt <- irt[order(NormalizedRetentionTime)]
irt$id <- 1:nrow(irt)
urine_files <- unique(pca_postDREurine.msms$`Raw file`)
pca_postDREurine.msms_sub <- lapply(urine_files, function(x){
  rep <- pca_postDREurine.msms[`Raw file` == x]
  rep <- subset(rep, select=c("Raw file","Sequence", "Modified sequence", "Charge",
                              "Retention time", "PEP", "Precursor Intensity", "Proteins", "Gene Names"))
  setnames(rep, colnames(rep), c("raw","PeptideSequence", "ModifiedPeptideSequence", "Charge",
                                 "RetentionTime", "pep", "Intensity", "Proteins", "Genes"))
  rep$ModifiedPeptideSequence <- reformat_mods(rep$ModifiedPeptideSequence)
  rep
})

names(pca_postDREurine.msms_sub) <- urine_files

############### Basic checking of the data ############

library(readxl)

patient_df <- list.files(paste0(basefolder, "data/clinical"), pattern = "xls", full.names = T)
patient_df <- patient_df %>% read_xls(sheet = "Direct-EPS(T.V.)+EPS-Urines(V.)")

col_id <- grep("Toronto", colnames(patient_df))

patient_df <- subset(patient_df, select = c("patient_id", colnames(patient_df)[col_id]))
patient_df <- patient_df %>% drop_na()
setnames(patient_df, colnames(patient_df), c("PatientID", "DirectEPS", "postDRE_urine"))
patient_df$postDRE_urine <- gsub("X", "S", patient_df$postDRE_urine)
patient_df <- patient_df %>% melt(id = "PatientID", value.name = "RunID")
pid <- patient_df$PatientID %>% unique()

raw <- pca_postDREurine.msms$`Raw file` %>% unique()
urine_id <- str_extract(raw, "S[0-9]*")
raw_id <- data.frame(raw_file = raw, run_id = urine_id, stringsAsFactors = F)
raw_id <- merge(raw_id, patient_df, by.x = "run_id", by.y = "postDRE_urine", all = T)
raw_id$DirectEPS <- NULL
pca_postDREurine.msms <- merge(pca_postDREurine.msms, raw_id, by.x = "Raw file", by.y = "raw_file")

msms_by_patient <- pca_postDREurine.msms[!is.na(PatientID)]

msms_irt_by_patient <- lapply(pid, function(x){
  rep <- msms_by_patient[PatientID == x]
  rep <- subset(rep, select = c("Retention time","ModifiedPeptideSequence", "Charge", "PEP", "run_id", "PatientID"))
  rep <- rep[order(PEP)]
  rep <- rep[!duplicated(ModifiedPeptideSequence)]
  rep$PEP <- NULL
  rep <- merge(irt, rep, by = c("ModifiedPeptideSequence"), all.x = T)
  rep$run_id <- NULL
  rep$Charge <- NULL
  rep$PrecursorMz <- NULL
  rep$iRT <- NULL
  setnames(rep, "Retention time", "postDRE_urine")
  rep
})


names(msms_irt_by_patient) <- pid


ggplot(msms_irt, aes(x = id, y = `Retention time`, group = id)) +
  geom_jitter(aes(color = `Raw file`), show.legend = F, width = 0.3, size = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.8) + theme_bw() + xlab("iRT_id")


urine.iRT <- lapply(urine_files, function(x){
  rep <- pca_postDREurine.msms_sub[[x]]
  rep <- rep[order(pep)]
  rep <- rep[!duplicated(rep, by = c("ModifiedPeptideSequence"))]
  rep <- subset(rep, select=c("ModifiedPeptideSequence", "RetentionTime", "raw", "Intensity"))
  rep <- merge(rep, irt, by = "ModifiedPeptideSequence", by.y = all)
  rep[, run:= x]
  rep
})


names(urine.iRT) <- urine_files

urine.iRT <- do.call('rbind', urine.iRT )
urine.iRT$iRT <- factor(urine.iRT$ModifiedPeptideSequence,levels = irt$ModifiedPeptideSequence)
urine.iRT$cohort <- "postDRE_urine"
ggplot(urine.iRT, aes(x = iRT, y = RetentionTime, color = run, group = iRT)) + 
  theme_bw(base_size = 16) +
  geom_jitter(show.legend = F, width = 0.3, size=0.2, alpha=0.8) +
  geom_boxplot(show.legend = F, width = 0.5, alpha = 0.8) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

iRT_all <- rbind(pca_eps.iRT, urine.iRT)
iRT_all <- iRT_all[run != "200117_postDREurine_R002_S030_1"]

ggplot(iRT_all, aes(x = RetentionTime, y = NormalizedRetentionTime, color = run)) + 
  geom_point(show.legend = F) +
  geom_line(show.legend = F) + facet_wrap(~cohort, nrow = 2) + theme_bw(base_size = 16) +
  scale_color_viridis_d(option = "plasma", end = 0.5)



#### Running checks between directEPS and urine ####

eps_raw_files <- pca_eps.msms$`Raw file` %>% unique()
eps_run_id <- str_extract(eps_raw_files, "E[0-9]*")

eps_id <- data.frame(raw_file = eps_raw_files, run_id = eps_run_id)
eps_id <- merge(eps_id, patient_df, by.x = "run_id", by.y = "DirectEPS", all = T)
eps_id$postDRE_urine <- NULL
eps_msms <- merge(pca_eps.msms, eps_id, by.x = "Raw file", by.y = "raw_file")
eps_msms$ModifiedPeptideSequence <- reformat_mods(eps_msms$`Modified sequence`)
eps_msms_irt <- lapply(pid, function(x){
  rep <- eps_msms[PatientID == x]
  rep <- subset(rep, select = c("Retention time","ModifiedPeptideSequence", "Charge", "PEP", "run_id", "PatientID"))
  rep <- rep[order(PEP)]
  rep <- rep[!duplicated(ModifiedPeptideSequence)]
  rep$PEP <- NULL
  rep <- merge(irt, rep, by = c("ModifiedPeptideSequence"), all.x = T)
  rep$Charge <- NULL
  rep$run_id <- NULL
  rep$PrecursorMz <- NULL
  rep$iRT <- NULL
  setnames(rep, "Retention time", "DirectEPS")
  rep
})

names(eps_msms_irt) <- pid

all_msms_irt <- lapply(pid, function(x){
  urine <- msms_irt_by_patient[[x]]
  eps <- eps_msms_irt[[x]]
  all <- merge(urine, eps, by = c("ModifiedPeptideSequence", "PatientID"))
  all <- drop_na(all)
  all
})

all_combined <- do.call('rbind', all_msms_irt)
ggplot(all_combined, aes(x = DirectEPS, y = postDRE_urine, color = PatientID)) + 
  geom_line(size = 1) +
  scale_color_viridis_d(option = "plasma") +
  theme_bw()


cor <- lapply(all_msms_irt, function(x) cor(x$DirectEPS, x$postDRE_urine)) %>% unlist()
cor
pid_cor <- data.table(PatientID = pid, correlation = cor)
pid_cor <- drop_na(pid_cor)
color <- viridis(10, option = "inferno")
ggplot(pid_cor, aes(x = PatientID, y = correlation)) + 
  geom_bar(stat ="identity", position = position_dodge(width = 0.9), width = 0.6, fill = color[2]) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle=30, vjust = 1)) + 
  coord_cartesian(ylim = c(0.92, 1.0))
############## Create library ##############

irt_counts <- lapply(urine_files, function(x){
  rep <- pca_postDREurine.msms[`Raw file` == x]
  rep <- rep[order(PEP)]
  rep <- subset(rep, select = c("Retention time","ModifiedPeptideSequence"))
  rep <- rep[!duplicated(rep, by = c("ModifiedPeptideSequence"))]
  rep <- merge(irt, rep, by = c("ModifiedPeptideSequence"))
  rep$PrecursorMz <- NULL
  setnames(rep, "Retention time", x)
  rep
})

irt_counts <- Reduce(function(x,y) merge(x,y, by = c("ModifiedPeptideSequence", "NormalizedRetentionTime"), all = T), irt_counts)

stats <- lapply(1:nrow(irt_counts), function(x){
  m <- irt_counts[x, 2:ncol(irt_counts)] %>% as.numeric() %>% mean(na.rm=T)
  sd <- irt_counts[x, 2:ncol(irt_counts)] %>% as.numeric() %>% sd(na.rm=T)
  na <- irt_counts[x, 2:ncol(irt_counts)] %>% is.na() %>% sum()
  c(m, sd, na)
})

irt_counts <- do.call('rbind', stats)
irt_counts <- as.data.table(irt_counts)
colnames(irt_counts) <- c("Mean", "SD", "Missing")
irt_counts$ModifiedPeptideSequence <- irt$ModifiedPeptideSequence
irt_counts$iRT <- irt$NormalizedRetentionTime
irt_counts <- irt_counts[order(iRT)]

pep <- irt_counts$ModifiedPeptideSequence
irt_counts$ModifiedPeptideSequence <- factor(irt_counts$ModifiedPeptideSequence, levels = pep)
urine.stats <- irt_counts
urine.stats$iRT <- NULL
urine.stats$cohort <- "postDRE_urine"
setnames(urine.stats, "Mean", "mean")
setnames(eps.stats, "Peptide", "ModifiedPeptideSequence")
eps.stats$iRT <- NULL

all_stats <- rbind(urine.stats, eps.stats)

ggplot(all_stats, aes(x = ModifiedPeptideSequence, y = Missing, fill = cohort)) + 
  geom_bar(stat="identity", position = position_dodge(width = 1)) + theme_bw(base_size = 16) +
  scale_fill_viridis_d(option = "plasma", end = 0.5, direction = -1) +
  theme(axis.text.x = element_text(angle=30, vjust = 1)) + xlab("iRT Peptides") + ylab("Number of Sample Missing")
#SD between 2 to 4 minutes

############# Perform pairwise ##########################

urine.seq <- lapply(urine_files, function(x){
  seq <- pca_postDREurine.msms[`Raw file` == x]$ModifiedPeptideSequence %>% unique()
  seq
})

pairwise <- matrix(nrow=length(urine.seq), ncol = length(urine.seq))

for (i in 1:length(urine_files)) {
  r <- urine.seq[[i]]
  for(j in 1:length(urine_files)){
    s <- urine.seq[[j]]
    pair <- intersect(r, s) %>% length()
    pairwise[i,j] <- pair
  }
}

library(viridisLite)
library(pheatmap)
runs <- gsub("(.*)postDREurine_", "", urine_files)
row.names(pairwise) <- runs
colnames(pairwise) <- runs
pheatmap(pairwise, color = viridis(20))

# Run R100_S083 seems to have very low pairwise match
grep("(.*)R100_S083(.*)", raw) # 198

length(urine.seq[[198]])
# 1004

# run S045 has the best pairwise
grep("(.*)_S045", raw) # 4

ref_count <- data.frame(S045_count = pairwise[4,], run = row.names(pairwise)) %>% as.data.table()
color <- viridis(2, end = 0.5, option = "plasma")
ggplot(ref_count[run != "R052_S045_1"], aes(x = as.integer(run), y = S045_count)) + 
  geom_bar(stat="identity", position = "dodge", width = 0.5, fill = color[1]) + theme_bw(base_size = 16) +
  xlab("Sample") + ylab("Peptide Overlap with S045") + ylim(0, 12000)

# Total number of peptides overlap
AllPeptides.urine <- pca_postDREurine.msms$`Modified sequence` %>% reformat_mods() %>% unique()
length(AllPeptides.urine)
# 57,494 unique peptides across ALL

names(pca_postDREurine.msms_sub) <- urine_files
peptide_count.urinee <- lapply(urine_files, function(x){
  rep <- pca_postDREurine.msms_sub[[x]]
  precursor <- rep[!duplicated(ModifiedPeptideSequence, Charge)]
  precursor <- precursor %>% nrow()
  peptide <- rep$PeptideSequence %>% unique() %>% length()
  protein <- rep$Proteins %>% unique() %>% length()
  count <- data.frame(raw=x, precursor=precursor, peptide=peptide, protein=protein)
  count
})
peptide_count.urinee <- do.call('rbind', peptide_count.urinee)

ggplot(peptide_count.urinee, aes(x = as.numeric(raw), y = peptide)) + 
  geom_bar(width = 0.5, stat = "identity") + theme_bw() + 
  xlab("Run") + ylab("Precursor Count")

###################### Use run S045 as reference ######################
msms_ref <- pca_postDREurine.msms[`Raw file` == raw[4]]
msms_ref <- subset(msms_ref, select = c("Raw file", "Sequence", "ModifiedPeptideSequence", "Proteins", 
                                        "Gene Names", "Charge", "m/z", "Mass", "Retention time",
                                        "PEP", "Intensities", "Masses", "Evidence ID"))
ev_ref <- subset(pca_postDREurine.ev, select = c("id", "Intensity"))

ms_ref <- merge(msms_ref, ev_ref, by.x = "Evidence ID", by.y = "id")
ms_ref <- ms_ref[order(PEP)]
ms_ref <- ms_ref[!duplicated(ms_ref, by = c("ModifiedPeptideSequence", "Charge"))]

ms_ref$`Evidence ID` %>% unique() %>% length()
# 19100 precursors
ms_ref$Sequence %>% unique() %>% length()
# 15450 peptides
ms_ref$`Gene Names` %>% unique() %>% length()
# 2715
ms_ref$Proteins %>% unique() %>% length()
# 3623 protein groups

ref_irt <- merge(irt, ms_ref, by = "ModifiedPeptideSequence")
ref_irt <- ref_irt[order(iRT)]
ggplot(ref_irt, aes(y = iRT, x = `Retention time`)) + geom_point() + geom_line() + geom_smooth(method = "lm")
setnames(ref_irt, "iRT", "NormalizedRetentionTime")
ref_irt_mod <- lm(NormalizedRetentionTime ~ `Retention time`, ref_irt)
outlier <- abs(ref_irt_mod$residuals) > 10
ref_irt <- ref_irt[abs(ref_irt_mod$residuals) <= 10]
ggplot(ref_irt, aes(y = NormalizedRetentionTime, x = `Retention time`)) + geom_point() + 
  geom_line() + geom_smooth(method = "lm")
ref_irt_mod <- lm(NormalizedRetentionTime ~ `Retention time`, ref_irt)

# Use this model to align the rest of this run
ms_ref$NormalizedRetentionTime <- predict(ref_irt_mod, newdata = ms_ref)
ggplot(ms_ref, aes(x = `Retention time`, y = `NormalizedRetentionTime`)) + geom_point() + 
  geom_line() + geom_smooth(method = "lm")

setnames(ms_ref, c("Evidence ID","Gene Names", "Sequence"), c("id","UniprotID", "PeptideSequence"))
setnames(ms_ref, "Sequence", "PeptideSequence")

write.table(ms_ref, file = paste0(basefolder, "data/irt/postDRE_urine_Alignment_hpiRT.tsv"),
            sep = "\t", quote = F, row.names = F)



