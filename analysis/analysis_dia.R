library(ggplot2)
library(data.table)
library(VennDiagram)
library(viridisLite)

#setwd("/home/klannie/Desktop/dia_pipeline/results/ankit_HEKqc/")
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

files <- list.files(path = '.', pattern = "newParams_pyprophet_export.tsv")
stats <- lapply(files, fread)

################################## Filtering and Data Processing ################################

filtered <- lapply(1:length(files), function(x){
  res <- stats[[x]]
  res <- res[decoy == 0 & peak_group_rank ==1]
  acq <- strsplit(files[x], "_")[[1]][4]
  rep <- strsplit(files[x], "_")[[1]][5]
  res[, acquisition := acq]
  res[, replicate := rep]
  res
})

fullstat <- do.call('rbind', filtered)
head(fullstat)
fullstat <- fullstat[m_score <= 0.01]

ggplot(fullstat, aes(x = delta_rt, fill = replicate)) + geom_density(alpha = 0.5) + 
  facet_wrap(~acquisition) + theme_classic(base_size = 16)
ggsave("delta_rt.png")
peptides <- lapply(filtered, function(x){
  res <- x[m_score_peptide_run_specific <= 0.01]
  p <- unique(res$Sequence)
  p
})

allPep <- Reduce(intersect, peptides)

static_id <- grep("(.*)_static_(.*)", files)
static <- lapply(static_id, function(x){
  res <- filtered[[x]]
  res <- res[m_score_peptide_run_specific <= 0.01]
  res
})

static_peps <- lapply(static, function(x){
  res <- x[m_score_peptide_run_specific <= 0.01]
  p <- unique(res$Sequence)
})

static_peps <- unlist(static_peps)
static_window <- fullstat[acquisition == "static"]
static_peptides <- Reduce(intersect, static_peps)
length(static_peptides)
static_peptides <- static_window$Sequence[duplicated(static_window[m_score_peptide_run_specific <= 0.01]$Sequence)]
length(unique(static_window$Sequence))

dynamic_window <- fullstat[acquisition == "dynamic" & m_score <= 0.01]
static_window <- fullstat[acquisition == "static" & m_score <= 0.01]

static_peptide_counts <- as.data.frame(table(static_peps))

ggplot(static_peptide_counts, aes(x = Freq)) + geom_bar()

static_peptides <- unique(static_window$Sequence)


####################################### Load DDA data ##################################################################

mqoutEv <- list.files('.', pattern="evidence.txt", full.names = T)
mqoutEv <- fread(mqoutEv)
mqoutAllPep <- list.files('.', pattern="peptides.txt", full.names = T)
mqoutAllPep <- fread(mqoutAllPep)
#mqoutPG <- list.files(dda, pattern = "proteinGroups.txt", full.names = T)
#mqoutPG <- fread(mqoutPG)

patt <- ".*(AS_1| AK_A11| AK_A18| AS_QC| AS_2).*"
col <- colnames(mqoutAllPep)
col <- col[-grep(patt, col)]
qc <- subset(mqoutAllPep, select =col )
colnames(qc)
nrow(qc)
nrow(qc[!is.na(Intensity)])
nrow(qc[!is.na(`LFQ intensity AK_QC`)])
sum(qc$`LFQ intensity AK_QC`==0)
#34638
sum(qc$`Intensity AK_QC` == 0)
#34638
sum(is.na(qc$`Experiment AK_QC`))
sum(is.na(qc$`Identification type AK_QC`))
#33585
sum(qc$`Identification type AK_QC`== "")
#33585
qcAK <- qc[!is.na(`Experiment AK_QC`)]
qcAK <- qcAK[PEP <= 0.01 & `Intensity AK_QC` != 0]
length(unique((qcAK$Sequence)))
#32526
# filtered 26867
qcAK <- qcAK[order(PEP)]

head(unique(qcAK$Sequence))

idx <- which(unique(static_peps) %in% unique(qcAK$Sequence))


#################### Load DIA lib #################################
lib <- fread("/home/klannie/Desktop/dia_pipeline/data/libraries/HEK_qc_mqout_SCX_10only_assaylib.tsv")

lib <- lib[!duplicated(PeptideSequence) & Decoy == 0]

head(qcAK[-(qcAK$Sequence %in% lib$PeptideSequence)])

########################################3 Plotting ##########################################3
ggplot(static_window, aes(x = sn_ratio)) + geom_density()
length(intersect(allPep, qcAK$Sequence))
# 9627 # 9010
length(intersect(static_peptides, unique(qcAK$Sequence)))
#11862
png("peptideOverlap_static_window.png")
vp <- venn.diagram(list(unique(qcAK$Sequence), static_peptides), filename = NULL,
                   category.names = c("DDA", "DIA"))
grid.draw(vp)
dev.off()

png("static_dynamic.png")
vp <- venn.diagram(list(unique(static_window$Sequence), unique(dynamic_window$Sequence)),
                   filename=NULL,
                   category.names = c("static", "dynamic"),main.col = NULL,
                   col=c(colors[1], colors[15]), fill = c(colors[1], colors[15]), alpha=0.6)
grid.draw(vp)
dev.off()

png("peptideOverlap.png")
vp <- venn.diagram(list(unique(qcAK$Sequence), allPep), filename = NULL,
             category.names = c("DDA", "DIA"))
grid.draw(vp)
dev.off()

png("libDDAOverlap.png")
vp <- venn.diagram(list(unique(qcAK$Sequence), lib$PeptideSequence), filename = NULL,
                   category.names = c("DDA", "Lib"))
grid.draw(vp)
dev.off()

unique_dia_pep <- static_window[-which(Sequence %in% unique(qcAK$Sequence))]
static_window$DDA <- static_window$Sequence %in% unique(qcAK$Sequence)

####### Find high med low scoring peptides #######
fullstat <- fullstat[order(m_score)]
head(fullstat$Sequence)
head(fullstat[which(unique(fullstat$Sequence) %in% unique(qcAK$Sequence))])
# Compare quality with the DDA runs

####### Re make some graphs ######################
tp_counts <- lapply(filtered, function(x){
  res <- x
  res <- res[m_score <= 0.01]
  prec_count <- nrow(res[!duplicated(FullPeptideName, Charge)])
  res <- res[m_score_peptide_run_specific <= 0.01]
  pep_count <- length(unique(res$Sequence))
  res <- res[m_score_protein_run_specific <= 0.01]
  pro_count <- length(unique(res$ProteinName))
  c(prec_count, pep_count, pro_count, unique(res$acquisition), unique(res$replicate))
})

tp_counts <- do.call('rbind', tp_counts)
tp_counts <- as.data.table(tp_counts)
colnames(tp_counts) <- c("Precursor", "Peptide", "Protein", "Acquisition", "Replicate")
tp_long <- melt(tp_counts, id.vars = c("Acquisition", "Replicate"))
typeof(tp_long$value)
tp_long$value <- as.numeric(tp_long$value)
ggplot(tp_long, aes(x = Acquisition, y = value, fill=Replicate, width=0.8)) + 
  geom_bar(stat="identity", position = position_dodge(1)) +
  facet_wrap(~variable, scales = "free") +
  theme_classic(base_size = 16)# + scale_y_continuous(limits = c(0, 25000))

# Library Generation
lib <- fread("../library/msms.txt")
irt <- fread("../library/irt/CiRT_ALL.tsv")
lib <- lib[order(PEP)]
lib <- lib[PEP <= 0.01]
lib <- subset(lib, select=c("Raw file", "Modified sequence", "Charge", "Retention time"))

colnames(lib) <- c('raw', 'sequence', 'charge', 'rt')
lib$sequence <- reformat_mods(lib$sequence)
raw <- unique(lib$raw)
lib <- lapply(raw, function(x){
  sublib <- lib[raw == x]
  sublib<- sublib[!duplicated(sequence, charge)]
  sublib
})
colors <- viridis(17, option = "C")
lib <- do.call('rbind', lib)
irt <- subset(irt, select=c("NormalizedRetentionTime", "ModifiedPeptideSequence", "PrecursorCharge"))
colnames(irt) <- c("irt", "sequence", "charge")
msms_irt <- merge(irt, lib, by = c("sequence", "charge"))
ggplot(msms_irt[raw == raw[4]], aes(x = rt, y = irt)) + geom_point(show.legend = F, color = colors[4]) +
  theme_classic(base_size = 16) + xlab("Retention time(min)") + ylab("iRT") +
  geom_smooth(method = "lm", se = F, show.legend = F, color=colors[4])

mqout <- fread("../../dia_pipline/library/mqout.tsv")
length(unique(mqout$PeptideSequence))
length(unique(mqout$transition_group_id))
length(unique(mqout$ProteinName))
