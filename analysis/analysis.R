library(ggplot2)
library(data.table)
library(VennDiagram)
library(viridisLite)

#setwd("/home/klannie/Desktop/dia_pipeline/results/ankit_HEKqc/")

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

ggplot(fullstat[m_score <= 0.01], aes(x = delta_rt, fill = replicate)) + geom_density(alpha = 0.5) + 
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

dynamic_window <- fullstat[acquisition == "dynamic"]

static_peptide_counts <- as.data.frame(table(static_peps))

ggplot(static_peptide_counts, aes(x = Freq)) + geom_bar()

static_peptides <- unique(static_window$Sequence)


####################################### Load DDA data ##################################################################

dda <- "/home/klannie/Desktop/dia_pipeline/data/dda/ankit_HEKqc/"
mqoutEv <- list.files(dda, pattern="evidence.txt", full.names = T)
mqoutEv <- fread(mqoutEv)
mqoutAllPep <- list.files(dda, pattern="peptides.txt", full.names = T)
mqoutAllPep <- fread(mqoutAllPep)
mqoutPG <- list.files(dda, pattern = "proteinGroups.txt", full.names = T)
mqoutPG <- fread(mqoutPG)

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

ggplot(static_window, aes(x = aggr_prec_Peak_Area, fill = DDA)) + geom_density(alpha = 0.5) +xlim(0, 1e+8)

write.table(static_window, file = "static_window_pyprophet.tsv", row.names = FALSE, sep = "\t")
