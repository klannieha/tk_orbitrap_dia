#!/bin/usr/ R
# module load gcc/7.3.0
# module load r/3.6.0
# Better to request for enough memory !!


library(data.table)
library(ggplot2)
library(tidyr)
#library(MSstats)
library(SWATH2stats)
library(aLFQ)
library(ggpubr)

GeneFormat <- function(x){
  if(grepl("Subgroup", x)){
    p <- unlist(strsplit(x, "\\|"))
    p <- p[grepl("^[A-Z][0-9](.*)", p)]
    p <- gsub("-.*", "", p)
    p <- unique(p)
    p <- sort(p)
    p <- paste0(p, collapse=";")
  } else{
    p <- strsplit(x, ";")[[1]]
    p <- gsub("-.*", "", p)
    p <- unique(p)
    p <- sort(p)
    p <- paste0(p, collapse=";")
  }
  return(p)
}

reduce_OpenSWATH_output <- function(data, column.names=NULL){
  if(is.null(column.names)){
    column.names <- c('ProteinName', 'FullPeptideName', 'Sequence', 'Charge', 'aggr_Fragment_Annotation', 'aggr_Peak_Area', 'filename', 'm_score', 'decoy', "Intensity", "RT", "run_id", "transition_group_id")
  }
  if(length(column.names) > length(column.names[column.names %in% colnames(data)])){
    col.names.missing <- column.names[!column.names %in% colnames(data)]
    warning("These columns are missing from the data:", paste(unlist(col.names.missing), collapse=", "))

  }
  # Keep only required columns for MSStats and mapDIA
  if(length(column.names) == length(column.names[column.names %in% colnames(data)])){
    data.filtered <- data[,..column.names]
    return(data.filtered)
  }
}

#data <- '/project/6002011/annieha/pca_urine_spectral_lib/data/openswath/200903_HEK_DIA'

setwd("D:/projects/pca_urine_spectral_lib/data/openswath/200929_Mstern_DIA/")

data <- list.files('.', pattern="pyprophet_export.tsv")

stats <- lapply(data, fread)

#files <- list.files(data, pattern="_0914lib_all_merged_pyprophet_export.tsv")

#stats <- lapply(files, function(x)  fread(x, sep = "\t"))

# load experiment conditions with corresponding filenames
annotation <- fread("annotation.txt")
#setnames(annotation, "Run", "filename")

lib.path <- '/project/6002011/annieha/pca_urine_spectral_lib/data/library/HEK/200914_tpp_HEK_assaylib_osw_target_decoy.tsv'
lib <- fread(lib.path)

# first fix the protein IDs by matching to the library protein IDs
#proteinIds <- lib$ProteinId %>% unique()

#dep_proteinIds <- gsub(";", "", proteinIds)

#lib.protein <- data.frame(proteinIds, dep_proteinIds, stringsAsFactors=F)

# Analysis: basic peptide and protein counts:

fixedWindows <- list.files("D:/projects/pca_urine_spectral_lib/results/2020Sept_HEK_dia/", pattern = "methodopt", full.names = T)
fixed <- lapply(fixedWindows, fread)

fixed <- lapply(fixed, function(x){
  rep <- x
  rep <- rep[decoy == 0]
  rep <- rep[m_score <= 0.01]
  rep <- rep[peak_group_rank == 1]
  rep <- rep[order(Intensity)]
  rep <- rep[!duplicated(rep, by = c("FullPeptideName", "Charge", "transition_group_id", "filename"))]
  #  rep <- merge(rep, annotation, by = "filename", all.x = T)
  rep})

fixed <- do.call('rbind', fixed)
fixed$GENE <- unlist(lapply(fixed$ProteinName, function(x) GeneFormat(x)))

cols <- c("UNIPROT", "SYMBOL")

genes <- select(org.Hs.eg.db, keys=unique(fixed$GENE), keytype = "UNIPROT", columns = cols)
fixed <- merge(fixed, genes, by.x="GENE", by.y="UNIPROT", all.x=T)

count <- lapply(annotation$filename, function(x){
  rep <- fixed[filename == x]
  rep <- rep[!duplicated(rep, by = c("FullPeptideName", "Charge", "transition_group_id"))]
  precursor <- nrow(rep)
  peptide <- rep$Sequence %>% unique() %>% length()
  rep <- rep[!grepl("DECOY", rep$ProteinName)]
  proteinGroups <- rep$ProteinName %>% unique() %>% length()
  genes <- rep[!grepl("DECOY|;", rep$GENE)]$SYMBOL %>% unique() %>% length()
  c <- c(precursor, peptide, proteinGroups, genes)
  c})

count <- do.call('rbind', count) %>% as.data.table()
colnames(count) <- c("Precursor", "Peptide", "Protein", "Genes")
count <- cbind(count, annotation)
count <- count[9:12]

count <- melt(count, id.vars=c("Condition", "BioReplicate", "filename"), value.name="count")
count$BioReplicate <- factor(count$BioReplicate, levels=c(1,2))



stats <- lapply(stats, function(x){
  rep <- x
  rep <- rep[decoy == 0]
  rep <- rep[m_score <= 0.01]
  rep <- rep[peak_group_rank == 1]
  rep <- rep[order(Intensity)]
  rep <- rep[!duplicated(rep, by = c("FullPeptideName", "Charge", "transition_group_id", "filename"))]
#  rep <- merge(rep, annotation, by = "filename", all.x = T)
  rep})

stats <- do.call('rbind', stats)
stats$GENE <- unlist(lapply(stats$ProteinName, function(x) GeneFormat(x)))

cols <- c("UNIPROT", "SYMBOL")

genes <- select(org.Hs.eg.db, keys=unique(stats$GENE), keytype = "UNIPROT", columns = cols)
stats <- merge(stats, genes, by.x="GENE", by.y="UNIPROT", all.x=T)

count <- lapply(annotation$filename, function(x){
  rep <- stats[filename == x]
  rep <- rep[!duplicated(rep, by = c("FullPeptideName", "Charge", "transition_group_id"))]
  precursor <- nrow(rep)
  peptide <- rep$Sequence %>% unique() %>% length()
  rep <- rep[!grepl("DECOY", rep$ProteinName)]
  proteinGroups <- rep$ProteinName %>% unique() %>% length()
  genes <- rep[!grepl("DECOY|;", rep$GENE)]$SYMBOL %>% unique() %>% length()
  c <- c(precursor, peptide, proteinGroups, genes)
  c})

count <- do.call('rbind', count) %>% as.data.table()
colnames(count) <- c("Precursor", "Peptide", "Protein", "Genes")
count <- cbind(count, annotation)
count <- count[1:8]

count <- melt(count, id.vars=c("Condition", "BioReplicate", "filename"), value.name="count")
count$BioReplicate <- factor(count$BioReplicate, levels=c(1,2))



g <- ggplot(count[variable == "Genes"], aes(x = Condition, y = count, fill = BioReplicate)) +
geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.6) +
geom_text(aes(label = count), position=position_dodge(width=0.8), vjust= -0.2, size=5) +
#facet_wrap(.~variable, scales="free_y")  
theme_light(base_size=20) + xlab("") + ylab("Genes") + guides(fill=guide_legend(title = "Replicate")) +
scale_fill_brewer(palette = "Pastel1")

g
png('/project/6002011/annieha/pca_urine_spectral_lib/results/2020Sept_HEKdia/count.png',
width=800)
print(g)
dev.off()


# Start with protein inference and quantification

stats$UniprotID <- stats$ProteinName
stats$ProteinName <- stats$GENE
stats$GENE <- NULL

stats <- merge(stats, annotation, by.x = "filename", by.y = "Filename")

Allstats <- reduce_OpenSWATH_output(stats) %>% as.data.frame(stringsAsFactors=F)
Allstats$FullPeptideName <- gsub("\\.", "", Allstats$FullPeptideName)
setnames(annotation, "filename", "Filename")
annotation$Run <- 1:nrow(annotation)
annotation <- annotation %>% as.data.frame(stringsAsFactors=F)
Allstats <- sample_annotation(Allstats, annotation)
Allstats <- Allstats[Allstats$filename != "/project/6002011/annieha/pca_urine_spectral_lib/data/dia/20200923_MStern_DIA_mzML/20200925_31Variable_03.mzML.gz",]

# reformat for aLFQ using SWATH2stats functions
disagg <- disaggregate(Allstats)
alfq.in <- convert4aLFQ(disagg, check_transitions = FALSE)

# run aFLQ for peptide inference and CV
peps <- PeptideInference(alfq.in, transition_topx = 6, transition_strictness = "loose", 
                         transition_summary = "sum", consensus_transitions = FALSE)

wide <- dcast(peps, peptide_sequence~run_id, drop = FALSE, fill = NaN, sum, value.var = "peptide_intensity")

# run aLFQ for protein inference

prots <- ProteinInference(alfq.in, peptide_method = "top", peptide_topx = 1, peptide_strictness = "strict",peptide_summary = "sum",
                          transition_topx = 5, transition_strictness = "loose", transition_summary = "sum",
                          fasta = NA, model = NA, combine_precursors = FALSE, combine_peptide_sequences = TRUE, consensus_proteins = FALSE, consensus_peptides = FALSE, consensus_transitions = FALSE)

wide <- dcast(prots, protein_id~run_id, drop=FALSE, fill=NaN, value.var="response")

# normalize
# todo: probably better to normalize in log space and revert to linear for CV
medNorm <- function(x){
  med <- apply(x, 2, median, na.rm=T)
  av.med <- mean(med)
  norm.factors.med <- av.med/med 
  y <- as.data.frame(t(t(x) * norm.factors.med))
  row.names(y) <- row.names(x)
  return(y)
}

# Separate experiments for CV calculation 

wide.exp <- lapply(unique(annotation$Condition), function(x){
   c <- colnames(wide)
   exp <- c[grepl(x, c)]
   rep <- subset(wide, select=c("protein_id", exp))
   rep})


# log transform

wide.log <- log(wide[, -1])
rownames(wide.log) <- wide$peptide_sequence
rownames(wide.log) <- wide$protein_id

# normalize in log space


wide.lognorm <- medNorm(wide.log)
wide.lognorm.median.abund <- apply(wide.lognorm, 1, median, na.rm=T) %>% as.data.table()# get normalized medians

rownames(wide.lognorm.median.abund) <- wide$protein_id

write.table(wide.lognorm.median.abund, file="D:/projects/pca_urine_spectral_lib/results/200929_Mstern_DIA/protein_median_lognorm.tsv", sep="\t", row.names = T,quote=F)

# check normalization

Alllog <- melt(wide.log, value.name="Log(Intensity)",
               variable.name="Method")

g <- ggplot(Alllog, aes(x = Method, y = `Log(Intensity)`)) +
geom_boxplot(width=0.6) + theme_light(base_size=14) +
theme(axis.text.x=element_text(angle=90))

png("Logintensity.png")
print(g)
dev.off()

Alllognorm <- melt(wide.lognorm, value.name="Log(Intensity)",
                   variable.name="Method")

g <- ggplot(Alllognorm, aes(x = Method, y = `Log(Intensity)`)) +
geom_boxplot(width=0.6) + theme_light(base_size=14) +
theme(axis.text.x=element_text(angle=90))
g
png("Normalized_Logintensity.png")
print(g)
dev.off()

# check protein abundance across experiments
exp <- colnames(wide.lognorm.median.abund)

wide.lognorm.allmedian.abund <- subset(wide.lognorm.median.abund,
select=c("14mz_steppedNCE_1_1", "16mz_staggered_1_3", "30Variable_1_5",
"31Variable_1_7"))

g <- ggplot(wide.lognorm.allmedian.abund, aes(x = `16mz_staggered_1_3`, y =
`14mz_steppedNCE_1_1`)) + geom_point() + geom_smooth(method="lm", se=F, formula = y ~ x, color="lightgrey", linetype="dashed") + stat_cor(label.y = 20) + ylab("Median log(Intensity) - stepped NCE") + xlab("Median log(Intensity) - 16mz Staggered") + theme_light(base_size=14)

png("Median_logInt_NCE_Staggered.png")
print(g)
dev.off()


# convert back to linear space for CV
wide.norm <- exp(wide.lognorm)

#wide.norm <- exp(wide.log)

# CV: sd / mean
getCV <- function(x){
  CV = apply(x, 1, function(y) 100*sd(y, na.rm=T)/mean(y, na.rm=T))
  mean = apply(x, 1, function(y) mean(y, na.rm=T))
  x$CV <- CV
  x$mean <- mean
  return(x)
}

#wide.cv <- getCV(wide[, -1]) # CV prior normalization
#wide.norm.cv <- getCV(wide.norm) # CV from linear trans
#wide.log.cv <- getCV(wide.lognorm) # CV from log tras

# CV prior normalization
wide.cv <- lapply(wide.exp, function(x){
  rep <- getCV(x[,-1])
  rownames(rep) <- x$protein_id
  rep})


# CV from linear trans after normalization
wide.norm <- lapply(unique(annotation$Condition), function(x){
   c <- colnames(wide.norm)
   exp <- c[grepl(x, c)]
   rep <- subset(wide.norm, select= exp)
   rep})

wide.norm.cv <- lapply(wide.norm, function(x) getCV(x))

NormCV <- lapply(1:4, function(x){
  exp <- unique(annotation$Condition)[x]
  rep <- wide.norm.cv[[x]]
  rep <- subset(rep, select=c("CV", "mean"))
  rep$Method <- exp
  rep})

NormCV <- do.call('rbind', NormCV) %>% as.data.table()
library(viridisLite)
library(tidyverse)

colours <- viridis(n=10, option="C", end=0.8)
g <- ggplot(NormCV, aes(x = Method, y = CV)) + geom_jitter(alpha=0.2, color=colours[1], width=0.1) + geom_boxplot(width=0.5, alpha=0.5) +
theme_light(base_size=14) + ylab("CV (%)")
g
png("NormCV_boxplot.png")
print(g)
dev.off()


g <- ggplot(NormCV, aes(y = log(mean), x = CV, color=Method)) + geom_point(alpha=0.5) + scale_color_viridis_d(option="C", end=0.8) 
g
png("CV_meanInt.png")
print(g)
dev.off()

NormCV_summary <- NormCV %>% group_by(Method) %>% summarize(mean=mean(CV, na.rm=T), median=median(CV, na.rm=T)) %>% as.data.table()
library(grid)
library(gridExtra)
grid.newpage()
grid.table(NormCV_summary)

ggplot(wide.norm[[4]], aes(x = log(`variable_1_10`), y = log(`variable_2_11`))) + geom_point() +
   theme_bw(base_size=14) + stat_cor(label.y = 20)

# CV from log trans normalization
wide.lognorm.exp <- lapply(unique(annotation$Condition), function(x){
   c <- colnames(wide.lognorm)
   exp <- c[grepl(x, c)]
   rep <- subset(wide.lognorm, select= exp)
   rep})



wide.lognorm.cv <- lapply(wide.lognorm.exp, function(x) getCV(x))

logNormCV <- lapply(1:4, function(x){
  exp <- unique(annotation$Condition)[x]
  rep <- wide.lognorm.cv[[x]]
  rep <- subset(rep, select=c("CV", "mean"))
  rep$Method <- exp
  rep})

logNormCV <- do.call('rbind', logNormCV) %>% as.data.table()

g <- ggplot(logNormCV, aes(x = Method, y = CV)) + geom_jitter(alpha=0.2, color=colours[1], width=0.1) + geom_boxplot(width=0.5, alpha=0.5) +
theme_light(base_size=14) + ylab("CV (%)")

png("logNormCV_boxplot.png")
print(g)
dev.off()



###################### check detection rate/ window #################

DIA_short <- lapply(stats, function(x){
  rep <- x
  rep <- merge(rep, annotation, by = "filename")
  rep <- subset(rep, select=c("RT", "mz", "filename", "Condition", "BioReplicate"))
  rep
})

# load windows

steppedNCE <- data.frame(start = seq(400, 1000-15, by=15), end = seq(415, 1000, by=15), stringsAsFactors = F )
steppedNCE$Centre <- steppedNCE$start + 15/2
Variable_30 <- fread("D:/documents/project_protocols/DIA/30_windows_constantIons.tsv")
Staggered <- read_xlsx("D:/documents/project_protocols/DIA/16mz_staggered.xlsx", sheet = 2)
Variable_31 <- fread("D:/documents/project_protocols/DIA/DIA_31Var.txt")

test <- DIA_short[[1]]

test$RT <- test$RT/60
test <- binColumn(test, "RT", 1)
test$RT <- round(test$RT)
test <- binColumn(test, "mz", 1)
set.seed(123)
n <- sample(1:500, 100)

library(reshape2)

t1 <- dcast(test, mz ~ RT, fun.aggregate = length)
t1.d <- melt(t1, id.vars = "mz", variable.name = "RT", value.name = "count")
str(t1)
t1.m <- as.matrix(t1[,-1])
row.names(t1.m) <- t1$mz %>% round()
colnames(t1.m) <- colnames(t1) %>% round()
library(pheatmap)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
h <- heatmap(t1.m, Rowv = NA,Colv = NA, xlab = "RT", ylab="mz", col = coul, cexRow = 1, cexCol = 1.2,margins = c(6,6))

png("31Var_heatmap.png", height=5.5,width=9,res=600,units="in")
heatmap(t1.m, Rowv = NA,Colv = NA, xlab = "RT", ylab="mz", col = coul, cexRow = 1, cexCol = 1.2,margins = c(6,6))
dev.off()
library(gplots)
library(viridisLite)
library(hrbrthemes)
t1.d$mz <- round(t1.d$mz)
heatmap.2(t1.m, Rowv = NA,Colv = NA, xlab = "RT", ylab="mz", col = coul, cexRow = 1, cexCol = 1.2)
h <- ggplot(t1.d, aes(x = RT, y = mz, fill = count, color = count)) + geom_tile(size=0.2) +
  scale_fill_distiller(direction = 1) + scale_color_distiller(direction=1) +
  theme_minimal() + scale_y_continuous(expand = c(0,0), breaks = seq(400, 1000, by=10)) +
  theme(axis.text = element_text(size = 8), axis.text.x = element_text(size=8, angle=90))
print(h)
dev.off()

heatmap.2(t1.m, Rowv = NA, Colv = NA, dendrogram = NULL )
vcols <- viridis(10, end=0.8)

df_long <- rbind(DIA_short[[2]], DIA_short[[3]]) %>% as.data.table()

sp <- ggplot(df_long, aes(x = RT, y = mz, color=Condition)) + 
  geom_point(size=1,alpha=0.2) + theme_bw()
sp
ggMarginal(sp, type = "histogram",groupColour = T, bins=300, groupFill = T, position=position_dodge())

png("SteppedNCE_Margin.png", width = 8.8, height = 7, res=600, units = "in")

ggMarginal(sp, type = "histogram", fill = "skyblue", bins=300, color="skyblue")
dev.off()


urine_14mz <- fread("D:/projects/pca_urine_spectral_lib/data/openswath/200929_Mstern_DIA/test_14mz_pyprophet_export.tsv")
head(urine_14mz)

urine_14mz <- urine_14mz[!duplicated(urine_14mz, by=c("Charge", "FullPeptideName"))]
tp <- nrow(urine_14mz)
tpep <- urine_14mz$Sequence %>% unique() %>% length()
urine_14mz$ProteinName <- lapply(urine_14mz$ProteinName, function(x) GeneFormat(x))
urine_14mz$ProteinName <- unlist(urine_14mz$ProteinName)
tpro <- urine_14mz$ProteinName %>% unique() %>% length()
tgene <- select(org.Hs.eg.db, columns = cols, keys = unique(urine_14mz$ProteinName), keytype = "UNIPROT")$SYMBOL %>% unique() %>% length()

count <- data.frame(variable=c("Precursor", "Peptide", "Protein", "Genes"), count = c(tp, tpep, tpro, tgene), stringsAsFactors = F)
count$variable <- factor(count$variable, levels = c("Precursor", "Peptide", "Protein", "Genes"))

ggplot(count, aes(x = variable ,y = count)) + 
  geom_bar(stat="identity", width=0.6, fill="lightgoldenrod") + 
  geom_text(aes(label=count), vjust=-0.2, size=5) + xlab("") + ylab("Count") +
  theme_light(base_size=16) 

############# Missing values #######################

missing_values <- lapply(wide.exp, function(x){
  rep <- x[, -1]
  ls <- rowSums(is.na(rep))
  ls <- ls[ls != 3]
  ls
})

missing_exp <- lapply(missing_values, function(x){
  miss1 <- sum(x == 1)
  miss2 <- sum(x == 2)
  c(miss1, miss2)
})

missing <- do.call('rbind', missing_exp) %>% as.data.table()
colnames(missing) <- c("Missing_1", "Missing_2")
missing$Method <- annotation$Condition %>% unique()
missing$Missing_2[4] <- NA
missing$Total <- c(1839, 1916, 1862, 1717)

missing$Missing_1_Percent <- missing$Missing_1 / missing$Total * 100
missing$Missing_2_Percent <- missing$Missing_2 / missing$Total * 100

df_long <- missing %>% subset(select=c("Missing_1_Percent", "Missing_2_Percent", "Method")) %>% melt(id.vars="Method")

ggplot(df_long, aes(x = Method, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = position_dodge(width=0.9)) + theme_bw(base_size = 14) +
  ylab("Percent Missing (Proteins)") + labs(fill = "Missing")
library(gridExtra)

grid.table(NormCV_summary)
alfq.in <- as.data.table(alfq.in)
# Check proteins with at least 2 sequences
StatsByRun <- lapply(unique(alfq.in$run_id), function(x){
  rep <- alfq.in[run_id == x]
  rep <- subset(rep, select=c("protein_id", "peptide_id", "run_id"))
  rep
})

Protein_byRun <- lapply(StatsByRun, function(x){
  rep <- x
  rep <- rep[order(protein_id)]
  rep <- rep[!duplicated(rep, by=c("protein_id", "peptide_id"))]
  t <- table(rep$protein_id) %>% as.data.table()
  t$run <- unique(x$run_id)
  t
})
Protein_byRun <- do.call('rbind', Protein_byRun)
Protein_byRun$Peptides <- ifelse(Protein_byRun$N == 1, "1 peptide", "More than 1 peptides")
Protein_byRun %>% ggplot(aes(x = run, fill=Peptides)) + geom_bar(stat = "count") + theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle=90))


