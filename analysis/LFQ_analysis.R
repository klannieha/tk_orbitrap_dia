#!/bin/usr/ R
# module load gcc/7.3.0
# module load r/3.6.0
# Better to request for enough memory !!


library(data.table)
library(ggplot2)
library(tidyr)
library(MSstats)
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

data <- '/project/6002011/annieha/pca_urine_spectral_lib/data/openswath/200903_HEK_DIA'

# loading the runs from Aug 27 first
# Notes: pyprophet output with protein IDs deprecated

files <- list.files(data, pattern="_0914lib_all_merged_pyprophet_export.tsv")

stats <- lapply(files, function(x)  fread(x, sep = "\t"))

# load experiment conditions with corresponding filenames
annotation <- fread("annotation.txt")
setnames(annotation, "Run", "filename")

lib.path <- '/project/6002011/annieha/pca_urine_spectral_lib/data/library/HEK/200914_tpp_HEK_assaylib_osw_target_decoy.tsv'
lib <- fread(lib.path)

# first fix the protein IDs by matching to the library protein IDs
#proteinIds <- lib$ProteinId %>% unique()

#dep_proteinIds <- gsub(";", "", proteinIds)

#lib.protein <- data.frame(proteinIds, dep_proteinIds, stringsAsFactors=F)

# Analysis: basic peptide and protein counts:

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

count <- lapply(annotation$filename, function(x){
  rep <- stats[filename == x]
  rep <- rep[!duplicated(rep, by = c("FullPeptideName", "Charge", "transition_group_id"))]
  precursor <- nrow(rep)
  peptide <- rep$Sequence %>% unique() %>% length()
  rep <- rep[!grepl("DECOY", rep$ProteinName)]
  proteinGroups <- rep$ProteinName %>% unique() %>% length()
  genes <- rep[!grepl("DECOY|;", rep$GENE)]$GENE %>% unique() %>% length()
  c <- c(precursor, peptide, proteinGroups, genes)
  c})

count <- do.call('rbind', count) %>% as.data.table()
colnames(count) <- c("Precursor", "Peptide", "Protein", "Genes")
count <- cbind(count, annotation)
count <- count[1:8]

count <- melt(count, id.vars=c("Condition", "BioReplicate", "filename"), value.name="count")
count$BioReplicate <- factor(count$BioReplicate, levels=c(1,2))

g <- ggplot(count, aes(x = Condition, y = count, fill = BioReplicate)) +
geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.6) +
geom_text(aes(label = count), position=position_dodge(width=0.8), vjust= -0.1) +
facet_wrap(.~variable, scales="free_y") + theme_light() +
scale_fill_viridis_d(option="C", end=0.6)

png('/project/6002011/annieha/pca_urine_spectral_lib/results/2020Sept_HEKdia/count.png',
width=800)
print(g)
dev.off()


# Start with protein inference and quantification

stats$UniprotID <- stats$ProteinName
stats$ProteinName <- stats$GENE
stats$GENE <- NULL

Allstats <- reduce_OpenSWATH_output(stats) %>% as.data.frame(stringsAsFactors=F)
#setnames(annotation, "filename", "Filename")
annotation$Run <- 1:nrow(annotation)
annotation <- annotation[1:8,] %>% as.data.frame(stringsAsFactors=F)
Allstats <- sample_annotation(Allstats, annotation)

# reformat for aLFQ using SWATH2stats functions
disagg <- disaggregate(Allstats)
alfq.in <- convert4aLFQ(disagg, check_transitions = FALSE)

# run aLFQ for protein inference

prots <- ProteinInference(alfq.in, peptide_method = "top", peptide_topx = 1, peptide_strictness = "strict",peptide_summary = "sum", transition_topx = 5, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, combine_peptide_sequences = TRUE, consensus_proteins = FALSE, consensus_peptides = FALSE, consensus_transitions = FALSE)

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

wide.exp <- lapply(annotation$Condition, function(x){
   c <- colnames(wide)
   exp <- c[grepl(x, c)]
   rep <- subset(wide, select=c("protein_id", exp))
   rep})


# log transform

wide.log <- log(wide[, -1])
rownames(wide.log) <- wide$protein_id
# normalize in log space


wide.lognorm <- medNorm(wide.log)
wide.lognorm.median.abund <- apply(wide.lognorm, 1, median, na.rm=T) # get normalized medians

rownames(wide.lognorm.median.abund) <- wide$protein_id

write.table(wide.lognorm.allmedian.abund, file="/project/6002011/annieha/pca_urine_spectral_lib/results/2020Sept_HEKdia/prot_abund_med_norm.tsv", sep="\t", row.names = T,quote=F)

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

png("Normalized_Logintensity.png")
print(g)
dev.off()

# check protein abundance across experiments
exp <- colnames(wide.lognorm.allmedian.abund)

wide.lognorm.allmedian.abund <- subset(wide.lognorm.allmedian.abund,
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

png("NormCV_boxplot.png")
print(g)
dev.off()


g <- ggplot(NormCV, aes(y = log(mean), x = CV, color=Method)) + geom_point(alpha=0.5) + scale_color_viridis_d(option="C", end=0.8) 
png("CV_meanInt.png")
print(g)
dev.off()

NormCV_summary <- NormCV %>% group_by(Method) %>% summarize(mean=mean(CV, na.rm=T), median=median(CV, na.rm=T)) %>% as.data.table()

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
