#!/usr/bin/env Rscript
#---------------- This script is to be called by bash for normalization ---------------
# load packages
library(data.table)
library(tidyverse)
library(tibble)
library(readxl)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(caret)


# load functions
source("./Utilities.R")
source("../library_generation/cross_validation_loess.R")

#------------- Load command line arguments as files ---------------------------------
args <- commandArgs(trailingOnly = TRUE)
# input: peptide text file (long format with Batch info)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  print("input file and output file required")
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}


#-------------Read data and check columns ------------------------------------------

peptides <- read.delim(args[1], sep = "\t", header = T)


if (!(colnames(peptides) %in% c("Batch", "Run", "log2Intensity"))){
    stop("Missing required columns: Batch, Run, log2Intensity", call.=FALSE)
    }

peptides$ID <- as.integer(str_extract(peptides$Run, "[0-9]+"))

medianBatches <- peptides %>% group_by(Batch, ID) %>%
	summarise(MedianIntensity = median(log2Intensity)) %>% as.data.table()


g <- ggplot(medianBatches, aes(x = ID, y = MedianIntensity)) +
	geom_point(aes(color = Batch) + stat_smooth(method = 'loess', se = F) +
	facet_wrap(~Batch, scales = 'free_x') + theme_light(base_size = 16)

ggsave(g, file = "./median_int_facet.png", width = 10, height = 9)

#----------- Optimize loess spans -------------------------------------------------
set.seed(123)

trainData <- lapply(1:7, function(x){
  df <- medianBatches[Batch == x]
  df
})

span_values <- seq(0, 1, 0.05)

cv_span_values <- lapply(seq_along(trainData), function(x){
  df <- trainData[[x]]
  setnames(df, c("ID", "MedianIntensity"), c("x", "y"))
  mod <- TrainLoessInt(df)
  mod
})


span_results <- lapply(seq_along(trainData), function(x){
  tb <- cv_span_values[[x]]$results
  tb <- as.data.table(tb)
  tb[, run := x]
  tb
})


span_results <- do.call('rbind', span_results)
span_results <- as.data.table(span_results)
setnames(span_results, "run", "Batch")

# check output

g <- span_results %>% mutate(Batch = factor(Batch)) %>%
	ggplot(aes(x = span, y = RMSE, color = Batch)) + geom_point() +
	geom_lin() + theme_classic(base_size = 16)

ggsave(g, file = "./span_results.png", width = 10, height = 8)

# choose the best tune span values

corFactor <- lapply(seq_along(cv_span_values), function(x){
  s <- cv_span_values[[x]]$bestTune$span
  df <- trainData[[x]]
  lo <- loess(y~x, df, span = s, control = loess.control(surface = "interpolate"))
  lo
})

medCorr <- lapply(seq_along(corFactor), function(x){
  mod <- corFactor[[x]]
  X <- mod$x
  y <- mod$y
  fitted <- mod$fitted
  data.frame(Batch = x, ID = X, Intensity = y, fitted = fitted, stringsAsFactors = F )
})

medCorr <- do.call('rbind', medCorr) %>% as.data.table()

# check fit
g <- medCorr %>% mutate(Batch = factor(Batch)) %>%
  ggplot(aes(x = x, y = fitted, color = Batch)) +
  geom_point() + facet_wrap(~Batch, scales = "free_x") + 
  theme_light()
ggsave(g, file = "./loess_fit.png", width = 10, height = 8) 

 
res <- lapply(corFactor, function(x){
  df <- data.frame(ID = x$x, residuals = x$residuals)
  df
})
res <- do.call('rbind', res) %>% as.data.table()
setnames(res, "x", "ID")

#------------------- Use the fitting for correction ---------------------------

peptides <- merge(peptides, res, by = "ID")

peptides$log2Intensity.corrected <- peptides$log2Intensity - peptides$residuals

peptides.int <- subset(peptides, select = c("log2Intensity.corrected", "SampleID", "PeptideSequence"))

peptides.int <- dcast(peptides.int, PeptideSequence~SampleID, value.var = 'log2Intensity.corrected')

#--------------- Global batch normalization with median centering -----------------------

peptides.norm <- medNorm(peptides.int[, 2:ncol(peptides.int)])

peptides.norm <- melt(peptides.norm, id.vars = "PeptideSequence", variable.name = "SampleID",
	value.name = "log2Intensity.Globalcorrected")

peptides.norm <- as.data.table(peptides.norm)
peptides.norm <- peptides.norm[!is.na(peptides.norm)]
peptides <- merge(peptides, peptides.norm, by = c("PeptideSequence", "SampleID"), all.x = T)

remove(peptides.int)
remove(peptides.norm)

save(peptides, file = args[2]) # save as RData format
