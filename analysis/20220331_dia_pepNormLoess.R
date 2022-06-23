#------------------- This is the normalization script ----------------------
# load packages
library(data.table)
library(tidyverse)
library(proBatch)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(tibble)
library(readxl)
source("D:/source/tk_orbitrap_dia/analysis/Utilities.R")

##################### Load data ###################################
dia.peptides <- load("D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/diapeptides_quant_annotated.RData")
load("D:/projects/pca_urine_spectral_lib/data/clinical/uEPS1_patients.RData")

# pseudo code:
# 1. get the peptide medians
# 2. fit the best span of loess within each batch
# 3. correct the bias by subtracting the residual
# 4. General global batch correction by median centering

################# Check data #######################################

head(dia.peptides)
dia.peptides$Log2Intensity.norm <- NULL
dia.peptides %>% 
  group_by(Batch, ID) %>%
  summarise(MedianInt = median(log2Intensity.cor)) %>%
  ggplot(aes(x = ID, y = MedianInt, color = Batch)) + 
  geom_point() +
  stat_smooth(method = 'loess', se = F) +
  scale_color_manual(values = color_list[["Batch"]]) +
  xlab("RunID") + ylab("Median log2Intensity") +
  theme_light(base_size = 16)


medianBatches <- dia.peptides %>% 
  group_by(Batch, ID) %>% 
  summarise(MedianIntensity = median(log2Intensity)) %>% as.data.table()

medianBatches %>%
  ggplot(aes(x = ID, y = MedianIntensity, color = Batch)) +
  geom_point() +
  stat_smooth(method = 'loess', se = F) +
  facet_wrap(~Batch, scales = "free_x") +
  theme_light(base_size = 16)

head(medianBatches)

medianBatches %>%
  arrange(ID) %>%
  group_by(Batch) %>%
  mutate(group_id = cur_group()) %>% View()

############## Optimize span value for each batch ################

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/cross_validation_loess.R")
library(tidyverse)
library(caret)
set.seed(123)


trainData <- lapply(1:7, function(x){
  df <- medianBatches[Batch == x]
  df
})

span_values <- seq(0, 1, 0.05)

cv_span_values <- lapply(1:7, function(x){
  df <- trainData[[x]]
  #setnames(df, c("ID", "MedianIntensity"), c("x", "y"))
  mod <- TrainLoessInt(df)
  mod
})


span_results <- lapply(1:7, function(x){
  tb <- cv_span_values[[x]]$results
  tb <- as.data.table(tb)
  tb[, run := x]
  tb
})


span_results <- do.call('rbind', span_results)
span_results <- as.data.table(span_results)
setnames(span_results, "run", "Batch")
span_results %>% 
  mutate(Batch = factor(Batch)) %>%
  ggplot(aes(x = span, y = Rsquared, color = Batch)) + geom_point() + geom_line() +
  scale_color_manual(values = color_list[["Batch"]])

cv_span_values[[6]]$results
cv_span_values[[6]]$bestTune$span <- 0.35
# use the best tune span value to get correction factor

corFactor <- lapply(1:7, function(x){
  s <- cv_span_values[[x]]$bestTune$span
  df <- trainData[[x]]
  lo <- loess(y~x, df, span = s, control = loess.control(surface = "interpolate"))
  lo
})
corFactor[[6]]$y

medCorr <- lapply(1:7, function(x){
  mod <- corFactor[[x]]
  X <- mod$x
  y <- mod$y
  fitted <- mod$fitted
  data.frame(Batch = x, ID = X, Intensity = y, fitted = fitted, stringsAsFactors = F )
})

medCorr <- do.call('rbind', medCorr) %>% as.data.table()

medCorr %>% mutate(Batch = factor(Batch)) %>%
  ggplot(aes(x = x, y = fitted, color = Batch)) +
  geom_point() + facet_wrap(~Batch, scales = "free_x") + 
  scale_color_manual(values = color_list[['Batch']]) + theme_light()
  
res <- lapply(corFactor, function(x){
  df <- data.frame(ID = x$x, residuals = x$residuals)
  df
})
res <- do.call('rbind', res) %>% as.data.table()
setnames(res, "x", "ID")

######################## Correct samples with the bias factor #################
dia.peptides <- merge(dia.peptides, res, by = "ID")

dia.peptides$log2Intensity.cor <- dia.peptides$log2Intensity - dia.peptides$residuals

dia.peptides %>%
  ggplot(aes(x = ID, y = log2Intensity, color = Batch, group = ID)) +
  geom_boxplot() +
  scale_color_manual(values = color_list[["Batch"]]) +
  theme_light()


dia.noNorm <- dcast(dia.peptides, peptide_id~SampleID, value.var = "log2Intensity")
head(dia.noNorm)
rownames(dia.noNorm) <- dia.noNorm$peptide_id
dia.noNorm$peptide_id <- NULL
dia.noNorm <- getCV(dia.noNorm)
head(dia.noNorm)
dia.noNorm$MAD <- rowMads(as.matrix(dia.noNorm[, 1:189]), na.rm = T)
dia.noNorm <- subset(dia.noNorm, select = c("CV", "mean", "MAD"))
dia.noNorm$peptide_id <- rownames(dia.noNorm)

dia.Norm <- dcast(dia.peptides, peptide_id~SampleID, value.var = "log2Intensity.cor")
head(dia.Norm)
rownames(dia.Norm) <- dia.Norm$peptide_id
dia.Norm$peptide_id <- NULL
dia.Norm <- getCV(dia.Norm)
head(dia.Norm)
dia.Norm$MAD <- rowMads(as.matrix(dia.Norm[, 1:189]), na.rm = T)
dia.Norm <- subset(dia.Norm, select = c("CV", "mean", "MAD"))
dia.Norm$peptide_id <- rownames(dia.Norm)

dia.stats <- merge(dia.noNorm, dia.Norm, by = "peptide_id", suffixes = c(".noNorm", '.Norm'))
head(dia.stats)
dia.stats %>% 
  pivot_longer(-peptide_id,
               names_sep = "\\.",
               names_to = c("Variable", "Normalization")
  ) %>%
  filter(Variable == "CV") %>%
  #filter(peptide_id %in% unique(dia.iRTs$peptide_id)) %>%
  ggplot(aes(x = Normalization, y = value) ) +
  #geom_jitter(width = 0.1) +
  geom_boxplot(alpha = 0.6) +
  ylab("CV") +
  theme_light(base_size = 16)


remove(dia.noNorm)
################### Global batch normalization with median centering ###################
dia.Norm <- medNorm(dia.Norm[, 2:190])
dia.Norm$peptide_id <- rownames(dia.Norm)
dia.Norm

dia.Norm <- melt(dia.Norm, id.vars = 'peptide_id', variable.name = "SampleID", 
                 value.name = "log2Intensity.medCor")
dia.Norm <- as.data.table(dia.Norm)
dia.Norm <- dia.Norm[!is.na(log2Intensity.medCor)]
head(dia.Norm)

dia.peptides <- merge(dia.peptides, dia.Norm, by = c("peptide_id", "SampleID"))
head(dia.peptides)
dia.peptides %>%
  ggplot(aes(x = ID, y= log2Intensity.medCor, color = Batch, group= ID)) +
  geom_boxplot() + 
  scale_color_manual(values = color_list[["Batch"]]) +
  theme_light()


save(dia.peptides, file = "D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/diapeptides_quant_annotated.RData")
dia.pep_medNorm.filtered <- dia.Norm[grepl("_2", peptide_id),]
dia.pep_medNorm.filtered <- dcast(dia.pep_medNorm.filtered,
                                  peptide_id~SampleID, value.var = "log2Intensity.medCor") %>% 
  as.data.frame.matrix()

plot_heatmap_diagnostic(dia.pep_medNorm.filtered, dia.annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE, sample_id_col = "SampleID",
                        heatmap_color = c(brewer.pal(11, "Spectral")),
                        fill_the_missing = 0, color_for_missing = "white")

dia.peptides$Intensity.medCor <- 2^dia.peptides$log2Intensity.medCor

# Perform the same set of normalization on the dEPS and combined DIA data


