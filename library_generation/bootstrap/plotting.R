#!/usr/env R

library(data.table)
library(ggplot2)

options(scipen=999)

# load data
files <- list.files('.', pattern = ".tsv")

df <- lapply(files, fread)
df <- do.call('rbind', df)

# format data

df$SampleSize <- gsub("_bootstrap(.*)", "", df$Bootstrap)
df$SampleSize <- as.integer(gsub("[A-Za-z]", "", df$SampleSize))

df$ID <- gsub("(.*)_bootstrap", "", df$Bootstrap)
df$ID <- as.integer(df$ID)

# check summmary of peptide counts
setDT(df)
df[, as.list(summary(`Number of peptides`)), by = SampleSize]
df[, as.list(summary(`Number of proteins`)), by = SampleSize]

pep.summary <- df[, as.list(summary(`Number of peptides`)), by = SampleSize]
pro.summary <- df[, as.list(summary(`Number of proteins`)), by = SampleSize]

pep.sd <- df[, as.list(c(Mean = mean(`Number of peptides`), SD = sd(`Number of peptides`))), by = SampleSize]
pro.sd <- df[, as.list(c(Mean = mean(`Number of proteins`), SD = sd(`Number of proteins`))), by = SampleSize]

# plot
# peptides
timestamp <- Sys.time()
timestamp <- gsub("[^0-9]", "_", timestamp)

pep.plt <- ggplot(df, aes(x = SampleSize, y = `Number of peptides`)) + 
  geom_point() + theme_light(base_size = 14) +
  xlab("Number of runs")

pro.plt <- ggplot(df, aes(x = SampleSize, y = `Number of proteins`)) + 
  geom_point() + theme_light(base_size = 14) +
  xlab("Number of runs")

ggsave(pep.plt, filename=paste0("peptideNumber_bootstrap_", timestamp, ".png"), width = 5, height = 6, unit = "in")

ggsave(pro.plt, filename=paste0("proteinNumber_bootstrap_", timestamp, ".png"), width = 5, height = 6, unit = "in")

# plot mean plus error
pep.sd$SampleSize <- factor(pep.sd$SampleSize, level = c(2,4,6,8,10, 12,14,16,18,20,40, 60, 80, 100, 120, 140, 160, 180))
pep.sumplt <- ggplot(pep.sd) + geom_col(aes(x = SampleSize, y = Mean), width = 0.6, fill = "lightgrey", color = "black") +
  geom_errorbar( aes(x = SampleSize, ymin = Mean-SD, ymax = Mean+SD), width = .2) +
  theme_light(base_size = 14) + xlab("Number of runs") + ylab("Number of peptides") +
  scale_y_continuous(expand = c(0,0), limits = c(0,70000), labels=scales::comma)

pro.sd$SampleSize <- factor(pro.sd$SampleSize, level = c(2,4,6,8,10,12,14,16,18, 20,40, 60, 80, 100, 120, 140, 160, 180))
pro.sumplt <- ggplot(pro.sd) + geom_col(aes(x = SampleSize, y = Mean), width = 0.6, fill = "lightgrey", color = "black") +
  geom_errorbar(aes(x = SampleSize, ymin = Mean-SD, ymax = Mean+SD), width = .2) +
  theme_light(base_size = 14) + xlab("Number of runs") + ylab("Number of proteins")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000), labels=scales::comma)

ggsave(pep.sumplt, filename=paste0("peptideNumber_barplot_bootstrap_", timestamp, ".png"), width = 5, height = 6, unit = 'in')
ggsave(pro.sumplt, filename=paste0("proteinNumber_barplot_bootstrap_", timestamp, ".png"), width = 5, height = 6, unit = 'in')

