##### This is the script for running Retention time alignment ########

library(data.table)
library(ggplot2)
library(tidyr)
library(viridisLite)

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")

basefolder <- "D:/projects/pca_urine_spectral_lib/"

#################### directEPS alignment ####################
nonlinear_lib <- list.files(paste0(basefolder, "data/irt"), pattern = "nonlinear_assaylib.tsv", full.names = T)
nonlinear_lib <- fread(nonlinear_lib)
# From run E156
count_lib <- nonlinear_lib[!duplicated(nonlinear_lib, by = c("ModifiedPeptideSequence", "PrecursorCharge"))]
#### Check stats of the non-linear library
count <- data.frame(variable=c( "Precursors", "Peptides", "ProteinGroups"),
                    counts = c(length(count_lib$ModifiedPeptideSequence %>% unique()),
                               length(count_lib$PeptideSequence %>% unique()),
                               length(count_lib$ProteinId %>% unique())),
                    cohort = "directEPS",
                    stringsAsFactors = FALSE)
count$variable <- factor(count$variable, levels = c("transitions", "precursors", "peptides", "proteins"))

ggplot(count, aes(x = variable, y = counts, label = counts)) + 
  geom_bar(stat="identity", fill=viridis(10)[1]) +
  geom_text(vjust = -1) + theme_bw(base_size = 14)

######## Load Data #####################

pca_eps.msms <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "msms.txt", full.names = T)
pca_eps.evidence <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "evidence.txt", full.names = T)
pca_eps.proteinGroups <- list.files(paste0(basefolder, "data/pca_dda/pca_directEPS"), 
                                    pattern = "proteinGroups.txt", full.names = T)
irt <- list.files(paste0(basefolder, "/data/irt"), pattern="irt_assay", full.names = T)

pca_eps.proteinGroups <- fread(pca_eps.proteinGroups) 
pca_eps.proteinGroups <- subset(pca_eps.proteinGroups,select=c("Protein IDs", "Majority protein IDs", "Gene names"))
pca_eps.msms <- fread(pca_eps.msms)
pca_eps.evidence <- fread(pca_eps.evidence)

irt <- fread(irt)

# 1. Remove the duplicated detections with lower scores
# 2. Perform the linear alignment
# 3. Remove outliers
# 4. Perform the linear alignment without the outliers
# 5. Align the peptides from each run using the linear model
irt <- subset(irt, select = c("PeptideSequence", "NormalizedRetentionTime"))
#setnames(irt, "NormalizedRetentionTime", "iRT")

irt <- irt[!duplicated(irt)]


#################### Run Alignments ######################################

## First check residuals and outliers

runs <- pca_eps.msms$`Raw file` %>% unique()

eps.ev <- subset(pca_eps.evidence, select=c("id", "Intensity"))

eps.msms <- subset(pca_eps.msms, select=c("Evidence ID", "Retention time", "Raw file", "Modified sequence", "Charge", "PEP"))
eps.msms$ModifiedPeptideSequence <- reformat_mods(eps.msms$`Modified sequence`)

eps.ms <- lapply(runs, function(x){
  m <- eps.msms[`Raw file` == x]
  m <- merge(m, eps.ev, by.x="Evidence ID", by.y = "id")
  m
})

# nonlinear matches/alignments
eps.ms_irt <- lapply(eps.ms, function(x){
  m <- x
  m <- m[order(PEP)]
  m <- m[!duplicated(m, by = c("Modified sequence", "Charge"))]
  m <- merge(x, nonlinear_lib, by.x = c("ModifiedPeptideSequence", "Charge"), by.y = c("ModifiedPeptideSequence", "PrecursorCharge"))
  m <- m[!duplicated(m, by = c("ModifiedPeptideSequence", "Charge"))]
  m
})
names(eps.ms_irt) <- runs
# Matched by fragment ions
# Check residuals

########## plot ########################
ggplot(eps.ms_irt[[12]], aes(x =`Retention time`, y = NormalizedRetentionTime )) + 
  geom_point(size=0.5, color = color[2]) + theme_bw(base_size = 16) + xlab("Experimental Retention time (min)") +
  ggtitle("Sample E041") + ylab("Empirical iRT Scale")

########### Continue alignment #########
eps.ms_lnmod <- lapply(runs, function(x){
  m <- eps.ms_irt[[x]]
  m <- subset(m, select=c("NormalizedRetentionTime", "Retention time"))
  lnmod <- lm(NormalizedRetentionTime ~ `Retention time`, m)
  lnmod
})

names(eps.ms_lnmod) <- runs

eps.ms_irt <- lapply(runs, function(x){
  rep <- eps.ms_irt[[x]]
  mod <- eps.ms_lnmod[[x]]
  rep[, residual := mod$residuals]
  rep
})
names(eps.ms_irt) <- runs

eps.residuals <- lapply(runs, function(x){
  mod <- eps.ms_lnmod[[x]]
  tb <- data.frame(residual = mod$residuals, stringsAsFactors = F)
  tb$run <- x
  tb
})

eps.residuals <- do.call('rbind', eps.residuals)
quantile(eps.residuals$residual, probs = c(0.05, 0.95))
outlier <- 5

ggplot(eps.residuals, aes(x = residual, fill = run)) + geom_density(show.legend = F)

pdf("irt_alignment_lnmod.pdf", onefile=T)
color <- viridis(length(eps.ms_irt))
for(run in ms_irt){
  peptides <- nrow(run)
  g <- ggplot(run, aes(x = NormalizedRetentionTime, y = `Retention time`)) +
    geom_point(size=1, color=sample(color, size=1)) + xlab("iRT") + ylab("Retention time (min)") +
    ggtitle(paste0("Alignment with ",peptides, "peptides"))
  print(g)
}
dev.off()


# Check the lnmod before removing residuals

pdf("irt_alignment_withloess.pdf", onefile=T)
for(run in ms_irt){
  peptides <- nrow(run)
  lo <- loess(`Retention time` ~ NormalizedRetentionTime, run)
  df <- data.frame(loess=predict(lo), NormalizedRetentionTime = run$NormalizedRetentionTime, stringsAsFactors = F)
  g <- ggplot(run, aes(x = NormalizedRetentionTime, y = `Retention time`)) +
    geom_point(size=1, color=sample(color, size=1)) +  geom_smooth(method="lm", color = "blue") +
    geom_line(data = df, aes(x = NormalizedRetentionTime, y = loess, color = "red")) +
    xlab("iRT") + ylab("Retention time (min)") +
    ggtitle(paste0("Alignment with ",peptides, " peptides"))
  print(g)
}
dev.off()

lo <- loess(`Retention time` ~ NormalizedRetentionTime, ms_irt[[55]])
df <- data.frame(loess=predict(lo), NormalizedRetentionTime = ms_irt[[55]]$NormalizedRetentionTime, stringsAsFactors = F)
ggplot(ms_irt[[55]], aes(x = NormalizedRetentionTime, y = `Retention time`)) + 
  geom_point(size = 1, color = sample(color, size = 1)) + geom_smooth(method = "lm",aes(color = "linear"), se = F) + 
  geom_line(data=df, aes(x = NormalizedRetentionTime, y = loess, color = "loess")) + 
  xlab("iRT") + ylab("Retention time (min)") + theme_classic() + labs(color = "Legend") + 
  scale_color_manual(values = viridis(2, direction = -1))

resid <- lapply(eps.ms_irt, function(x){
  s <- lm(`Retention time` ~ NormalizedRetentionTime, x)
  res <- s$residuals
  res
})

ms_irt <- lapply(1:length(eps.ms_irt), function(x){
  rep <- eps.ms_irt[[x]]
  rep[,residuals := resid[[x]]]
  rep
})

g <- ggplot(ms_irt[[83]], aes(x = NormalizedRetentionTime, y = residuals)) + geom_point(color = sample(color, 1)) + theme_bw()
g1 <- ggplot(ms_irt[[83]], aes(x = residuals)) + geom_density(fill = sample(color, 1)) + theme_bw()
grid.arrange(g, g1)

pdf("irt_lnmod_residuals.pdf", onefile = T, compress = T)
for (run in ms_irt) {
  g <- ggplot(run)
  g1 <- g + geom_point(aes(x = NormalizedRetentionTime, y = residuals), color=sample(color, 1)) + theme_bw()
  g2 <- g + geom_density(aes(x = residuals), fill = sample(color, 1)) + theme_bw()
  grid.arrange(g1, g2, nrow = 2, newpage = T)
  remove(list = c('g', 'g1', 'g2'))
}
dev.off()
############## check loess span value by cross validation #####
source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/cross_validation_loess.R")
library(tidyverse)
library(caret)
set.seed(123)


trainData <- lapply(runs, function(x){
  df <- eps.ms_irt[[x]]
  df <- df[residual <= 5]
  df <- subset(df, select = c("Retention time", "NormalizedRetentionTime", "ModifiedPeptideSequence", "Charge"))
  setnames(df, c("Retention time", "NormalizedRetentionTime"), c("rt", "irt"))
  df
})

#trainData <- trainData[1:2]
cv_span_values <- lapply(trainData, function(x){
  mod <- TrainLoess(x)
  mod
})

names(cv_span_values) <- runs

library(lattice)
a <- cv_span_values[[1]]

a$results
a$times
a$resample

ggplot(a$results, aes(x = span, y = Rsquared)) + geom_point() + scale_x_reverse() + theme_bw()





eps.all_span_results <- lapply(runs, function(x){
  tb <- cv_span_values[[x]]$results
  tb <- as.data.table(tb)
  tb[, run := x]
  tb
})
eps.all_span_results <- do.call('rbind', eps.all_span_results)
g <- ggplot(eps.all_span_results, aes(x = span, y = Rsquared, color = run)) + geom_point(show.legend = F) +
  geom_line(show.legend = F)  + theme_bw() + scale_x_reverse() + coord_cartesian(xlim = c(0, 0.1))
g1 <- ggplot(eps.all_span_results, aes(x = span, y = MAE, color = run)) + geom_point(show.legend = F) +
  geom_line(show.legend = F) + scale_x_reverse() + theme_bw() + coord_cartesian(xlim = c(0, 0.1))

grid.arrange(g, g1)
############## Check average residual range ###############

# Give maximum 10 minutes span

# Filter the model

outlier <- 10

ms_irt <- lapply(1:length(ms_irt), function(x){
  res <- ms_irt[[x]]
  res <- res[abs(residuals) <= 10]
  res
})

# perform loess alignment
ms_nonlinear <- lapply(1:length(ms_irt), function(x){
  res <- ms_irt[[x]]
  res <- subset(res, select = c("Modified sequence", "Charge", "Retention time", "NormalizedRetentionTime"))
  setnames(res, c("Retention time"), "RetentionTime")
  span <- cv_span_values[[x]]$bestTune$span
  nlmod <- loess(NormalizedRetentionTime ~ RetentionTime, span = span, data = res, surface = "direct")
  nlmod
})



################## Start Alignment #####################

# Extract required information: 
#   Evidence ID, Intensities, Masses, Retention time, Charge, Modified sequence, Sequence,
#   Intensity, m/z

ms <- subset(pca_eps.msms, select = c("Evidence ID", "Raw file", "Modified sequence", "Sequence", 
                             "Retention time", "m/z", "Charge", "Intensities", "Masses", 
                             "PEP", "Proteins", "Gene Names"))
ev <- subset(pca_eps.evidence, select = c("id", "Intensity"))

ms <- merge(ms, ev, by.x = "Evidence ID", by.y = "id")
setnames(ms, c("Evidence ID", "Raw file", "Modified sequence", "Sequence", "Retention time", "Gene Names"),
         c("id", "raw", "ModifiedPeptideSequence", "PeptideSequence", "RetentionTime", "UniProtID"))

ms <- ms[order(PEP)]
ms$ModifiedPeptideSequence <- reformat_mods(ms$ModifiedPeptideSequence)
ms <- ms[!duplicated(ms, by = c("ModifiedPeptideSequence", "Charge"))]
# 62519

############# Count Number of library #################
df <- data.frame(count = numeric(length = 3), variable =c("Precursor", "Peptide", "Protein"), stringsAsFactors = F)
df$count[1] <- nrow(ms)
df$count[2] <- ms$PeptideSequence %>% unique() %>% length()
df$count[3] <- ms$Proteins %>% unique() %>% length()

ggplot(df, aes(x = variable, y = count)) + geom_bar(stat="identity")

library(gridExtra)
png("directEPS_libcount.png", height = 50*nrow(df), width = 200*ncol(df))
grid.table(df)
dev.off()
############# Start Alignment ########################

# Align each fraction using the previously computed loess model

ms <- lapply(1:length(runs), function(x){
  res <- ms[raw == runs[x]]
  mod <- ms_nonlinear[[x]]
  nRT <- predict(mod, newdata = res, se = TRUE)
  res$NormalizedRetentionTime <- nRT$fit
  res$StandardError <- nRT$se.fit
  res
})

################ calculate R square ###############

eps_rsquare <- lapply(1:length(runs), function(x){
  res <- ms[raw == runs[x]]
  r2 <- cor(res$RetentionTime, res$NormalizedRetentionTime)^2
  r2
})

r2 <- data.frame(run = runs, R_square = unlist(eps_rsquare), stringsAsFactors = F)
r2$cohort <- "directEPS"
ggplot(r2, aes(x = R_square)) + geom_density()

############## Combine results ###################

ms <- do.call('rbind', ms)
ggplot(ms, aes(x = RetentionTime, y = NormalizedRetentionTime, color = raw)) + geom_point(show.legend = F)
ggplot(ms, aes(x = as.numeric(factor(raw)), y = StandardError, group = raw)) +
  geom_boxplot(width = 0.5, show.legend = F)

ggplot(ms, aes(x = StandardError)) + geom_density()
ms[max(ms$StandardError)]
sum(is.na(ms$NormalizedRetentionTime))
quantile(ms$StandardError, c(0.025, 0.975))
# Remove points over 1

ms <- ms[StandardError < 1.0]
nrow(ms)
# 62,448

write.table(ms, file = "D:/projects/pca_urine_spectral_lib/data/library/directEPS_nonlinear_lib_se_filtered_20200402.tsv", sep = "\t", row.names = F, quote = F)

##################### Urine Alignment #####################

# Load data
urine_irt <- list.files(paste0(basefolder, "/data/irt/"), pattern = "osw_assaylib", full.names = T)
urine_irt <- urine_irt %>% fread()

urine_msms <- list.files(paste0(basefolder, "/data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og"), 
                         pattern = "msms.txt", full.names = T) %>% fread()

urine_ev <- list.files(paste0(basefolder, "/data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og"), 
                         pattern = "evidence", full.names = T) %>% fread()

urine_ev <- subset(urine_ev, select = c("id", "Intensity"))

urine_ms <- merge(urine_msms, urine_ev, by.x = "Evidence ID", by.y = "id")

urine_files <- urine_ms$`Raw file` %>% unique()

########## Normalization with hpirt ###################

urine_ms <- urine_ms[order(PEP)]

urine_ms$ModifiedPeptideSequence <- reformat_mods(urine_ms$`Modified sequence`)
urine_irt <- urine_irt[!duplicated(urine_irt, by = c("ModifiedPeptideSequence", "PrecursorCharge"))]
urine_irt <- subset(urine_irt, select = c("PrecursorMz", "PrecursorCharge", "LibraryIntensity", 
                                          "NormalizedRetentionTime", "PeptideSequence", 
                                          "ModifiedPeptideSequence", "TransitionGroupId"))
urine_msms_irt <- lapply(urine_files, function(x){
  rep <- urine_ms[`Raw file` == x]
  rep <- rep[order(PEP)]
  rep <- rep[!duplicated(rep, by = c("ModifiedPeptideSequence", "Charge"))]
  rep <- subset(rep,select = c("ModifiedPeptideSequence", "Charge", "id", "Retention time", "m/z", "PEP"))
  rep <- merge(rep, urine_irt, by.x = c("ModifiedPeptideSequence", "Charge"), 
               by.y = c("ModifiedPeptideSequence", "PrecursorCharge"))
  rep
})

names(urine_msms_irt) <- urine_files

pdf(paste0(basefolder, "results/2020Feb_pca_urine/urine_msms_irt.pdf"), onefile = T)

for (r in urine_files) {
  g <- ggplot(urine_msms_irt[[r]], aes(x = `Retention time`, y = NormalizedRetentionTime)) + geom_point(size=1)
  g <- g + theme_light(base_size = 12) + 
    ggtitle(r, subtitle = paste0(nrow(urine_msms_irt[[r]]), " peptide for alignment"))
  print(g)
}

dev.off()


############### check residuals and outliers ##################
urine_lnmod <- lapply(urine_files, function(x){
  rep <- urine_msms_irt[[x]]
  rep <- subset(rep, select=c("NormalizedRetentionTime", "Retention time"))
  lnmod <- lm(NormalizedRetentionTime ~ `Retention time`, rep)
  lnmod
})

names(urine_lnmod) <- urine_files

urine_absd <- lapply(urine_files, function(x){
  rep <- urine_lnmod[[x]]
  r <- rep$residuals
  df <- data.table(residual = r)
  df[, run := x]
  df
})

urine_absd <- do.call('rbind', urine_absd)

quantile(urine_absd$residual, probs = c(0.05, 0.95))

#5%       95% 
#-3.486216  2.199035 

urine_msms_irt <- lapply(urine_files, function(x){
  ms <- urine_msms_irt[[x]]
  r <- urine_lnmod[[x]]$residual
  ms$residual <- r
  ms
})

names(urine_msms_irt) <- urine_files

color <- viridis(option = "plasma", n = length(urine_files))
pdf(paste0(basefolder, "results/2020Feb_pca_urine/urine_irt_residuals.pdf"), onefile = T)
for (r in urine_files) {
  g <- ggplot(urine_msms_irt[[r]])
  g1 <- g + geom_point(aes(x = NormalizedRetentionTime, y = r), color=sample(color, 1)) + theme_bw() +
    ggtitle(r)
  g2 <- g + geom_density(aes(x = r), fill = sample(color, 1)) + theme_bw()
  grid.arrange(g1, g2, nrow = 2, newpage = T)
  remove(list = c('g', 'g1', 'g2'))
}
dev.off()

ggplot(urine_absd[run %in% urine_files[1:2]], aes(x = residual, fill = run)) + geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(-10, 10)) + scale_fill_viridis_d() + theme_bw()

# R002, R003 -> wrong gradient
# R050  -> column spitting
# R100 -> only 1000 peptides (sparse)

# Contaminents: 362.15 mz, 389.16 mz # take them out

#################### Cross Validation of span values ######################
library(tidyverse)
library(caret)
set.seed(123)

trainData_urine <- lapply(urine_files, function(x){
  df <- urine_msms_irt[[x]]
  df <- df[residual <= 5]
  df <- subset(df, select = c("Retention time", "NormalizedRetentionTime", "ModifiedPeptideSequence", "Charge"))
  setnames(df, c("Retention time", "NormalizedRetentionTime"), c("rt", "irt"))
  df
})
names(trainData_urine) <- urine_files

#trainData_urine <- trainData_urine[1:2]

#mod <- TrainLoess(trainData_urine[[urine_files[9]]])

#urine_span_mod <- lapply(urine_files[5:10], function(x){
#  df <- trainData_urine[[x]]
#  mod <- TrainLoess(irt_data = df)
#  mod
#})
#names(urine_span_mod) <- urine_files[5:10]
#results <- lapply(urine_files[5:10], function(x){
##  r <- urine_span_mod[[x]]$results %>% as.data.table()
#  r[, run := x]
#  r
#})

#results <- do.call('rbind', results)
#g <- ggplot(results[run == urine_files[5]], aes(x = span, y = Rsquared, color = run)) + geom_point(show.legend = F) +
#  geom_line(show.legend = F)

#g1 <- ggplot(results[run == urine_files[6]], aes(x = span, y = Rsquared, color = run)) + geom_point(show.legend = F) +
#  geom_line(show.legend = F)

#grid.arrange(g, g1)

# Used the actual msms_irt, originally used the irt only which explains why it looked really nice

# plan: output the best span value from each model, and then each fraction is aligned with the specific model

cv_span_values <- lapply(trainData_urine, function(x){
  mod <- TrainLoess(x)
  mod
})

names(cv_span_values) <- urine_files


urine_all_span_results <- lapply(urine_files, function(x){
  tb <- cv_span_values[[x]]$results
  tb <- as.data.table(tb)
  tb[, run := x]
  tb
})

urine_all_span <- lapply(urine_files, function(x){
  s <- cv_span_values[[x]]$bestTune$span
  c(x, s)
})

urine_all_span_results <- do.call('rbind', urine_all_span_results)
urine_span <- do.call('rbind', urine_all_span) %>% as.data.table()
colnames(urine_span) <- c("run", "span")

g <- ggplot(urine_all_span_results, aes(x = span, y = Rsquared, color = run)) + geom_point(show.legend = F) +
  geom_line(show.legend = F)  + theme_bw() + scale_x_reverse() + coord_cartesian(xlim = c(0, 0.1))
g1 <- ggplot(urine_all_span_results, aes(x = span, y = MAE, color = run)) + geom_point(show.legend = F) +
  geom_line(show.legend = F) + scale_x_reverse() + theme_bw() + coord_cartesian(xlim = c(0, 0.1))

grid.arrange(g, g1)


# all replicates showed 0.02 is the best span parameter

##################### Run alignment #################

outlier <- 5

urine_msms_irt <- lapply(urine_files, function(x){
  rep <- urine_msms_irt[[x]]
  rep <- rep[abs(residual) <= 5]
  rep
})

names(urine_msms_irt) <- urine_files

urine_nonlnmod <- lapply(urine_files, function(x){
  rep <- urine_msms_irt[[x]]
  s <- cv_span_values[[x]]$bestTune$span
  rep <- subset(rep, select=c("Retention time", "NormalizedRetentionTime", 
                              "ModifiedPeptideSequence", "PeptideSequence","Charge", "m/z"))
  setnames(rep, "Retention time", "RetentionTime")
  mod <- loess(NormalizedRetentionTime ~ RetentionTime, data = rep, 
               span = s, control = loess.control(surface = "direct"))
  mod
})

names(urine_nonlnmod) <- urine_files

# Predict the runs
urine_ms <- urine_ms[order(PEP)]
#urine_ms <- urine_ms[!duplicated(urine_ms, by = c("ModifiedPeptideSequence", "Charge"))]
urine_ms_sub <- lapply(urine_files, function(x){
  res <- urine_ms[`Raw file` == x]
  res <- res[!duplicated(res, by = c("ModifiedPeptideSequence", "Charge"))]
  setnames(res, "Retention time", "RetentionTime")
  res
})

names(urine_ms_sub) <- urine_files

prediction <- lapply(urine_files, function(x){
  mod <- urine_nonlnmod[[x]]
  sub <- urine_ms_sub[[x]]
  newRT <- predict(mod, newdata = sub, se = TRUE)
  newRT
})

names(prediction) <- urine_files
urine_ms_sub <- lapply(urine_files, function(x){
  ms <- urine_ms_sub[[x]]
  newRT <- prediction[[x]]
  ms[, NormalizedRetentionTime := newRT$fit]
  ms[, StandardError := newRT$se.fit]
  ms
})

names(urine_ms_sub) <- urine_files
colors <- viridis(length(urine_files), option = "plasma")

pdf(paste0(basefolder, "results/2020Feb_pca_urine/fit_se.pdf"), onefile = T)

for (run in seq_along(urine_files)) {
  sub <- urine_ms_sub[[run]]
  g <- ggplot(sub, aes(x = RetentionTime, y = StandardError)) + 
    geom_point(color = colors[run]) + theme_bw(base_size = 12) + ylim(0, 1.0)
  g1 <- ggplot(sub, aes(x = StandardError)) + geom_density(color = colors[run]) + theme_bw(base_size = 12) +
    geom_vline(xintercept = quantile(sub$StandardError, probs = c(0.05, 0.95))) + xlim(0, 1.0)
  grid.arrange(g, g1, top = urine_files[run], newpage = T)
  remove(list=c('g', 'g1'))
}

dev.off()

############# Combine alignment results #####################

# no outlier removal:

urine_ms_nrt <- do.call('rbind', urine_ms_sub)
# 2421846 rows
urine_ms_nrt <- urine_ms_nrt[order(PEP)]
urine_ms_nrt <- urine_ms_nrt[!duplicated(urine_ms_nrt, by = c("ModifiedPeptideSequence", "Charge"))]
# 78179 rows
setnames(urine_ms_nrt, "Raw file", "raw")
ggplot(urine_ms_nrt, aes(x = NormalizedRetentionTime, y = `RetentionTime`, color = `raw`)) +
  geom_point(show.legend = F) + scale_color_viridis_d()
ggplot(urine_ms_nrt, aes(x = as.numeric(factor(raw)), y = StandardError, group = raw)) + 
  geom_boxplot(width = 0.5, show.legend = F)
quantile(urine_ms_nrt$StandardError, c(0.025, 0.975))
sum(is.na(urine_ms_nrt$NormalizedRetentionTime))
# with outlier removal:

urine_ms_nrt <- lapply(urine_ms_sub, function(x){
  rep <- x
  rep <- rep[StandardError <= 1]
  rep
})


urine_ms_nrt <- do.call('rbind', urine_ms_nrt)
urine_ms_nrt <- urine_ms_nrt[order(PEP)]
urine_ms_nrt <- urine_ms_nrt[!duplicated(urine_ms_nrt, by = c("ModifiedPeptideSequence", "Charge"))]
# 78092 rows

ggplot(urine_ms_nrt, aes(x = NormalizedRetentionTime, y = `RetentionTime`, color = `Raw file`)) +
  geom_point(show.legend = F) + scale_color_viridis_d()

urine_ms_nrt$ModifiedPeptideSequence %>% unique() %>% length()
# 56546
urine_ms_nrt$Sequence %>% unique() %>% length()
# 52027

urine_ms_nrt$Proteins %>% unique() %>% length()

# 6962

urine_ms_nrt$`Gene Names` %>% unique() %>% length()
# 4625

############## Write table ##########################
urine_lib <- subset(urine_ms_nrt, select = c("Evidence ID", "Sequence", "ModifiedPeptideSequence",
                                             "Proteins", "Gene Names", "Charge", "m/z", "PEP",
                                             "Intensities", "Masses", "NormalizedRetentionTime"))

setnames(urine_lib, c("Evidence ID", "Sequence", "Gene Names"), c("id", "PeptideSequence", "UniProtID"))

write.table(urine_lib, file = paste0(basefolder, "data/library/urine_nonlinear_lib_se_filtered_20200403.tsv"), sep = "\t",
            quote = F, row.names = F)

