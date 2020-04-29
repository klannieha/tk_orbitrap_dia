#!/bin/env R

# This is a script for anaalyzing the differences between runs within 1
# experimment. Identifying the spans between irt peptides and the other common
# peptides
library(viridisLite)
library(data.table)
library(ggplot2)
library(tidyr)
library(VennDiagram)
setwd("D:/projects/pca_urine_spectral_lib")
#source("/home/annieha/source/tk_orbitrap_dia/analysis/Utilities.R")

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")

basefolder <- "D:/projects/pca_urine_spectral_lib"

# test out the data first

pca_exo <- list.files(paste0(basefolder, "/data/pca_dda/pca_urine_exo"),
  full.names=T)

pca_exo.ev <- fread(pca_exo[1])
pca_exo.msms <- fread(pca_exo[2])
pca_exo.pep <- fread(pca_exo[3])
pca_exo.protein <- fread(pca_exo[4])


### Identify iRT peptides and the SUC2 peptides

irt <- list.files(paste0(basefolder, "/data/irt"), pattern="assay", 
  full.names=T)

irt <- fread(irt)

pca_exo.submsms <- subset(pca_exo.msms, select=c("Raw file", "Sequence",
	"Modified sequence", "Proteins", "Charge", "Retention time",
	"Precursor Intensity", "m/z"))

setnames(pca_exo.submsms, colnames(pca_exo.submsms), c("raw",
  "PeptideSequence", "ModifiedPeptideSequence", 
  "Proteins", "PrecursorCharge", "RetentionTime", 
  "Intensity", "PrecursorMz"))

pca_exo.submsms$ModifiedPeptideSequence <- reformat_mods(pca_exo.submsms$ModifiedPeptideSequence)

irt <- subset(irt, select=c("ModifiedPeptideSequence", "NormalizedRetentionTime"))
irt <- irt[!duplicated(ModifiedPeptideSequence)]
pca_exo.msms_irt <- merge(pca_exo.submsms, irt, by="ModifiedPeptideSequence")

raw_files <- pca_exo.submsms$raw %>% unique()

#pdf(paste0(basefolder,"/results/pca_urine_exo_irt_all.pdf"),onefile = TRUE)

for(r in raw_files){
  msms <- pca_exo.msms_irt[raw == r]
#  msms <- msms[order(PEP)]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  p <- ggplot(msms, aes(x = NormalizedRetentionTime, y = RetentionTime)) + geom_point(size=1, alpha=0.6)
	geom_abline(slope=1)+ theme_bw()
  print(p)
}
dev.off()


png(paste0(basefolder,"/results/pca_urine_exo_irt_all.png"))
p <- ggplot(pca_exo.msms_irt, aes(x = NormalizedRetentionTime,
  y = RetentionTime, color = raw)) + geom_point(size=1, alpha=0.6, show.legend=F) + 
  geom_line(show.legend=F) +theme_bw()
print(p)
dev.off()

########## Calcute the rsq values between the iRT and the exp #########

rsq <- unlist(lapply(raw_files, function(x){
  msms <- pca_exo.msms_irt[raw == x]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  ln <- lm(NormalizedRetentionTime ~ RetentionTime, msms)
  r <- summary(ln)$r.squared
  r
}))

rsq_tb <- data.frame(r_squared = rsq, run = 1:155)

p <- ggplot(rsq_tb, aes(x = run, y = r_squared)) + geom_point() + theme_bw()
+ ylim(0, 1)
ggsave(file=paste0(basefolder, "/results/pca_urine_exo_irt_rsq.png"), plot =print(p))


lnmod <- lapply(raw_files, function(x){
  msms <- pca_exo.msms_irt[raw == x]
  msms <- msms[!duplicated(ModifiedPeptideSequence, PrecursorCharge)]
  ln <- lm(iRT ~ RetentionTime, msms)
  ln
})

# Look at the distribution of RT diff

#pdf(paste0(basefolder, "/results/pca_urine_exo_irt_resid.pdf"), onefile=TRUE)

p <- ggplot(pca_exo.msms_irt, aes(x = PeptideSequence, y = RetentionTime)) +
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 30))
ggsave("pca_urine_exo_irt_box.png", plot = print(p))


setnames(pca_exo.msms_irt, "NormalizedRetentionTime", "iRT")

pca_exo.msms_irt_filtered <- lapply(raw_files, function(x){
  sub <- pca_exo.msms_irt[raw == x]
  sub <- sub[order(-Intensity)]
  sub <- sub[!duplicated(PeptideSequence)]
  sub
})

pca_exo.msms_irt_filtered <- do.call('rbind', pca_exo.msms_irt_filtered)
pca_exo.msms_irt <- pca_exo.msms_irt_filtered %>% subset(select = c("ModifiedPeptideSequence", 
  "raw", "RetentionTime"))

pca_exo.msms_irt <- reshape(pca_exo.msms_irt,idvar = "ModifiedPeptideSequence", timevar="raw",
  direction="wide")

colnames(pca_exo.msms_irt) <- c("Peptide", raw_files)

sd_na <- lapply(1:11, function(x){
  s <- sd(pca_exo.msms_irt[x, 2:156], na.rm = T)
  na <- sum(is.na(pca_exo.msms_irt[x, 2:156]))
  c(s,na)
})

sd_na <- do.call('rbind', sd_na) %>% as.data.frame()

pca_exo.msms_irt$sd <-sd_na$V1
pca_exo.msms_irt$na <- sd_na$V2


############## Look at the common peptides ###########################

pca_exo.msms_sub <- lapply(raw_files, function(x){
  df <- pca_exo.msms[`Raw file` == x]
  df <- df[order(PEP)]
  df <- subset(df, select=c("Raw file", "Proteins", "Gene Names", "Modified sequence",
	 "m/z", "Charge", "Retention time", "PEP", "Precursor Intensity", "Sequence", 
	"Masses", "Intensities"))
  setnames(df, colnames(df), c("raw", "proteins", "genes", "ModifiedPeptideSequence", 
	"PrecursorMz", "Charge", "RetentionTime", "pep", "Intensity", "PeptideSequence", 
	"masses", "intensities"))
  df <- df[!duplicated(ModifiedPeptideSequence, Charge)]
  df
})

pep_id_names <- colnames(pca_exo.pep)
pep_id_names <- pep_id_names[grep("Identification*", pep_id_names)]
pca_exo.pep_id <- subset(pca_exo.pep, select=c("Sequence", pep_id_names))

pca_pep_count <- lapply(pca_exo.msms_sub, function(x){
  df <- x
  count <- df$PeptideSequence %>% unique() %>% length()
  c(count, unique(df$raw))
})

pep_count_run <- do.call('rbind', pca_pep_count)
pep_count_run <- pep_count_run %>% as.data.frame(stringsAsFactors=FALSE)
colnames(pep_count_run) <- c("count", "run") 
pep_count_run$count <- as.numeric(pep_count_run$count)
 

pca_pep <- lapply(pca_exo.msms_sub, function(x){
  df <- x
  sequence <- df$PeptideSequence %>% unique()
  sequence
})

names(pca_pep) <- raw_files

common <- Reduce(intersect, pca_pep)
q <- quantile(pep_count_run$count)
g <- ggplot(pep_count_run, aes(x = count)) + geom_density() + geom_vline(xintercept = q, linetype = "dotted") + annotate(geom="text", x=q, y=0, label=names(q)) + theme_bw(base_size=16)

# quantile

g <- ggplot(pep_count_run, aes(x = count)) + stat_ecdf(geom = "step") + theme_bw(base_size=16)

# Upper 50% quantile: 7750 peptides
##### Further filter by upper 50% quantile:
pep_count_run <- pep_count_run %>% as.data.table()
run_id <- pep_count_run[count >= q[3]]$run # 78 runs

names(pca_exo.msms_sub) <- raw_files

pca_exo.msms_upper <- pca_exo.msms_sub[run_id]

pca_upper_pep <- pca_pep[run_id]
common <- Reduce(intersect, pca_upper_pep)


upper_pepdf <- lapply(pca_exo.msms_upper, function(x){
  run <- x
  run <- run[PeptideSequence %in% common]
  run <- run[order(pep)]
  run <- run[order(-Intensity)]
  run$ModifiedPeptideSequence <- reformat_mods(run$ModifiedPeptideSequence)
  run <- run[!duplicated(ModifiedPeptideSequence, Charge)]
  run <- subset(run, select=c("raw", "proteins", "genes", "ModifiedPeptideSequence", "PrecursorMz", "Charge", "RetentionTime", "PeptideSequence", "Intensity"))
  run
})

statsdf <- lapply(1:length(run_id), function(x){
  run <- upper_pepdf[[x]]
  df <- subset(run, select=c( "ModifiedPeptideSequence","RetentionTime"))
  setnames(df, "RetentionTime", run_id[x])
  df
})

pep_stats <- Reduce(function(x, y) merge(x, y, by = "ModifiedPeptideSequence", all=TRUE), statsdf)

sd_na <- lapply(1:nrow(pep_stats), function(x){
  s <- sd(pep_stats[x, 2:ncol(pep_stats)], na.rm = T)
  na <- sum(is.na(pep_stats[x, 2:ncol(pep_stats)]))
  c(s,na)
})

sd_na <- do.call('rbind', sd_na) %>% as.data.frame()

pep_sd <- data.frame(ModifiedPeptideSequence = pep_stats$ModifiedPeptideSequence, SD = sd_na$V1, na = sd_na$V2)

g <- ggplot(pep_sd, aes(x = SD)) + geom_density() + theme_bw(base_size=16) +scale_fill_viridis_c(option="plasma")


pdf("RT_alignment_common_peptides.pdf",onefile = TRUE)

exclude <- run_id[-6]

colors <- viridis(option="plasma", n = 77)

for(i in 1:length( exclude)){
  r <- exclude[i]
  cols <- c(run_id[6], r)
  df <- subset(pep_stats, select=cols)
  g <- ggplot(df, aes_string(x = cols[1], y = cols[2])) + geom_point(size = 1, col = colors[i]) + geom_abline(slope = 1) + theme_bw(base_size=16) + xlim(0, 150) + ylim(0, 150)
  print(g)
}
dev.off()


r <- exclude[1]
cols <- c(run_id[6], r)
df <- subset(pep_stats, select = cols)
g <- ggplot(df, aes_string(x = cols[1], y = cols[2])) + geom_point(size=1, col = colors[1]) + xlim(0, 150) + ylim(0, 150)

for(i in 2:length(exclude)){
  r <- exclude[i]
  cols <- c(run_id[6], r)
  df <- subset(pep_stats, select=cols)
  g <- g + geom_point(data=df, aes_string(x = cols[1], y = cols[2]), size = 1, col = colors[i])
}

g <- g + geom_abline(slope = 1) + theme_bw(base_size=16)

ggsave("allRT_common_peptides.png", plot = print(g))

########################## Take runs with 75% quant cutoff peptide count###############
run_id <- pep_count_run[count >= 10000]$run

msms_runs <- pca_exo.msms_sub[run_id] # 46 runs only

msms_runs <- lapply(msms_runs, function(x){
  seq <- x$ModifiedPeptideSequence
  seq <- reformat_mods(seq)
  x$ModifiedPeptideSequence <- seq
  x
})

seqlst <- lapply(msms_runs, function(x) x$ModifiedPeptideSequence)
common <- Reduce(intersect, seqlst)

msms_irt <- lapply(1:length(msms_runs), function(x){
  run <- msms_runs[[x]]
  run <- run[ModifiedPeptideSequence %in% common]
  run <- subset(run, select=c("ModifiedPeptideSequence", "RetentionTime"))
  setnames(run, "RetentionTime", run_id[x])
  run
})

msms_irt <- Reduce(function(x,y) merge(x,y, by="ModifiedPeptideSequence", all = TRUE), msms_irt)

sd_na <- lapply(1:nrow(msms_irt), function(x){
  s <- sd(msms_irt[x, 2:ncol(msms_irt)], na.rm = T)
  na <- sum(is.na(msms_irt[x, 2:ncol(msms_irt)]))
  c(s,na)
})

sd <- do.call("rbind", sd_na)
sd <- sd %>% as.data.table()

g <- ggplot(sd, aes(x = SD)) + geom_density(alpha=0.8, fill = sample(colors, size=1)) + theme_bw(base_size=16)
ggsave("SD_75quant_filtered.png", plot = print(g))
 

########### Look at the residuals of each run to find the span value #####################
# run number 137 has the lowest sd
# Use CiRT

# ONLY 8 PEPTIDES OVERLAPPED ----- ABORT

msms_runs_irt <- lapply(msms_runs, function(x){
  rep <- x
  rep <- subset(rep, select=c("raw", "RetentionTime", "PeptideSequence", "Charge", "Intensity"))
  rep <- merge(rep, irt, by.x = "PeptideSequence", by.y = "ModifiedPeptideSequence", all.y = TRUE)
  rep <- rep[order(NormalizedRetentionTime)]
  rep
})

col <- sample(colors, size=46)
pdf("top75_runs_msms_irt.pdf", onefile=TRUE)
for(run in 1:length(msms_runs_irt)){
  rep <- msms_runs_irt[[run]]
  g <- ggplot(rep, aes(x = NormalizedRetentionTime, y = RetentionTime)) + geom_point(size=1, color=col[run]) + theme_bw(base_size=16) + xlab("iRT")
  print(g)
}

dev.off()

msms_runs_norm <- do.call('rbind', msms_runs_irt)
msms_runs_norm <- as.data.table(msms_runs_norm)
msms_runs_norm <- msms_runs_norm[order(NormalizedRetentionTime)]

tips2$day <- factor(tips2$day,levels = c("Fri", "Sat", "Sun", "Thur"))
irt <- irt[order(NormalizedRetentionTime)]
msms_runs_norm$irt_id <- factor(msms_runs_norm$PeptideSequence, levels = irt$ModifiedPeptideSequence)


g <- ggplot(msms_runs_norm, aes(x = irt_id, y = RetentionTime)) + geom_boxplot(width=0.5) + theme_bw(base_size=15) + theme(axis.text.x = element_text(angle = 30))

ggsave("top75runs_irt.png")

######### Add these peptides to the iRT calibration like Bruderer et al #######
lnmod <- lapply(msms_runs_irt, function(x){
  rep <- x
  ln <- lm(NormalizedRetentionTime ~ RetentionTime, rep)
  ln})

# Take the mean? of the data?

hp_irt <- lapply(1:length(msms_runs), function(x){
  rep <- msms_runs[[x]]
  mod <- lnmod[[x]]
  normRT <- predict(mod, rep$RetentionTime)
})
  
# Should take the mean first before fitting the linear regression!!!!!

mean_irt <- lapply(1:length(run_id), function(x){
  rep <- subset(msms_runs_irt[[x]], select=c("RetentionTime", "PeptideSequence"))
  setnames(rep, "RetentionTime", run_id[x])
  rep
})


mean_irt <- Reduce(function(x, y) merge(x, y, by = "PeptideSequence", all=TRUE), mean_irt)
mean_irt$mean <- lapply(1:nrow(mean_irt), function(x){
  m <- mean_irt[x, 2:ncol(mean_irt)] %>% as.numeric()
  m <- m %>% mean(na.rm=T)
  m
}) %>% unlist()


mean_irt <- subset(mean_irt, select=c("PeptideSequence", "mean"))
setnames(mean_irt, "mean", "mRT")

mean_irt <- merge(mean_irt, irt, by.x = "PeptideSequence", by.y = "ModifiedPeptideSequence")

col <- sample(colors, size=1)
ggplot(mean_irt, aes(x = NormalizedRetentionTime, y = mRT)) + geom_point(size=1, color = col) + theme_bw(base_size=16)
#ggsave("exo_mean_RT_iRT.png", plot = print(g))

######### Use the mean RT to convert to iRT space ###########

lnmod <- lm(NormalizedRetentionTime ~ mRT, mean_irt)

mean_hpiRT <- lapply(1:nrow(msms_irt), function(x){
  r <- msms_irt[x, 2:ncol(msms_irt)] %>% as.numeric()
  r <- mean(r, na.rm = T)
  r
}) %>% unlist()

hpiRT <- data.frame(Peptide=msms_irt$ModifiedPeptideSequence, mRT=mean_hpiRT, stringsAsFactors = F)
hpiRT

ggplot(hpiRT, aes(x = mRT)) + geom_density(fill = col)
ggsave("/results/exo_mRT_density.png")

############### Use the lnmod to convert the mRT to iRT space ############

newRT <- predict(lnmod, hpiRT)
hpiRT$iRT <- newRT

ggplot(hpiRT, aes(x = iRT)) + geom_density(fill = col)
names(pca_exo.msms_sub) <- raw_files

testRT <- pca_exo.msms_sub[run_id[1]]
testRT <- testRT[[1]]
testRT$ModifiedPeptideSequence <- reformat_mods(testRT$ModifiedPeptideSequence)

testRT <- merge(testRT, hpiRT, by.x = "ModifiedPeptideSequence", by.y = "Peptide")
ggplot(testRT, aes(x = iRT, y = RetentionTime)) + geom_point(size=1)

setnames(hpiRT, c("Peptide", "iRT"), c("ModifiedPeptideSequence", "iRT"))

highPrecisioniRT <- subset(hpiRT, select = c("ModifiedPeptideSequence", "iRT"))
write.table(highPrecisioniRT, file="./data/irt/pca_exo_high_precision_iRT.tsv", row.names = F, quote = F, sep = "\t")
