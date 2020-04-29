library(viridisLite)
library(data.table)
library(ggplot2)
library(tidyr)
library(VennDiagram)

#source("/home/annieha/source/tk_orbitrap_dia/analysis/Utilities.R")
basefolder <- "D:/projects/pca_urine_spectral_lib"
source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")

irt <- list.files(paste0(basefolder, "/data/irt"), pattern = "irt_assaylib", full.names=T) %>% fread()

pca_eps.ev <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "evidence", full.names=T) %>% fread()
pca_eps.msms <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), 
  pattern = "msms.txt", full.names = T) %>% 
  fread()

pca_eps.allPep <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "peptides", full.names=T) %>% fread()

eps_raw_files <- pca_eps.msms$`Raw file` %>% unique()

length(eps_raw_files) # 148 runs
 
irt <- subset(irt, select=c('ModifiedPeptideSequence', 'PrecursorMz', 'NormalizedRetentionTime'))

pca_eps.msms_sub <- lapply(eps_raw_files, function(x){
  rep <- pca_eps.ms[`Raw file` == x]
  rep <- subset(rep, select=c("Raw file","Sequence", "Modified sequence", "Charge",
                              "Retention time", "PEP", "Intensity", "Proteins", "Gene Names"))
  setnames(rep, colnames(rep), c("raw","PeptideSequence", "ModifiedPeptideSequence", "Charge",
                                 "RetentionTime", "pep", "Intensity", "Proteins", "Genes"))
  rep$ModifiedPeptideSequence <- reformat_mods(rep$ModifiedPeptideSequence)
  rep
})

names(pca_eps.msms_sub) <- eps_raw_files
pca_eps.ev <- subset(pca_eps.evidence, select = c("id", "Calibrated retention time", "Intensity"))
pca_eps.ms <- merge(pca_eps.msms, pca_eps.ev, by.x = "Evidence ID", by.y = "id")
################# compare iRT peptides ############################

irt <- irt[!duplicated(ModifiedPeptideSequence, NormalizedRetentionTime)]
irt <- irt[order(NormalizedRetentionTime)]
setnames(irt, "NormalizedRetentionTime", "iRT")

pca_eps.iRT <- lapply(eps_raw_files, function(x){
  rep <- pca_eps.msms_sub[[x]]
  rep <- rep[order(-Intensity)]
  rep <- subset(rep, select=c("ModifiedPeptideSequence", "RetentionTime", "raw", "Intensity"))
  rep <- merge(rep, irt, by = "ModifiedPeptideSequence", all.y = T)
  rep[, run:= x]
  rep
})


names(pca_eps.iRT) <- eps_raw_files

pca_eps.iRT <- do.call('rbind', pca_eps.iRT )
pca_eps.iRT$iRT <- factor(pca_eps.iRT$ModifiedPeptideSequence,levels = irt$ModifiedPeptideSequence)
pca_eps.iRT$cohort <- "dEPS"
ggplot(pca_eps.iRT, aes(x = iRT, y = RetentionTime, color = run, group = iRT)) + 
  theme_bw(base_size = 16) +
  geom_jitter(show.legend = F, width = 0.3, size=0.2, alpha=0.8) +
  geom_boxplot(show.legend = F, width = 0.5, alpha = 0.8) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

ggplot(pca_eps.iRT, aes(x = iRT)) + geom_boxplot(aes(y = RetentionTime)) + 
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 


pca_eps.msms_irt <- lapply(pca_eps.msms_sub, function(x){
  rep <- x
  rep <- subset(rep, select = c("ModifiedPeptideSequence", "RetentionTime", "raw", "Charge"))
  rep <- rep[!duplicated(ModifiedPeptideSequence)]
  run <- unique(x$raw)
  setnames(rep, "RetentionTime", run)
  rep <- merge(rep, irt, by = "ModifiedPeptideSequence")
  rep$NormalizedRetentionTime <- NULL
  rep$PrecursorMz <- NULL
  rep$Charge <- NULL
  rep$raw <- NULL
  rep
})

############### pca analysis on RT and Intensity ###############

irt_pca.eps <- copy(pca_eps.iRT)
irt_pca.eps <- lapply(irt_pca.eps, function(x){
  rep <- x
  rep <- subset(rep, select = c("ModifiedPeptideSequence", "RetentionTime", "Intensity", "run"))
  rep <- melt(rep, id.vars = c("ModifiedPeptideSequence", "run"))
  rep$variable_names <- paste(rep$ModifiedPeptideSequence, rep$variable, sep = ".")
#  rep <- reshape(rep, idvar = "run", timevar = "variable_names",
#                 drop = c("ModifiedPeptideSequence", "variable"),
#                 direction = "wide")
  rep
})

irt_pca.eps <- lapply(irt_pca.eps, function(x){
  t <- x
  t <- reshape(t, idvar = "run", timevar = "variable_names", drop = c("ModifiedPeptideSequence", "variable"),
               direction = "wide")
  t
})

irt_pca.eps <- do.call('rbind', irt_pca.eps)


irt_pca <- prcomp(irt_pca.eps[, 2:23])
library(ggfortify)

autoplot(irt_pca)
########### continue iRT analysis ##################
msms_irt <- Reduce(function(x,y) merge(x,y, by = "ModifiedPeptideSequence", all = T), pca_eps.msms_irt)

stats <- lapply(1:nrow(msms_irt), function(x){
  m <- msms_irt[x, 2:ncol(msms_irt)] %>% as.numeric() %>% mean(na.rm=T)
  sd <- msms_irt[x, 2:ncol(msms_irt)] %>% as.numeric() %>% sd(na.rm=T)
  na <- msms_irt[x, 2:ncol(msms_irt)] %>% is.na() %>% sum()
  c(m, sd, na)
})

stats <- do.call('rbind', stats) %>% as.data.frame(stringsAsFactors=F)
colnames(stats) <- c("mean", "SD", "Missing")
stats$Peptide <- msms_irt$ModifiedPeptideSequence
irt_peptides <- irt$ModifiedPeptideSequence
stats <- merge(stats, irt, by.x = "Peptide", by.y = "ModifiedPeptideSequence")
stats$iRT <- factor(stats$Peptide, levels = irt_peptides)
stats <- stats[order(iRT)]
stats$Peptide <- factor(stats$Peptide)
eps.stats <- stats
eps.stats$cohort <- "direct EPS"
colors <- viridis(10)

ggplot(stats, aes(x = SD)) + geom_density(fill = sample(colors, size=1)) + theme_bw()
ggplot(stats, aes(x = iRT, y = Missing)) + geom_bar(stat = "identity", fill = sample(colors, size=1), width = 0.5) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.y = element_text(size = 16)) 


################### Find common peptides #############################
# Total number of identified peptides:
AllPeptides <- pca_eps.msms$`Modified sequence` %>% reformat_mods() %>% unique()
length(AllPeptides)
# 49,999 unique peptides across ALL

names(pca_eps.msms_sub) <- eps_raw_files
peptide_count <- lapply(names(pca_eps.msms_sub), function(x){
  rep <- pca_eps.msms_sub[[x]]
  precursor <- rep[!duplicated(ModifiedPeptideSequence, Charge)]
  precursor <- precursor %>% nrow()
  peptide <- rep$PeptideSequence %>% unique() %>% length()
  protein <- rep$Proteins %>% unique() %>% length()
  count <- data.frame(raw=x, precursor=precursor, peptide=peptide, protein=protein)
  count
})
peptide_count <- do.call('rbind', peptide_count)


ggplot(peptide_count, aes(x = as.numeric(raw), y = peptide)) + 
  geom_bar(width = 0.5, stat = "identity") + theme_bw() + 
  xlab("Run") + ylab("Precursor Count")

a <- peptideCountOverObservations(pca_eps.allPep)
a$Total <- a$Found + a$Matched
colors <- viridis(10)
df <- data.frame(count = a$Found, category="MS/MS", stringsAsFactors = F)
df1 <- data.frame(count = a$Matched, category = "Matching", stringsAsFactors = F)
df <- rbind(df, df1)

ggplot(a, aes(x = Found)) + geom_bar(stat="count", color = sample(colors, size=1) ) + 
  geom_bar(aes(x = Matched), position = "stack")+ scale_x_reverse()
ggplot(df, aes(x = count, fill = category)) + geom_bar(position="stack") + scale_x_reverse() +
  theme_bw() + scale_y_continuous(breaks = seq(0, 10000, by = 400))

common <- lapply(pca_eps.msms_sub, function(x) unique(x$PeptideSequence))
common <- Reduce(intersect, common)
#109

pca_eps_common <- pca_eps.msms[Sequence %in% common]

pca_eps_common <- pca_eps_common[order(-`Intensity coverage`)]
pca_eps_common_proteins <- pca_eps_common$`Gene Names` %>% table()
pca_eps_common_proteins <- pca_eps_common_proteins %>% as.data.table(rownames=names(pca_eps_common_proteins))
colnames(pca_eps_common_proteins) <- c("Gene Names", "Counts")
pca_eps_common_proteins <- pca_eps_common_proteins[!is.na(`Gene Names`)]
pca_eps_common_proteins$`Gene Names` <- gsub(";(.*)", "_ProteinGroup", pca_eps_common_proteins$`Gene Names`)

ggplot(pca_eps_common_proteins, aes(x = `Gene Names`, y = Counts)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=45, hjust = 1))


################## start alignment ####################################
irt <- irt[!duplicated(ModifiedPeptideSequence)]
colors <- viridis(length(pca_eps.msms_irt))
#pdf("msms_irt.pdf", onefile = T)

pca_eps.msms_irt <- do.call('rbind', pca_eps.msms_irt)

ggplot(pca_eps.msms_irt, aes(x = NormalizedRetentionTime, y = RetentionTime, color = raw)) + geom_line(show.legend = F) +
  geom_point(size=1, show.legend = F)  + theme_bw()


###### Use pairwise comparison to find the best alignment ############

pca_eps.sequence <- lapply(pca_eps.msms_sub, function(x){
  seq <- x$PeptideSequence %>% unique()
  seq
})



pairwise <- matrix(nrow=length(pca_eps.sequence), ncol = length(pca_eps.sequence))

for (i in 1:148) {
  r <- pca_eps.sequence[[i]]
  for(j in 1:148){
    s <- pca_eps.sequence[[j]]
    pair <- intersect(r, s) %>% length()
    pairwise[i,j] <- pair
  }
}


row.names(pairwise) <- eps_raw_files
colnames(pairwise) <- eps_raw_files

pheatmap(pairwise, color = viridis(20))

length(pca_eps.sequence[[105]])
# 12945

ref_count <- data.frame(E156_count = pairwise[105,], run = row.names(pairwise)) %>% as.data.table()
color <- viridis(2, end = 0.5, option = "plasma")
ggplot(ref_count[-105], aes(x = as.integer(run), y = E156_count)) + 
  geom_bar(stat="identity", position = "dodge", width = 0.5, fill = color[2]) + theme_bw(base_size = 16) +
  xlab("Sample") + ylab("Peptide Overlap with E156") + ylim(0, 12000)
ggplot(alignment, aes(x = RetentionTime)) + geom_density(fill = viridis(1), alpha=0.6) + theme_bw(base_size = 16)

E156_irt <- pca_eps.iRT[[105]]
E156_irt <- E156_irt[!duplicated(ModifiedPeptideSequence)]


ggplot(E156_irt, aes(x = RetentionTime , y= NormalizedRetentionTime)) + geom_point() + geom_line()

setnames(E156_irt, "NormalizedRetentionTime", "iRT")
E156_lnmod <- lm(iRT ~ RetentionTime, E156_irt)
E156_irt$predicted <- E156_lnmod$fitted.values
ggplot(E156_irt, aes(x = RetentionTime, y = iRT)) + geom_point() + geom_line() + geom_smooth(method = "lm")
E156_irt$residuals <- E156_lnmod$residuals

# tolerance = 6 minutes
tol <- 6

E156_irt <- E156_irt[abs(residuals) <= 6]
ggplot(E156_irt, aes(x = RetentionTime, y = iRT)) + geom_point() + 
  geom_line() + geom_smooth(method = "lm")

E156_lnmod <- lm(iRT ~ RetentionTime, E156_irt)
E156_irt$NormalizedRetentionTime <- predict(E156_lnmod)
E156_irt <- E156_irt[order(NormalizedRetentionTime)]

E156_Alignment <- pca_eps.msms_sub[["20180613_E156"]]

E156_Alignment$NormalizedRetentionTime <- predict(E156_lnmod, newdata = E156_Alignment)
ggplot(E156_Alignment, aes(x = NormalizedRetentionTime, y = RetentionTime)) + geom_line() + geom_point() +
  theme_bw(base_size = 16)
ggplot(E156_Alignment, aes(x = NormalizedRetentionTime)) + geom_density(fill = viridis(1))


############# Write the alignment ######################

E156_Alignment <- pca_eps.msms[`Raw file` =="20180613_E156"]
E156_Alignment <- subset(E156_Alignment, select = c("Raw file", "Sequence", "Modified sequence", "Proteins", 
                                                    "Gene Names", "Charge", "m/z", "Mass", "Retention time",
                                                    "PEP", "Intensities", "Masses", "Evidence ID", 'Reverse'))

E156_Alignment <- E156_Alignment[Reverse != "+"]
E156_Alignment <- merge(E156_Alignment, pca_eps.ev, by.x = "Evidence ID", by.y = "id")

setnames(E156_Alignment, c("Raw file", "Sequence", "Modified sequence", "Gene Names", "m/z", "Retention time"),
         c("raw", "PeptideSequence", "ModifiedPeptideSequence","Genes", "m/z", "RetentionTime"))

setnames(E156_Alignment, "Evidence ID", "id")

E156_Alignment$ModifiedPeptideSequence <- reformat_mods(E156_Alignment$ModifiedPeptideSequence)
#E156_Alignment <- E156_Alignment[order(-Intensity)]
E156_Alignment <- E156_Alignment[order(PEP)]
E156_Alignment <- E156_Alignment[!duplicated(ModifiedPeptideSequence, Charge)]
E156_Alignment$NormalizedRetentionTime <- predict(E156_lnmod, newdata = E156_Alignment)
E156_Alignment$`Calibrated retention time` <- NULL

write.table(E156_Alignment, file = "D:/projects/pca_urine_spectral_lib/data/irt/direct_EPS_pairwise_E156_hpirt.tsv",
            sep = "\t", row.names = F, quote = F)

E156 <- pairwise[, "E156"]

E156 <- as.data.table(E156)
ggplot(E156, aes(x = 1:148, y = E156)) + geom_bar(stat = "identity", fill = viridis(10)[4]) + 
  theme_bw(base_size = 11) + xlab("Runs") + ylab("Peptide Intersection Count with E156")
min(E156)
# 1602

# Use this table to align for the rest
# Use lowess alignment and remove outliers that are over 7 minutes
pdf("RT_alignment_E156_nonlinear.pdf", onefile = T)

E156_Alignment$id %>% unique() %>% length()
# 13392
