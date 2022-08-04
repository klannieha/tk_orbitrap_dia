########### this script is to do statistical tests ##############
library(data.table)
library(ggplot2)
library(tidyverse)
library(aLFQ)
library(SWATH2stats)
library(ggpubr)
library(ggrepel)
library(ggpmisc)
library(matrixStats)
library(BoutrosLab.plotting.general)
library(matrixStats)
library(ggrepel)

source("D:/source/tk_orbitrap_dia/analysis/Utilities.R")
######################## Load data #############################
load("D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/diapeptides_quant_annotated.RData")
load("D:/projects/pca_urine_spectral_lib/data/openswath/20220328_dia_osw_filtered.RData", environment())

tissue.pep <- read.table("D:/projects/pca_urine_spectral_lib/data/openswath/CPCG_mpMRI_all_peptides.csv", 
                         header = T, sep = ",")
load("D:/projects/pca_urine_spectral_lib/data/clinical/uEPS1_patients.RData")
# data is annotated & is normalized

# filter out peptides based on peptide
# 1. charges 2 & 3
# 2. Observed in > 25% of runs
# 3. No miscleavages
# 4. AA length 7-25
# 5. present in tissue

head(dia.peptides)
dia.peptides.count <- dia.peptides %>%
  group_by(peptide_id) %>%
  summarise(Count = n()) %>% 
  as.data.table()

# filter data
dia.peptideFilter <- dia.peptides %>%
  filter(precursor_charge %in% c(2,3)) %>% # 48922
  filter(peptide_id %in% dia.peptides.count[Count >= 189*0.25]$peptide_id) %>% #15849
  mutate(Subseq = substr(peptide_sequence, 1, nchar(peptide_sequence)-1)) %>%
  filter(!grepl("K|R", Subseq)) %>% # 11587
  filter(nchar(peptide_sequence) >= 7 & nchar(peptide_sequence) <25) %>% # 11346
  filter(peptide_sequence %in% tissue.pep$sequence) %>% #distinct(peptide_sequence) # 5947
  as.data.table()

#################### sum up peptide quant ########################3
head(dia.peptides)

head(dia.peptideFilter)

dia.peptides.wide <- dcast(dia.peptides, peptide_sequence ~ SampleID, sum,
                           value.var = "Intensity.medCor", fill = NA)

write.table(dia.peptides.wide, file = "D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/20220406_diaPeptides_uEPSlib_Normalized_peptides.tsv",
            sep = "\t", row.names = F, quote = F)

head(dia.peptides.wide)

pep.de <- melt(dia.peptides.wide, id.vars = "peptide_sequence", value.name = "Intensity", variable.name = "SampleID") %>% as.data.table()
head(pep.de)

pep.de <- pep.de[!is.na(Intensity)]

pep.de <- merge(pep.de, dia.peptideFilter[, c("SampleID", "ID", "protein_id", "peptide_sequence",
                                              "Run", "Batch", "Age-DX", "PrePSA", "ClinStageT", 
                                              "cISUP", "PathStageT", "pISUP", "cGrade")],
                by = c("SampleID", "peptide_sequence"))

pep.de %>% distinct(peptide_sequence) # 5947
dia.peptideFilter %>% distinct(peptide_sequence)

head(pep.de)
pep.de$log2Intensity <- log2(pep.de$Intensity)

pep.de <- pep.de %>%
  mutate(pGrade = ifelse(pISUP == 1, "1", "2+")) %>%
  as.data.table()

#################### Statistical tests ###########################

head(dia.peptideFilter)
head(pep.de)

FC <- pep.de %>%
  group_by(cGrade, peptide_sequence, protein_id) %>% # CHANGE categories
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "protein_id"),
              names_from = "cGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>%
  left_join(peptideMap %>% distinct(protein_id, peptide_sequence, GeneName),
            by = c("protein_id", "peptide_sequence"))


FC.pISUP <- pep.de %>%
  group_by(pGrade, peptide_sequence, protein_id) %>% # CHANGE categories
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "protein_id"),
              names_from = "pGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>%
  left_join(peptideMap %>% distinct(protein_id, peptide_sequence, GeneName),
            by = c("protein_id", "peptide_sequence"))
FC.pISUP
FC.noNorm <- dia.peptideFilter %>%
  group_by(cGrade, peptide_sequence, peptide_id, protein_id) %>%
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "peptide_id", "protein_id"),
              names_from = "cGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>%
  left_join(peptideMap, by = c("protein_id", "peptide_id", "peptide_sequence"))

FC %>% ggplot(aes(x = log2FC )) +
  geom_histogram(fill = NA, color = 'black', bins = 100) +
  theme_classic(base_size = 16)

wt <- lapply(seq_along(FC$peptide_sequence), function(x){
  p <- FC$peptide_sequence[x]
  t <- Htest(p, peptide_sequence, "1", "2+", cGrade, pep.de, intensity_col = log2Intensity) # change it!!
  t
})

FC$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC$p.adj.value <- p.adjust(FC$p.value, method = "fdr")

FC.pISUP$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.pISUP$p.adj.value <- p.adjust(FC.pISUP$p.value, method = "fdr")

FC %>% ggplot(aes( x = pISUP.pvalue)) +
  geom_histogram(bins = 100) +
  ylab("Peptides (n)") + xlab(expression(p*"-"*value[pISUP])) +
  plot_theme() 
FC
FC.filter %>% # CHANGE dataset!!!
  mutate(LABEL = paste(peptide_sequence, GeneName, sep = "*"),
         Sig = ifelse(pISUP.p.adj.value < 0.05, TRUE, FALSE)) %>%
  mutate(LABEL = ifelse(Sig, LABEL, NA),
         Enriched = ifelse(Sig & log2FC >0 , "UP", NA)) %>%
  mutate(Enriched = ifelse(Sig & log2FC < 0, "DOWN", Enriched)) %>%
  ggplot(aes(x = log2FC, y = -log10(pISUP.p.value))) +
  geom_point(aes(color = Enriched), show.legend = F) +
  #geom_label_repel(aes(label = LABEL), max.overlaps = 30) +
  scale_color_manual(breaks = c("UP", "DOWN"), values = c("red", "blue")) +
  scale_x_continuous(limits = c(-4, 6)) +
  xlab(expression(log[2]*"FC (cISUP"*('2+'/'1'*")"))) +
  ylab(expression("-log"[10]*"p-value")) +
  plot_theme()

FC %>% filter(p.adj.value < 0.05) %>% distinct(peptide_sequence) %>% nrow()# 23 targets
# pISUP statistically sigfig: 9 targets

FC %>%
  full_join(FC.pISUP, by = c("peptide_sequence", "protein_id", "GeneName"),
            suffix = c(".cISUP", ".pISUP")) %>%
  mutate(Significance = case_when((p.adj.value.cISUP < 0.05 & p.adj.value.pISUP < 0.05) ~ "Both",
                                  (p.adj.value.cISUP < 0.05 & p.adj.value.pISUP > 0.05) ~ "cISUP",
                                  (p.adj.value.cISUP > 0.05 & p.adj.value.pISUP < 0.05) ~ "pISUP")) %>% #View()
  mutate(Label = paste(peptide_sequence, GeneName, sep = " * ")) %>%
  mutate(Label = ifelse(Significance == "Both", Label, NA)) %>%
  ggplot(aes(x = log2FC.cISUP, y = log2FC.pISUP)) +
  geom_point(aes(color = Significance, alpha = Significance)) +
  stat_cor(cor.coef.name = 'rho') +
  geom_label_repel(aes(label = Label), max.overlaps = 25, box.padding = 1) +
  scale_color_manual(values = c("green", "dodgerblue2", "darkorange1"), 
                     breaks = c("Both", "cISUP", "pISUP")) +
  scale_alpha_manual(breaks = c("Both", "cISUP", "pISUP"), na.value = 0.2, values = c(1,1,1)) +
  scale_x_continuous(limits = c(-2, 6)) +
  scale_y_continuous(limits = c(-2,6)) +
  ggtitle(expression(beta*' coefficient')) +
  xlab("cISUP") + ylab("pISUP") +
  plot_theme()


# Linear model
# Combine the gg 4 + 5 -> linear model
# Use Kruskal-Wallis between each groups
# 1. Regular linear model (cISUP / pISUP ~ log2Intensity)
# 2. Combined 4+5 linear model

head(pep.de)

wt <- lapply(seq_along(FC.pISUP$peptide_sequence), function(x){ # change dataset!!
  p <- FC.pISUP$peptide_sequence[x]
  lnmod <- lm(log2Intensity ~ pISUP, pep.de[peptide_sequence == p])
  lnmod <- summary(lnmod)
})

FC$beta1 <- lapply(wt, function(x) x$coefficients[2,1]) %>% unlist()
FC$p.value.lm <- lapply(wt, function(x) x$coefficients[2,4]) %>% unlist()
FC$p.adj.lm <- p.adjust(FC$p.value.lm , method = "fdr")

FC.pISUP$beta1 <- lapply(wt, function(x) x$coefficients[2,1]) %>% unlist()
FC.pISUP$p.value.lm <- lapply(wt, function(x) x$coefficients[2,4]) %>% unlist()
FC.pISUP$p.adj.lm <- p.adjust(FC.pISUP$p.value.lm , method = "fdr")


FC%>% ggplot(aes(x = p.adj.lm.comb)) + 
  geom_histogram(bins = 100) +
  xlab(expression("FDR"[lm])) +
  ylab("Peptides (n)") +
  plot_theme()

FC %>% filter(p.adj.lm.comb < 0.05) %>% distinct(peptide_sequence) # 25 targets

FC %>%
  mutate(Sig = ifelse(p.adj.lm.comb < 0.05, T, F)) %>%
  mutate(Enriched = case_when((Sig & beta1 > 0) ~ "Up", (Sig & beta1 < 0) ~ "Down")) %>%
  mutate(Label = ifelse(Sig, paste(peptide_sequence, GeneName, sep = " * "), NA)) %>%
  ggplot(aes(x = beta1.comb, y = -log10(p.value.lm.comb))) +
  geom_point(aes(color = Enriched), show.legend = F) +
  geom_label_repel(aes(label = Label), max.overlaps = 40, box.padding = 2) +
  scale_x_continuous(limits = c(-3, 4)) +
  scale_color_manual(values = c("red", "blue"), breaks = c("Up", "Down")) +
  xlab(expression(beta[1])) + ylab(expression("-log"[10]*"p-value")) +
  plot_theme()

pep.de %>%
  filter(peptide_sequence %in% FC[p.adj.lm.comb < 0.05 & beta1.comb > 0 ]$peptide_sequence) %>%
  inner_join(FC, by = "peptide_sequence") %>%
  mutate(Label = paste(peptide_sequence, GeneName, sep = " * "),
         cGradeCombined = as.integer(cGradeCombined)) %>%
  ggplot() +
  geom_point(aes(x = cGradeCombined, y = log2Intensity)) + 
  geom_boxplot(aes(x = cGradeCombined, y = log2Intensity, group = cGradeCombined),alpha = 0.6, outlier.shape = NA) +
  geom_smooth(aes(x = cGradeCombined, y = log2Intensity), method = "lm", se = F, color = "darkblue") +
  geom_text(data = FC[p.adj.lm.comb < 0.05 & beta1.comb > 0] %>% 
              mutate(Label = paste(peptide_sequence, GeneName, sep = " * "),
                     b.value = round(beta1.comb, digits = 4),
                     p.adj.value = round(p.adj.lm.comb, digits = 5)),
            x = 3, y = 42, hjust = 0.1,
            aes(label = paste('beta[1]==', b.value, '~~~', "FDR==", p.adj.value)), parse = T) +
  facet_wrap(~Label) +
  theme_classic(base_size = 16) +
  ylab(expression(log[2]*"(Intensity)"))

FC %>% filter(p.adj.lm < 0.05) %>% distinct(peptide_sequence) # 18 targets # 4 targets in pISUP

## combine the isup 4+5s

head(pep.de)
pep.de$pISUP <- as.integer(pep.de$pISUP)
pep.de <- pep.de %>%
  mutate(cGradeCombined = ifelse(cISUP == 4 | cISUP == 5, 4, cISUP ),
         pGradeCombined = ifelse(pISUP == 4 | pISUP == 5, 4, pISUP)) %>%
  as.data.table()
# redo linear models
wt <- lapply(seq_along(FC.pISUP$peptide_sequence), function(x){ # change dataset!!
  p <- FC.pISUP$peptide_sequence[x]
  lnmod <- lm(log2Intensity ~ pGradeCombined, pep.de[peptide_sequence == p])
  lnmod <- summary(lnmod)
})

FC$beta1.comb <- lapply(wt, function(x) x$coefficients[2,1]) %>% unlist()
FC$p.value.lm.comb <- lapply(wt, function(x) x$coefficients[2,4]) %>% unlist()
FC$p.adj.lm.comb <- p.adjust(FC$p.value.lm.comb, method = "fdr")

FC.pISUP$beta1.comb <- lapply(wt, function(x) x$coefficients[2,1]) %>% unlist()
FC.pISUP$p.value.lm.comb <- lapply(wt, function(x) x$coefficients[2,4]) %>% unlist()
FC.pISUP$p.adj.lm.comb <- p.adjust(FC.pISUP$p.value.lm.comb, method = "fdr")


# Further filter the peptides to intersection between dEPS/comb/uEPS library searched
# Test out Kruskal Wallis

p <- "QHMDSDSSPSSSSTYCNQMMR"
p
kruskal.test(log2Intensity~cISUP, pep.de[peptide_sequence == p])
pairwise.wilcox.test(pep.de[peptide_sequence == p]$log2Intensity, pep.de[peptide_sequence == p]$cISUP,
                     p.adjust.method = "BH")

pep.de %>%
  filter(peptide_sequence == p) %>%
  distinct(`Age-DX`, cISUP, log2Intensity) %>%
  ggplot(aes(x = `Age-DX`, y = log2Intensity, color = factor(cISUP))) +
  geom_point(fill = NA) + theme_classic()

pep.de %>% filter(peptide_sequence == p) %>%
  ggplot(aes(x = cISUP, y = log2Intensity, group = cISUP)) +
  geom_boxplot()

### Try out KW test
wt <- lapply(seq_along(FC.pISUP$peptide_sequence), function(x){ # CHANGE dataset!!
  p <- FC.pISUP$peptide_sequence[x]
  t <- kruskal.test(log2Intensity ~ pISUP, pep.de[peptide_sequence == p])
  t
})

FC$p.value.kw <- lapply(wt, function(x) x$p.value) %>% unlist()
FC$p.adj.kw <- p.adjust(FC$p.value.kw, method = 'fdr')

FC.pISUP$p.value.kw <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.pISUP$p.adj.kw <- p.adjust(FC.pISUP$p.value.kw, method = 'fdr')

remove(wt)
FC %>% ggplot(aes(x = p.adj.kw)) + geom_histogram()
FC %>% filter(p.adj.kw < 0.05) %>% ungroup %>% summarise(n = n()) # 34 peptides
FC <- as.data.table(FC)
FC %>% filter(p.adj.kw < 0.05 & log2FC < 0) %>% ungroup %>% summarise(n = n()) # 20 up, 14 down
FC.pISUP %>% filter(p.adj.value < 0.05 & p.adj.lm < 0.05 & p.adj.kw < 0.05)


pep.de %>%
  filter(peptide_sequence %in% FC[p.adj.lm < 0.05]$peptide_sequence) %>%
  inner_join(peptideMap, by = 'peptide_sequence') %>%
  mutate(Label = paste(peptide_sequence, GeneName, sep = " * ")) %>%
  ggplot(aes(x = cISUP, y = log2Intensity)) +
  geom_point(alpha = 0.3) +
  geom_boxplot(aes(group = cISUP), alpha = 0.6) +
  stat_smooth(method = 'lm', se = F) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho') +
  facet_wrap(~Label) +
  theme_classic()

FC <- as.data.table(FC)
FC.pISUP <- as.data.table(FC.pISUP)

create_venn(list(FC[p.adj.value < 0.05]$peptide_sequence, 
                 FC[p.adj.lm < 0.05]$peptide_sequence,
                 FC[p.adj.kw < 0.05]$peptide_sequence),
            category_names = c("uTest", "lnmod", "KW"),
            category_colours = brewer.pal(3, "Set2"))



################# Filter more #########################

osw <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/",
                  pattern = "wide_intensity.tsv", full.names = T)

osw
dEPS <- fread(osw[2])
comb <- fread(osw[1])
head(dEPS)
head(comb)

create_venn(list(dEPS$peptide_sequence,
                 comb$peptide_sequence, 
                 unique(dia.peptideFilter$peptide_sequence)),
            category_names = c('dEPS', "comb", "uEPS.filtered"), 
            category_colours = dataset_colours[c("dEPS1", "CPCG", "uEPS1")])


# filter further: take only filtered ones in all 3 libraries (highest confidence)

peptides <- Reduce(intersect, list(dEPS$peptide_sequence, comb$peptide_sequence, 
                                   unique(dia.peptideFilter$peptide_sequence)))

pep.de.filter <- pep.de[peptide_sequence %in% peptides]

FC.filter <- pep.de.filter %>%
  group_by(cGrade, peptide_sequence, protein_id) %>% # CHANGE categories
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "protein_id"),
              names_from = "cGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>%
  left_join(peptideMap %>% distinct(protein_id, peptide_sequence, GeneName),
            by = c("protein_id", "peptide_sequence"))

FC.pISUP.filter <- pep.de.filter %>%
  group_by(pGrade, peptide_sequence, protein_id) %>% # CHANGE categories
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "protein_id"),
              names_from = "pGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>%
  left_join(peptideMap %>% distinct(protein_id, peptide_sequence, GeneName),
            by = c("protein_id", "peptide_sequence"))

wt <- lapply(seq_along(FC.filter$peptide_sequence), function(x){
  p <- FC.pISUP.filter$peptide_sequence[x]
  t <- Htest(p, peptide_sequence, "1", "2+", pGrade, pep.de, intensity_col = log2Intensity) # change it!!
  t
})

FC.filter$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.filter$p.adj.value <- p.adjust(FC.filter$p.value , method = 'fdr')
FC.filter %>% filter(p.adj.value < 0.05) %>% summarise(Count = n()) # 21 targets

FC.pISUP.filter$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.pISUP.filter$p.adj.value <- p.adjust(FC.pISUP.filter$p.value, method = "fdr")

plot_volcano(FC.filter, label = F, `p.value`, `p.adj.value`, `log2FC`) +
  scale_x_continuous(limits = c(-4, 6)) +
  xlab(expression(log[2]*"FC(cISUP 2+/1)")) +
  ylab(expression(log[10]*'(p-value)'))


dia.peptides.count
dia.sampleCount <- dia.peptides %>%
  group_by(SampleID) %>%
  distinct(peptide_sequence) %>%
  summarise(PeptideCount = n())
dia.sampleCount


dia.sampleInt <- dia.peptides %>%
  group_by(SampleID) %>%
  summarise(MedianInt = median(log2Intensity))


dia.sampleCount <- merge(dia.sampleCount, dia.sampleInt, by = "SampleID")
dia.sampleCount %>% ggplot(aes(x = PeptideCount, y = MedianInt)) + geom_point()
