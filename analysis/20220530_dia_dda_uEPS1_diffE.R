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
######################## Load data #############################
load("D:/projects/pca_urine_spectral_lib/data/openswath/20211014_TPPLibs/diapeptides_quant_annotated.RData")
tissue.pep <- read.table("D:/projects/pca_urine_spectral_lib/data/openswath/CPCG_mpMRI_all_peptides.csv", 
                         header = T, sep = ",")
load("D:/projects/pca_urine_spectral_lib/data/clinical/uEPS1_patients.RData")
# data is annotated & is normalized

dda.peptide <- read.table('D:/projects/pca_urine_spectral_lib/data/pca_dda/pca_postDRE_urine/uEPS1_Normalized/loess_uEPS1/loess_uEPS1/normInt_optimalSpan.txt',
                          sep = "\t", header = T)

mq.peptide <- read.table("D:/projects/pca_urine_spectral_lib/data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og/peptides.txt",
                         sep = "\t", header = T)

mq.peptide <- subset(mq.peptide,
                     select = c("Sequence", "Length", "Missed.cleavages", "Charges", "Gene.names", "Protein.names"))

mq.peptide <- as.data.table(mq.peptide)
# check data

dda.peptide %>% distinct(peptide) %>% summarise(Count = n()) # 52042

create_venn(list(unique(mq.peptide$Sequence),
                 dda.peptide$peptide),
            category_names = c("MQ", "DDA"),
            category_colours = c("#EAFF80", "#CCF7D2"))


dda.peptide %>%
  pivot_longer(cols = -"peptide", names_to = "Sample", values_to = "Intensity", values_drop_na = T) %>%
  ggplot(aes(x = Sample, y = Intensity)) +
  geom_boxplot()

######################### Filter peptides ########################################
# filter out peptides based on peptide
# 1. charges 2 & 3
# 2. Observed in > 25% of runs
# 3. No miscleavages
# 4. AA length 7-25
# 5. present in tissue
dda.peptides.long <- dda.peptide %>%
  pivot_longer(cols = -"peptide", names_to = "Sample", values_to = "Intensity", values_drop_na = T) %>% 
  mutate(Sample = str_extract(Sample, "UP[0-9]+")) %>%
  as.data.table()
dda.peptides.count <- dda.peptides.long %>%
  group_by(peptide) %>% summarise(Count = n()) %>% 
  as.data.table()
head(dda.peptides.long)

head(dda.peptides.count)
dda.peptides.count %>% 
  ggplot(aes(x = Count)) +
  geom_histogram(binwidth = 5)
dia.peptides.count %>%
  ggplot(aes(x = Count)) + geom_histogram(binwidth = 5)

dda.peptides.long %>%
  filter(Sample %in% unique(dia.peptides$SampleID)) %>%
  distinct(Sample)

dda.peptideFilter <- dda.peptides.long %>%
  filter(peptide %in% mq.peptide[grepl("2|3", mq.peptide$Charges)]$Sequence) %>% # 43217
  filter(peptide %in% dda.peptides.count[Count >= 199*0.25]$peptide) %>% #17567
  filter(peptide %in% mq.peptide[Missed.cleavages == 0]$Sequence) %>% # 10935
  filter(nchar(peptide) >= 7 & nchar(peptide) <25) %>% # 10380
  filter(peptide %in% tissue.pep$sequence) %>%  # 8611
  inner_join(mq.peptide[, c("Sequence", "Gene.names", "Protein.names")], 
             by = c("peptide" = "Sequence")) %>%
  as.data.table()

# check data
ggplot(dda.peptideFilter, aes(x = Sample, y = Intensity)) +
  geom_boxplot()

peptides.lst <- list(dda = unique(dda.peptideFilter$peptide),
                     dia = unique(pep.de$peptide_sequence))
fit <- euler(peptides.lst, control = list(extraopt = T))

plot(fit,
     fills = c("#FAE5A1", "#99C19A"),
     #col = c("gold", "dodgerblue", "#99d8c9"),
     edges = F,
     quantities = list(fontsize = 18), adjust_labels = T, lwd = 4, 
     legend = list(fontsize = 30, side = "bottom", nrow = 1, ncol = 4))


######################### Statistical tests #################################
head(dda.peptideFilter)

dda.pep.de <- merge(dda.peptideFilter, patients, by.x = "Sample", by.y = "Sample ID", all.x = T)

head(dda.pep.de)

FC.dda <- dda.pep.de %>%
  group_by(cGrade, peptide, Protein.names, Gene.names) %>% # CHANGE categories
  summarise(mean = mean(Intensity)) %>%
  pivot_wider(id_cols = c("peptide", "Gene.names", "Protein.names"),
              names_from = "cGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>% as.data.table()

FC.dda %>% ggplot(aes(x = log2FC)) + geom_histogram()

# check with dia
FC.dda %>%
  full_join(FC, by = c("peptide" = "peptide_sequence"), suffix = c(".DDA", ".DIA")) %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_hex(binwidth = 0.1) +
  stat_cor(cor.coef.name = "rho") + 
  expand_limits(x = c(-2, 6), y = c(-2, 6)) +
  scale_fill_distiller(palette = "RdPu", direction = -1) +
  plot_theme() +
  xlab(expression(log[2]*"FC"[DDA])) +
  ylab(expression(log[2]*"FC"[DIA])) +
  ggtitle("Effect size (cISUP2+/1)")

wt <- lapply(seq_along(FC.dda$peptide), function(x){ # change dataset!!
  p <- FC.dda$peptide[x]
  t <- Htest(p, peptide, "1", "2+", cGrade, dda.pep.de, intensity_col = Intensity) # change it!!
  t
})

FC.dda$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.dda$p.adj.value <- p.adjust(FC.dda$p.value, method = "fdr")
head(FC.dda)

FC.dda %>% filter(p.adj.value < 0.05) # 16

FC.dda %>% ggplot(aes(x = p.adj.value)) + geom_histogram(bins = 100) +
  plot_theme() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1500), breaks = seq(0, 1500, 200)) +
  ylab("Peptides (n)") +
  xlab("FDR")

FC.dda %>%
  mutate(Sig = ifelse(p.adj.value < 0.05, TRUE, FALSE)) %>%
  mutate(Enrich = case_when(Sig & log2FC < 0 ~ "DOWN",
                            Sig & log2FC > 0 ~ "UP")) %>%
  ggplot(aes(x = log2FC, y = -log10(p.value), color = Enrich)) +
  geom_point(show.legend = F) +
  scale_color_manual(breaks = c("UP", "DOWN"), values = c("red", "blue")) +
  plot_theme() +
  xlab(expression(log[2]*"FC (cISUP 2+/1)")) +
  ylab(expression("-log"[10]*"(p-value)"))

# compare signifance with dia
peptides.lst <- list(DIA = FC[p.adj.value < 0.2]$peptide_sequence,
                     DDA = FC.dda[p.adj.value < 0.2]$peptide)

fit <- euler(peptides.lst, control = list(extraopt = T))

plot(fit,
     fills = c("#FAE5A1", "#99C19A"),
     #col = c("gold", "dodgerblue", "#99d8c9"),
     edges = F,
     quantities = list(fontsize = 18), adjust_labels = T, lwd = 4, 
     legend = list(fontsize = 30, side = "bottom", nrow = 1, ncol = 4))

Reduce(intersect, peptides.lst)
 # 1 overlap TF
FC.dda %>%
  mutate(Sig = ifelse(p.adj.value < 0.05, TRUE, FALSE)) %>%
  full_join(FC %>% 
              mutate(Sig = ifelse(p.adj.value < 0.05, T, F)),
            by = c("peptide" = "peptide_sequence"), suffix = c(".DDA", ".DIA")) %>%
  mutate(SigFull = case_when(Sig.DDA & Sig.DIA ~ "Both",
                             Sig.DDA & !Sig.DIA ~ "DDA", 
                             Sig.DIA & !Sig.DDA ~ "DIA")) %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_point(aes(color = SigFull, alpha = SigFull)) +
  stat_cor(cor.coef.name = "rho") + 
  expand_limits(x = c(-2, 6), y = c(-2, 6)) +
  scale_alpha_manual(breaks = c("DDA", "DIA", "Both"), values = c(0.9, 0.9, 0.9), na.value = 0.2,
                     name = "Significance") +
  scale_color_manual(breaks = c("DDA", "DIA", "Both"), values = c("#FAE5A1", "#99C19A", "red"),
                     name = 'Significance') +
  plot_theme() +
  xlab(expression(log[2]*"FC"[DDA])) +
  ylab(expression(log[2]*"FC"[DIA])) +
  ggtitle("Effect size (cISUP2+/1)")


# try linear model

wt <- lapply(seq_along(FC.dda$peptide), function(x){ # change dataset!!
  p <- FC.dda$peptide[x]
  lnmod <- lm(Intensity ~ cISUP, dda.pep.de[peptide == p])
  lnmod <- summary(lnmod)
})

FC.dda <- as.data.table(FC.dda)

FC.dda$beta1 <- lapply(wt, function(x) x$coefficients[2,1]) %>% unlist()
FC.dda$p.value.lm <- lapply(wt, function(x) x$coefficients[2,4]) %>% unlist()
FC.dda$p.adj.lm <- p.adjust(FC.dda$p.value.lm, method = "fdr")

FC.dda %>%
  ggplot(aes(x = p.adj.lm)) + geom_histogram(bins = 100)

FC.dda %>%
  ggplot(aes(x = p.adj.lm)) +
  geom_histogram(bins = 100) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1000)) +
  plot_theme() +
  ylab("Peptides (n)") + xlab("FDR")


FC.dda %>%
  mutate(Sig = ifelse(p.adj.lm < 0.05, TRUE, FALSE)) %>%
  mutate(Enrich = case_when(Sig & beta1 > 0 ~ "UP",
                            Sig & beta1 < 0 ~ "DOWN")) %>%
  mutate(Label = ifelse(Sig, paste(peptide, Gene.names, sep = " * "), NA )) %>%
  ggplot(aes(x = beta1, y = -log10(p.value.lm))) +
  geom_point(aes(color = Enrich), show.legend = F) +
  geom_label_repel(max.overlaps = 20, box.padding = 1, aes(label = Label)) +
  scale_color_manual(breaks = c("UP", "DOWN"), values = c("red", "blue")) + plot_theme() +
  xlab(expression(beta[1])) + ylab(expression("-log"[2]*"(p-value)"))


FC.dda %>% 
  mutate(Sig = ifelse(p.adj.lm < 0.05, T, F)) %>%
  full_join(FC %>% 
                       mutate(Sig = ifelse(p.adj.lm < 0.05, T, F)),
                     by = c("peptide" = "peptide_sequence"), suffix = c(".DDA", ".DIA")) %>%
  mutate(SigFull = case_when(Sig.DDA & Sig.DIA ~ "Both",
                             Sig.DDA & !Sig.DIA ~ "DDA", 
                             Sig.DIA & !Sig.DDA ~ "DIA")) %>%
  ggplot(aes(x = beta1.DDA, y = beta1.DIA)) +
  geom_point(aes(color = SigFull, alpha = SigFull)) +
  stat_cor(cor.coef.name = "rho") + 
  expand_limits(x = c(-2, 6), y = c(-2, 6)) +
  scale_alpha_manual(breaks = c("DDA", "DIA", "Both"), values = c(0.9, 0.9, 0.9), na.value = 0.2,
                     name = "Significance") +
  scale_color_manual(breaks = c("DDA", "DIA", "Both"), values = c("#FAE5A1", "#99C19A", "red"),
                     name = 'Significance') +
  plot_theme() +
  xlab(expression(beta[1*"DDA"])) +
  ylab(expression(beta[1*"DIA"]))


# dot map

dotmap <- FC.dda %>%
  filter(p.adj.lm < 0.2) %>%
  select(peptide, beta1, p.adj.lm, Gene.names) %>%
  full_join(FC %>% filter(p.adj.lm < 0.2) %>%
              mutate(Label = paste(peptide_sequence, GeneName, sep = " * ")) %>%
              select(peptide_sequence, beta1, p.adj.lm, GeneName), 
            by = c("peptide" = "peptide_sequence", "Gene.names" = "GeneName"),
            suffix = c(".DDA", ".DIA")) %>%
  mutate(Label = paste(peptide, Gene.names, sep = " * ")) %>%
  as.data.frame()
dotmap 
spot.size.function <- function(x) { 0.1 + (2 * abs(x)); }
spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1];
  colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2];
  return(colours);
}

dotmap <- as.data.table(dotmap)

head(peptideList.filter)
# NEED TO SCALE TO NORMALIZED EFFECT SIZE -1 to 1!!!!
dotmap.subset <- dotmap %>% filter(p.adj.lm.DDA < 0.05) %>% arrange(p.adj.lm.DDA) %>% top_n(10)

# normalize log2FC to z score -1 to 1

peptideList.filter <- peptideList.filter %>%
  mutate(Zscore.DDA = (log2FC.DDA - mean(log2FC.DDA, na.rm = T))/sd(log2FC.DDA, na.rm = T),
         Zscore.DIA = (log2FC.DIA - mean(log2FC.DIA, na.rm = T))/sd(log2FC.DIA, na.rm = T)) %>%
  as.data.table()


create.dotmap(peptideList.filter[p.adj.value.DIA < 0.05, c("Zscore.DDA", "Zscore.DIA")],
              yaxis.cex = 1.5,
              xaxis.cex = 1.5,
              na.pch = 1,
              na.spot.size = 0,
              yaxis.lab = paste(peptideList.filter[p.adj.value.DIA < 0.05]$peptide_sequence,
                                peptideList.filter[p.adj.value.DIA < 0.05]$GeneName, sep = " * "),
              xaxis.lab = c("DDA", "DIA"), 
              spot.size.function = spot.size.function,
              spot.colour.function = spot.colour.function,
              key = list(
                space = "right",
                points = list(
                  cex = spot.size.function(seq(-1, 1, 0.2)),
                  col = spot.colour.function(seq(-1, 1,0.2)),
                  pch = 19
                ),
                text = list(
                  lab = as.character(seq(-1, 1, 0.2)),
                  cex = 1.5,
                  adj = 1.0,
                  fontface = "bold"
                )
              ),
              # control spacing at top of key
              key.top = 1,
              pch = 21,
              pch.border.col = "white",
              # add the background
              bg.data = peptideList.filter[p.adj.value.DIA < 0.05, c("p.adj.value.DDA", "p.adj.value.DIA")],
              # add a colourkey
              colourkey = TRUE,
              # set colour scheme for background data
              colour.scheme = c("black", "white"),
              # make bg colour scheme a discrete colour scheme, with breaks at these places
              at = rev(seq(0, 0.1, 0.01))
              )


################ Trends ############################



head(dotmap)
all.FC <- merge(FC, FC.dda, by.x = c("peptide_sequence", "GeneName"), by.y = c("peptide", "Gene.names"),
                all = T, suffixes = c(".DIA", ".DDA"))
peptideList <- all.FC %>%
  filter(p.adj.value.DDA < 0.2 | p.adj.value.DIA < 0.2)

head(peptideList)

dia.countPeptidesList <- pep.de %>%
  filter(peptide_sequence %in% peptideList$peptide_sequence) %>%
  group_by(peptide_sequence) %>%
  distinct(SampleID) %>%
  summarise(Observations = n()) %>%
  mutate(ObservationPercent = Observations/189) %>%
  mutate(Method = "DIA")

dda.countPeptidesList <- dda.pep.de %>%
  filter(peptide %in% peptideList$peptide_sequence) %>%
  group_by(peptide) %>%
  distinct(Sample) %>%
  summarise(Observations = n()) %>%
  mutate(ObservationPercent = Observations/199) %>%
  mutate(Method = "DDA")
setnames(dda.countPeptidesList, "peptide", "peptide_sequence")

head(peptideList)

peptideList <- merge(peptideList, dia.countPeptidesList %>% select(-Method), by = "peptide_sequence")
setnames(peptideList, c("Observations", "ObservationPercent"), c("Observations.DIA", "ObservationPercent.DIA"))
peptideList <- merge(peptideList, dda.countPeptidesList %>% select(-Method),
                     by = "peptide_sequence")
setnames(peptideList, c("Observations", "ObservationPercent"), c("Observations.DDA", "ObservationPercent.DDA"))

obs <- rbind(dia.countPeptidesList, dda.countPeptidesList)
head(obs)
ggplot(obs, aes(x = ObservationPercent)) +
  geom_histogram(aes(fill = Method), position = position_dodge())
obs %>%
  pivot_wider(id_cols = peptide_sequence, names_from = Method, values_from = ObservationPercent) %>%
  ggplot(aes(x = DDA, y = DIA)) +
  geom_point() +
  plot_theme()

peptideList %>%
  ggplot(aes(x = ObservationPercent.DDA, ObservationPercent.DIA)) +
  geom_point(aes(color = p.adj.value.DDA))

peptideList %>% filter(!is.na(log2FC.DDA) & !is.na(log2FC.DIA)) # 252 total with an intensity

peptideList %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) + 
  geom_point(aes(color = p.adj.value.DDA)) +
  scale_x_continuous(limits = c(-2, 6)) +
  scale_y_continuous(limits = c(-2, 6)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_cor(cor.coef.name = "rho") +
  xlab(expression(log[2]*"FC "[DDA])) +
  ylab(expression(log[2]*"FC "[DIA])) +
  scale_color_continuous(name = expression(FDR[DDA])) +
  plot_theme()


peptideList %>%
  filter(log2FC.DDA * log2FC.DIA > 0) %>% View() # 213 peptides

peptideList %>%
  filter(log2FC.DDA * log2FC.DIA > 0) %>%
  filter(p.adj.value.DIA < 0.05 | p.adj.value.DDA < 0.05) %>% # 32
  View()


peptideList %>%
  pivot_longer(cols = c("p.adj.value.DDA", "p.adj.value.DIA"), names_to = "Method", values_to = "FDR") %>%
  select(peptide_sequence, Method, FDR) %>%
  mutate(Method = gsub("p\\.adj\\.value\\.", "", Method)) %>%
  ggplot(aes(x = FDR)) +
  geom_histogram(bins = 100, position = position_dodge(), aes(fill = Method)) +
  scale_fill_manual(breaks = c("DDA", "DIA"), values = c("orange", "chartreuse4")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 42)) +
  plot_theme() +
  ylab("Peptides (n)") +
  theme(legend.position = c(0.9, 0.8))


peptideList %>%
  filter(log2FC.DDA * log2FC.DIA > 0) %>%
  mutate(Significant = case_when((p.adj.value.DDA < 0.05 & p.adj.value.DIA < 0.05) ~ "Both",
                                 p.adj.value.DDA < 0.05 ~ "DDA", 
                                  p.adj.value.DIA < 0.05 ~ "DIA")) %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_point(aes(color = Significant, alpha = Significant), size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_cor(cor.coef.name = "rho") +
  xlab(expression(log[2]*"FC "[DDA])) +
  ylab(expression(log[2]*"FC "[DIA])) +
  scale_alpha_manual(breaks = c("Both", "DDA", "DIA"), values = c(0.9, 0.9, 0.9), na.value = 0.5) +
  scale_color_manual(breaks = c("Both", "DDA", "DIA"), values = c("violetred3", "orange", "chartreuse4")) +
  plot_theme()


str(peptideList)
peptideList %>%
  mutate(Direction = log2FC.DDA * log2FC.DIA) %>%
  filter(Direction < 0) # 38 peptides are opposite


dda.pep.de %>% 
  filter(peptide %in% peptideList[(log2FC.DDA * log2FC.DIA) < 0]$peptide_sequence) %>%
  ggplot(aes(x = cISUP, y = Intensity, group = cISUP)) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6)

pep.de %>%
  filter(peptide_sequence %in% peptideList[(log2FC.DDA*log2FC.DIA) < 0]$peptide_sequence) %>%
  distinct(SampleID, log2Intensity, peptide_sequence, .keep_all = T) %>%
  ggplot(aes(x = cISUP, y = log2Intensity, group = cISUP)) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6)



peptide.oppo <-  dda.pep.de %>% 
  filter(peptide %in% peptideList[(log2FC.DDA * log2FC.DIA) < 0]$peptide_sequence) %>%
  distinct(Sample, Intensity, peptide) %>%
  mutate(Method = "DDA")

temp <- pep.de %>%
  filter(peptide_sequence %in% peptideList[(log2FC.DDA*log2FC.DIA) < 0]$peptide_sequence) %>%
  distinct(SampleID, log2Intensity, peptide_sequence) %>%
  mutate(Method = "DIA")

setnames(temp, c("SampleID", "log2Intensity", "peptide_sequence"),
         c("Sample", "Intensity", "peptide"))


peptide.oppo <- rbind(peptide.oppo, temp)
remove(temp)

peptide.oppo %>%
  ggplot(aes(x = Method, y = Intensity)) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.3) +
  plot_theme()
head(peptide.oppo)
peptide.oppo %>%
  pivot_wider(id_cols = c("Sample", "peptide"), names_from = "Method", values_from = "Intensity") %>%
  ggplot(aes(x = DDA, y = DIA)) +
  geom_point() +
  stat_cor()

############### Look at the trends of these peptides ############################


peptideList <- peptideList %>% filter(!(peptide_sequence %in% peptide.oppo$peptide))


peptideList.filter <- peptideList %>%
  filter(log2FC.DDA * log2FC.DIA > 0) %>%
  filter(p.adj.value.DIA < 0.05 | p.adj.value.DDA < 0.05) # 32 peptides


peptideList.filter %>%
  filter(p.adj.value.DDA < 0.05) %>%# 13
  filter(log2FC.DDA > 0)

pep.de %>% inner_join( (peptideList.filter %>%
  filter(p.adj.value.DIA < 0.05)  %>%  # 20
  filter(log2FC.DIA > 0)) , by = "peptide_sequence") %>%
  distinct(peptide_sequence, SampleID, .keep_all = T) %>%
  mutate(Label = paste(peptide_sequence, GeneName, sep = " * ")) %>%
  ggplot(aes(x = cISUP, y = log2Intensity, group = cISUP)) +
  geom_point() +
  #geom_jitter(width = 0.1) +
  stat_smooth(method = "loess", aes(group = 1)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~Label, scale = 'free_x') +
  theme_classic()

peptides.dia.DOWN <- peptideList.filter %>%
  filter(p.adj.value.DIA < 0.05)  %>%  # 9
  filter(log2FC.DIA < 0)


peptides.dia.UP <- peptideList.filter %>%
  filter(p.adj.value.DIA < 0.05)  %>%  # 20
  filter(log2FC.DIA > 0)


peptides.dda.UP <- peptideList.filter %>%
  filter(p.adj.value.DDA < 0.05)

pdf("D:/projects/pca_urine_spectral_lib/results/20220327_proBatch/20220608_dda_FDR005_Log2FC_TrendsinDIA.pdf",
    onefile = T) # change filename!!

for(i in seq_along(peptides.dda.UP$peptide_sequence)){ # change dataset
  p <- peptides.dda.UP$peptide_sequence[i]  # change dataset
  g <- peptides.dda.UP$GeneName[i]  # change dataset
  g<- pep.de %>%
    filter(peptide_sequence == p) %>%
    distinct(peptide_sequence, SampleID, .keep_all = T) %>%
    ggplot(aes(x = cISUP, y = log2Intensity, group = cISUP)) +
    geom_jitter(width = 0.1, alpha= 0.8) +
    stat_smooth(method = "loess", aes(group = 1)) +
    geom_boxplot(width = 0.8, alpha = 0.8, outlier.shape = NA) +
    plot_theme() +
    theme(title = element_text(size = 20)) +
    ylab(expression(log[2]*"Intensity")) +
    ggtitle(paste(g, ":", p, sep = " "))
  print(g)
}

dev.off()

pdf("D:/projects/pca_urine_spectral_lib/results/20220327_proBatch/20220608_dda_FDR005_Log2FC_Trends.pdf",
    onefile = T) # change filename!!

for(i in seq_along(peptides.dda.UP$peptide_sequence)){ # change dataset
  p <- peptides.dda.UP$peptide_sequence[i]  # change dataset
  g <- peptides.dda.UP$GeneName[i]  # change dataset
  g<- dda.pep.de %>%
    filter(peptide == p) %>%
    distinct(peptide, Sample, .keep_all = T) %>%
    ggplot(aes(x = cISUP, y = Intensity, group = cISUP)) +
    geom_jitter(width = 0.1, alpha= 0.8) +
    stat_smooth(method = "loess", aes(group = 1)) +
    geom_boxplot(width = 0.8, alpha = 0.8, outlier.shape = NA) +
    plot_theme() +
    theme(title = element_text(size = 20)) +
    ylab(expression(log[2]*"Intensity")) +
    ggtitle(paste(g, ":", p, sep = " "))
  print(g)
}

dev.off()

########################### Differential expression for pISUP #########################

head(pep.de)

FC.pISUP <- pep.de %>%
  inner_join(peptideMap[, c("peptide_sequence", "GeneName", "protein_id")], 
             by = c("peptide_sequence", "protein_id")) %>%
  group_by(pGrade, peptide_sequence, protein_id, GeneName) %>%
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "GeneName", "protein_id"),
              names_from = "pGrade", values_from = 'mean') %>%
  mutate(log2FC = `2+` - `1`) %>% as.data.table()

dda.pISUP <- dda.pep.de %>%
  mutate(pGrade = ifelse(pISUP == 1, "1", "2+")) %>%
  group_by(pGrade, peptide, Gene.names, Protein.names) %>%
  summarise(mean = mean(Intensity)) %>%
  pivot_wider(id_cols = c("peptide", "Gene.names", "Protein.names"),
              names_from = "pGrade", values_from = "mean") %>%
  mutate(log2FC = `2+` - `1`) %>% as.data.table()


wt <- lapply(seq_along(FC.pISUP$peptide_sequence), function(x){
  p <- FC.pISUP$peptide_sequence[x]
  t <- Htest(p, peptide_sequence, "1", "2+", pGrade, pep.de, intensity_col = log2Intensity) # change it!!
  t
})

FC.pISUP$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.pISUP$p.adj.value <- p.adjust(FC.pISUP$p.value, method = "fdr")

wt <- lapply(seq_along(dda.pISUP$peptide), function(x){
  p <- dda.pISUP$peptide[x]
  t <- Htest(p, peptide, "1", "2+", pGrade, dda.pep.de, intensity_col = Intensity) # change it!!
  t
})

dda.pISUP$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
dda.pISUP$p.adj.value <- p.adjust(dda.pISUP$p.value, method = 'fdr')
head(dda.pISUP)
FC.pISUP <- merge(FC.pISUP, dda.pISUP, by.x = c("peptide_sequence", "GeneName"),
                  by.y = c("peptide", "Gene.names"), suffixes = c(".DIA", ".DDA"), all = T)

head(FC.pISUP)
colnames(FC.pISUP)
FC.pISUP$p.adj.value.DIA <- p.adjust(FC.pISUP$p.value.DIA, method = 'fdr')
FC.pISUP$p.adj.value.DDA <- p.adjust(FC.pISUP$p.value.DDA, method = "fdr")

FC.pISUP %>% 
  ggplot(aes(x = p.adj.value.DIA)) + geom_histogram(bins = 100) + plot_theme()

FC.pISUP %>%
  filter(p.adj.value.DDA < 0.05)
# 9 DIA peptides ssignificant

plot_volcano(FC.pISUP, pvalue_column = `p.value.DIA`, padj_column = p.adj.value.DIA, 
             FC_column = log2FC.DIA, legends = T, label = F) +
  xlab(expression(log[2]*"FC (pISUP 2+ / 1)"[DIA])) +
  ylab(expression("-log"[2]*"(p-value)"))


plot_volcano(FC.pISUP, pvalue_column = `p.value.DDA`, padj_column = p.adj.value.DDA, 
             FC_column = log2FC.DDA, legends = T, label = F) +
  xlab(expression(log[2]*"FC (pISUP 2+ / 1)"[DDA])) +
  ylab(expression("-log"[2]*"(p-value)")) +
  scale_x_continuous(limits = c(-2, 6)) +
  scale_y_continuous(limits = c(0, 6))


FC.pISUP %>%
  mutate(Significant = case_when(p.adj.value.DIA < 0.05 & p.adj.value.DDA < 0.05 ~ "Both",
                                 p.adj.value.DIA < 0.05 & !(p.adj.value.DDA < 0.05) ~ "DIA" ,
                                 !(p.adj.value.DIA < 0.05) & p.adj.value.DDA < 0.05 ~ "DDA")) %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_point(aes(color = Significant)) +
  stat_cor(cor.coef.name = "rho") +
  scale_color_manual(breaks = c("DIA"), values = "red", na.value = "grey") +
  scale_x_continuous(limits = c(-2.2, 6.2)) +
  scale_y_continuous(limits = c(-2.2, 6.2)) +
#  geom_abline(slope = 1, lty = 2, lwd = 1, alpha = 0.5) +
  plot_theme() +
  ylab(expression(log[2]*"FC"[DIA])) +
  xlab(expression(log[2]*"FC"[DDA]))

FC.pISUP %>%
  ggplot(aes(x = log2FC.DIA)) + geom_histogram(bins = 100) + plot_theme()


head(FC)

dia.dotmap <- merge(FC.pISUP[p.adj.value.DIA < 0.05], FC[p.adj.value < 0.05], 
                    by = c("peptide_sequence", "protein_id", "GeneName"), all = T)

head(dia.dotmap)
all.FC <- merge(all.FC, FC.pISUP, by = c("peptide_sequence", "GeneName"),
                all = T, suffixes = c(".cISUP", ".pISUP"))

head(all.FC)
all.FC %>%
  filter(!is.na(log2FC.DIA.cISUP)) %>%
  mutate(Sig = case_when(p.adj.value.DIA.cISUP < 0.05 & p.adj.value.DIA.pISUP < 0.05 ~ "ISUP",
                         !(p.adj.value.DIA.cISUP < 0.05) & p.adj.value.DIA.pISUP < 0.05 ~ "pISUP",
                         p.adj.value.DIA.cISUP < 0.05 & !(p.adj.value.DIA.pISUP < 0.05) ~ "cISUP")) %>%
  ggplot(aes(x = log2FC.DIA.cISUP, y = log2FC.DIA.pISUP)) +
  geom_point(aes(color = Sig)) +
  stat_cor(cor.coef.name = "rho") +
  scale_color_manual(breaks = c("ISUP", "cISUP", 'pISUP'), values = c("red", "blue", "purple"),
                     name = "FDR (<0.05)") +
  plot_theme() +
  xlab(expression(log[2]*"FC"[cISUP])) +
  ylab(expression(log[2]*"FC"[pISUP]))


pISUP.lst <- FC.pISUP %>%
  #filter(!is.na(log2FC.DIA) & !is.na(log2FC.DDA)) %>%
  filter(p.adj.value.DIA < 0.2 | p.adj.value.DDA < 0.2)


create.dotmap(pISUP.lst[p.adj.value.DDA < 0.2, c("log2FC.DDA", "log2FC.DIA")],
yaxis.cex = 1.5,
xaxis.cex = 1.5,
na.pch = 1,
na.spot.size = 0,
yaxis.lab = paste(pISUP.lst$peptide_sequence,
                  pISUP.lst$GeneName, sep = " * "),
xaxis.lab = c("DDA", "DIA"), 
spot.size.function = spot.size.function,
spot.colour.function = spot.colour.function,
key = list(
  space = "right",
  points = list(
    cex = spot.size.function(seq(-1, 1, 0.2)),
    col = spot.colour.function(seq(-1, 1,0.2)),
    pch = 19
  ),
  text = list(
    lab = as.character(seq(-1, 1, 0.2)),
    cex = 1.5,
    adj = 1.0,
    fontface = "bold"
  )
),
# control spacing at top of key
key.top = 1,
pch = 21,
pch.border.col = "white",
# add the background
bg.data = pISUP.lst[p.adj.value.DDA < 0.2, c("p.adj.value.DDA", "p.adj.value.DIA")],
# add a colourkey
colourkey = TRUE,
# set colour scheme for background data
colour.scheme = c("black", "white"),
# make bg colour scheme a discrete colour scheme, with breaks at these places
at = rev(seq(0, 0.5, 0.05))
)


pISUP.lst


pep.de %>% filter(peptide_sequence %in% pISUP.lst$peptide_sequence[1]) %>%
  ggplot(aes(x = pISUP, y = log2Intensity, group = pISUP)) +
  geom_jitter(width = 0.1, alpha = 0.9) +
  stat_smooth(method = 'loess', aes(group = 1)) +
  geom_boxplot(width = 0.8, alpha = 0.8, outlier.shape = NA) + plot_theme() +
  annotate("text", label = "FDR", x = 4 , y = 30)

pdf("D:/projects/pca_urine_spectral_lib/results/20220327_proBatch/20220622_dia_FDR02_pISUP_Log2FC_TrendsInDDA.pdf",
    onefile = T) # change filename!!

for(i in seq_along(pISUP.lst[p.adj.value.DIA < 0.05]$peptide_sequence)){ # change dataset
  p <- pISUP.lst[p.adj.value.DIA < 0.05]$peptide_sequence[i]  # change dataset
  g <- pISUP.lst[p.adj.value.DIA < 0.05]$GeneName[i]  # change dataset
  g<- dda.pep.de %>%
    filter(peptide == p) %>%
    distinct(peptide, Sample, .keep_all = T) %>%
    ggplot(aes(x = pISUP, y = Intensity, group = pISUP)) +
    geom_jitter(width = 0.1, alpha= 0.8) +
    stat_smooth(method = "loess", aes(group = 1)) +
    geom_boxplot(width = 0.8, alpha = 0.8, outlier.shape = NA) +
    plot_theme() +
    theme(title = element_text(size = 20)) +
    ylab(expression(log[2]*"Intensity")) +
    ggtitle(paste(g, ":", p, sep = " "))
  print(g)
}

dev.off()

patients %>% filter(cISUP %in% c(4,5)) %>%
  select(`Sample ID`, cISUP, CGG1, CGG2, pISUP, PGG1, PGG2) %>% 
  filter(pISUP != cISUP) %>% View()

patients %>% filter(pISUP %in% c(4,5))

pep.de %>% distinct(SampleID, peptide_sequence, cISUP, pISUP) %>%
  group_by(cISUP,pISUP,  SampleID) %>%
  summarise(PeptideCount = n()) %>%
  ggplot(aes(x = pISUP, y = PeptideCount, group = pISUP)) + 
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  theme_classic()

############################ Differential expression for UPGRADED ####################

pep.de %>% filter(cISUP == 1) %>%
  distinct(SampleID, cISUP, pISUP) %>% # 81 cISUP 1 patients
  mutate(pISUP = factor(pISUP, levels = seq(1, 5))) %>%
  ggplot(aes(x = pISUP)) +
  geom_bar(stat = 'count', fill = "lightgrey", color = "black", size = 1.2, width = 0.8) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, size = 6) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 60)) +
  plot_theme() + xlab('pISUP') + ylab("Samples (n)") +
  ggtitle("cISUP1 patients (81)")
  
# 28 that upgraded

FC.UP <- pep.de %>%
  filter(cISUP == 1) %>%
  mutate(Upgraded = ifelse(pISUP > 1, TRUE, FALSE)) %>%
  inner_join(peptideMap[, c("peptide_sequence", "GeneName", "protein_id")], 
             by = c("peptide_sequence", "protein_id")) %>%
  group_by(Upgraded, peptide_sequence, protein_id, GeneName) %>%
  summarise(mean = mean(log2Intensity)) %>%
  pivot_wider(id_cols = c("peptide_sequence", "GeneName", "protein_id"),
              names_from = "Upgraded", values_from = 'mean') %>%
  mutate(log2FC = `TRUE` - `FALSE`) %>% as.data.table()

FC.UP %>% ggplot(aes(x = p.adj.value)) + geom_histogram(bins = 500) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 110)) +
  plot_theme() +
  xlab(expression(log[2]*"FC(Upgraded vs Not upgraded)")) +
  ylab("Peptides (n)")

pep.UP <- pep.de %>% filter(cISUP == 1) %>% mutate(Upgraded = ifelse(pISUP > 1, TRUE, FALSE)) %>%
  as.data.table()

wt <- lapply(seq_along(FC.UP$peptide_sequence), function(x){
  p <- FC.UP$peptide_sequence[x]
  t <- Htest(p, peptide_sequence, "TRUE", "FALSE", Upgraded, pep.UP, intensity_col = log2Intensity) # change it!!
  t
})

FC.UP$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
FC.UP$p.adj.value <- p.adjust(FC.UP$p.value, method = "fdr")

plot_volcano(FC.UP, pvalue_column = `p.value`, padj_column = `p.adj.value`, FC_column = `log2FC`) +
  xlab(expression(log[2]*"FC(Upgraded vs ISUP 1)"))

dda.UP <- dda.pep.de %>%
  filter(cISUP == 1) %>%
  mutate(Upgraded = ifelse(pISUP > 1, TRUE, FALSE)) %>%
  group_by(Upgraded,peptide, Gene.names, Protein.names) %>%
  summarise(mean = mean(Intensity)) %>%
  pivot_wider(id_cols = c("peptide", "Gene.names", "Protein.names"),
              names_from = "Upgraded", values_from = 'mean') %>%
  mutate(log2FC = `TRUE` - `FALSE`) %>% as.data.table()
dda.pep.UP <- dda.pep.de %>%
  filter(cISUP == 1) %>%
  mutate(Upgraded = ifelse(pISUP > 1, TRUE, FALSE))

wt <- lapply(seq_along(dda.UP$peptide), function(x){
  p <- dda.UP$peptide[x]
  t <- Htest(p, peptide, "TRUE", "FALSE", Upgraded, dda.pep.UP, intensity_col = Intensity) # change it!!
  t
})

dda.UP$p.value <- lapply(wt, function(x) x$p.value) %>% unlist()
dda.UP$p.adj.value <- p.adjust(dda.UP$p.value, method = "fdr")

dda.UP %>% ggplot(aes(x = p.adj.value)) + geom_histogram(bins = 500) + plot_theme() 

plot_volcano(dda.UP, pvalue_column = `p.value`, padj_column = `p.adj.value`, FC_column = `log2FC`) +
  xlab(expression(log[2]*"FC(Upgraded vs ISUP 1)"))

all.UP <- merge(FC.UP, dda.UP, by.x = c("peptide_sequence", "GeneName"),
                by.y = c("peptide", "Gene.names"), suffixes = c(".DIA", ".DDA"), all = T)

all.UP %>% ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_hex(bins = 100) +
  stat_cor(cor.coef.name = 'rho') +
  scale_fill_viridis_c(option = "magma", begin = 0.2) +
  plot_theme() +
  xlab(expression(log[2]*"FC"[DDA])) +
  ylab(expression(log[2]*"FC"[DIA])) +
  scale_x_continuous(limits = c(-4, 6))

all.UP %>% 
  mutate(Significant = case_when(p.adj.value.DDA < 0.25 & p.adj.value.DIA < 0.25 ~ "Both",
                                 p.adj.value.DDA < 0.25 & !(p.adj.value.DIA < 0.25) ~ "DDA",
                                 !(p.adj.value.DDA < 0.25) & p.adj.value.DIA < 0.25~ "DIA")) %>%
  ggplot(aes(x = log2FC.DDA, y = log2FC.DIA)) +
  geom_point(aes(color = Significant, alpha = Significant)) +
  stat_cor(cor.coef.name = 'rho') +
  scale_color_manual(breaks = c("Both", "DDA", "DIA"), values = c("red", "blue", "purple"), 
                     name = "FDR < 0.2" ) +
  scale_alpha_manual(breaks = c("Both", "DDA", "DIA"), values = c(1,1,1), na.value = 0.6, 
                     name = "FDR < 0.2") +
  geom_vhlines(xintercept = 0, yintercept = 0) +
  plot_theme() +
  xlab(expression(log[2]*"FC"[DDA])) +
  ylab(expression(log[2]*"FC"[DIA])) +
  scale_x_continuous(limits = c(-4, 6))



