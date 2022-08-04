############## Script for comparing libraries ####################
library(data.table)
library(ggplot2)
library(tidyverse)
library(BoutrosLab.plotting.general)
library(ggpubr)
library(eulerr)
library(aLFQ)
library(SWATH2stats)

source("D:/projects/pca_urine_spectral_lib/src/Amanda_source.R")

###################### load library data ######################

folder <- list.dirs("D:/projects/pca_urine_spectral_lib/data/library/202207_library_benchmark_uEPS", recursive = F)

folder

fragpipe.lib <- list.files(folder[4], pattern = "library.tsv", full.names = T)
fragpipe.lib
fragpipe.lib <- fread(fragpipe.lib[1])

mq.files <- list.files(paste0(folder[3], "/combined/txt/"), pattern = ".txt", full.names = T)
mq.files
mq.peptides <- read.table(mq.files[11], header = T, sep = "\t")
mq.allPeptides <- read.table(mq.files[1], header = T, sep = "\t")
mq.evidence <- read.table(mq.files[2], header = T, sep= "\t")

head(mq.allPeptides)
mq.lib <- list.files(paste0(folder[3], "/combined_dia/txt/"), full.names = T)
mq.lib
mq.lib <- read.table(mq.lib[2], header = T, sep = "\t")

tpp.lib <- list.files(folder[6], pattern = ".tsv", full.names = T)
tpp.lib
tpp.lib <- fread(tpp.lib[2])

######################### peptide ###############################

head(fragpipe.lib)

fragpipe.lib %>% filter(grepl("DECOY", ProteinId)) # no weird accessions check

fragpipe.lib %>% distinct(PeptideSequence) %>% summarise(Count = n())
fragpipe.pr <- fragpipe.lib %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>%
  summarise(Count = n())

head(tpp.lib)

tpp.lib %>% distinct(PeptideSequence) %>% summarise(Count = n())
tpp.lib %>% distinct(FullUniModPeptideName, PrecursorCharge) %>% summarise(Count = n())
tpp.pr <- tpp.lib %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>% 
  summarise(Count = n())
tpp.pr
tpp.tr <- tpp.lib %>% distinct(TransitionGroupId) %>% summarise(Count = n())

head(tpp.lib)
head(mq.peptides)
tail(mq.allPeptides)
head(mq.lib)
mq.lib %>% distinct(Sequence) %>% summarise(Count = n())
mq.lib %>% distinct(Library.entry.index) %>% summarise(Count = n())
mq.peptides %>% distinct(Sequence) %>% summarise(Count = n())
mq.allPeptides %>% filter(!is.na(Sequence)) %>% distinct(Charge, Modified.sequence)
head(mq.evidence)

mq.evidence %>% distinct(Modified.sequence, Charge) %>% summarise(Count = n())

mq.pr <- mq.lib %>% distinct(Modified.sequence, Precursor.charge) %>% summarise(Count = n())
mq.pr


# plot peptide count

libs.pepCount <- data.frame(Library = c("Fragpipe", "MQ", "TPP"),
                            Peptides = c(length(unique(fragpipe.lib$PeptideSequence)),
                                         length(unique(mq.lib$Sequence)),
                                         length(unique(tpp.lib$PeptideSequence))),
                            Precursors = c(fragpipe.pr$Count, mq.pr$Count, tpp.pr$Count),
                            stringsAsFactors = F)

create.barplot(Peptides~Library, libs.pepCount, col = "darkgrey", 
               border.lwd = 3, ylimits = c(0, 30000), text.above.bars = libs.pepCount$Peptides)

ggplot(libs.pepCount, aes(x = Library, y = Peptides)) +
  geom_col(width = 0.6, fill = 'darkgrey', color = 'black', size = 1.2) +
  geom_text(aes(label = Peptides), vjust = -0.4, size = 5 ) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 31000)) +
  plot_theme()

ggplot(libs.pepCount, aes(x = Library, y = Precursors)) +
  geom_col(width = 0.6, fill = 'darkgrey', color = 'black', size = 1.2) +
  geom_text(aes(label = Precursors), vjust = -0.4, size = 5 ) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 45000)) +
  plot_theme()

######################### compare between the 3 libraries ########

head(fragpipe.lib)
head(mq.peptides)
head(tpp.lib)
head(mq.lib)
set.seed(123)
peptides.lst <- list(
                     TPP = unique(tpp.lib$PeptideSequence),
                     Fragpipe = unique(fragpipe.lib$PeptideSequence),
                     MQ = unique(mq.lib$Sequence)
                     )

fit <- euler(peptides.lst, shape = "ellipse")

plot(fit,
     fills = c("#FAE5A1", "lightcoral","purple4"),
     col = 'white',
     #col = c("gold", "dodgerblue", "#99d8c9"),
     edges = T,
     quantities = list(fontsize = 18), adjust_labels = T, lwd = 2, 
     legend = list(fontsize = 24, side = "bottom", nrow = 1, ncol = 4), rotation = 1)

tpp.lib$Tr_recalibrated <- as.numeric(tpp.lib$Tr_recalibrated)
fragpipe.lib$NormalizedRetentionTime <- as.numeric(fragpipe.lib$NormalizedRetentionTime)

tpp.lib %>% distinct(FullUniModPeptideName, PrecursorCharge, FragmentSeriesNumber, 
                     LibraryIntensity) %>%
  group_by(FullUniModPeptideName, PrecursorCharge) %>%
  summarise(MedianIntensity = median(LibraryIntensity))

tpp.lib$NormalizedRetentionTime <- tpp.lib$Tr_recalibrated

fragpipe.lib %>% arrange(desc(LibraryIntensity)) %>%
  distinct(PeptideSequence, PrecursorCharge, NormalizedRetentionTime)

mq.lib %>% distinct(Sequence, Precursor.charge, Retention.time)

fragpipe.lib %>% distinct(PeptideSequence, PrecursorCharge, 
                          NormalizedRetentionTime) %>% 
  inner_join(mq.lib %>% 
               distinct(Sequence, Precursor.charge, Retention.time), 
            by = c("PeptideSequence" = "Sequence", 
                   "PrecursorCharge" = "Precursor.charge"), 
            suffix = c(".fragpipe", ".MQ")) %>%
  ggplot(aes(x = NormalizedRetentionTime, y = Retention.time)) +
  geom_point() +
  stat_cor(cor.coef.name = "rho")

fragpipe.lib %>%
  ggplot(aes(x = NormalizedRetentionTime)) +
  geom_density()

tpp.lib %>%
  ggplot(aes(x = NormalizedRetentionTime)) +
  geom_density()

fragpipe.lib %>%
  distinct(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime, PeptideSequence) %>%
  inner_join(mq.lib %>%
               distinct(Sequence, Precursor.charge, Retention.time),
             by = c("PrecursorCharge" = "Precursor.charge", "PeptideSequence" = "Sequence")) %>%
  # inner_join(tpp.lib %>% distinct(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime),
  #            by = c("ModifiedPeptideSequence", "PrecursorCharge"),
  #            suffix = c(".fragpipe", ".tpp")) %>%
  ggplot(aes(x = NormalizedRetentionTime, y = Retention.time)) +
  geom_point(alpha = 0.8) + 
  stat_cor(cor.coef.name = "rho") +
  plot_theme() +
  xlab("iRT - Fragpipe") + ylab("iRT - MQ")


head(fragpipe.lib)


########################## Compare DIA analysis output ###############
  
folder

dia.folder <- list.dirs("D:/projects/pca_urine_spectral_lib/data/openswath/20220714_libray_benchmark_osw/",
                        full.names = T)
fragpipe.14mz <- list.files(dia.folder[2], pattern = "diann-output.tsv", recursive = T, full.names = T)

dia.folder
fragpipe.dia <- list.files(folder[5], pattern = "output.tsv",  full.names = T, recursive = T)
mq.dia <- list.files(paste0(folder[3], "/combined_dia/txt/"), full.names = T)
mq.dia

fragpipe.dia <- lapply(fragpipe.dia, fread)
fragpipe.dia <- do.call('rbind', fragpipe.dia) %>% as.data.table()

fragpipe.dia$Run <- str_sub(fragpipe.dia$Run, -2 , 1)

head(fragpipe.dia)

fragpipe.dia %>%
  group_by(Run) %>%
  distinct(Stripped.Sequence) %>%
  summarise(Peptides = n()) %>%
  mutate(Run = str_sub(Run,-2, -1)) %>% 
  ggplot(aes(x = Run, y = Peptides)) +
  geom_col(width = 0.6, fill = 'darkgrey', colour = 'black', size = 1.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12000)) +
  plot_theme()

head(fragpipe.dia)
fragpipe.dia$Q.Value %>% max()
fragpipe.dia %>%
  filter(Stripped.Sequence %in% iRTs$Sequence) %>%
  filter(Precursor.Charge == 2) %>%
  ggplot(aes(x = Stripped.Sequence, y = log2(Precursor.Normalised))) +
  geom_jitter(width = 0.1) + 
  geom_boxplot(width = 0.8, alpha = 0.7) + plot_theme() +
  ylab(expression(log[2]*"Intensity")) + xlab("Sequence") +
  theme(axis.text.x = element_text(angle = 90))

fragpipe.dia %>% 
  mutate(Run = str_sub(Run,-2, -1)) %>% 
  mutate(PeakWidth = (RT.Stop - RT.Start)*60) %>% 
  ggplot(aes(x = Run, y = PeakWidth)) +
  geom_jitter(width = 0.1 , alpha = 0.7) +
  geom_boxplot(width = 0.8, alpha = 0.9, outlier.shape = NA) +
  plot_theme()

fragpipe.dia %>% mutate(PeakWidth = (RT.Stop - RT.Start)*60,
                        Run = str_sub(Run, -2, -1)) %>%
  filter(Run == "01") %>%
  ggplot(aes(x = RT, y = PeakWidth)) +
  geom_point()

fragpipe.dia %>%
  mutate(PeakWidth = (RT.Stop - RT.Start)*60,
         Run = str_sub(Run, -2, -1)) %>%
  filter(Run == "01") %>% 
  arrange(RT) %>% 
  select(Run, Genes, Modified.Sequence, Stripped.Sequence, RT, RT.Start,
         RT.Stop, PeakWidth, iRT, Predicted.RT) %>%
  View()

fragpipe.dia %>% 
  mutate(Run = str_sub(Run,-2, -1)) %>% 
  mutate(PeakWidth = (RT.Stop - RT.Start)*60) %>% 
  ggplot(aes(x = RT, y = PeakWidth)) +
  geom_point()

fragpipe.dia %>%
  mutate(Run = str_sub(Run, -2, -1)) %>%
  distinct(Stripped.Sequence, Run) %>%
  group_by(Stripped.Sequence) %>%
  summarise(Detection = n()) %>%
  ggplot(aes(x = Detection)) + 
  geom_bar(stat = "count", fill = 'lightgrey', color = 'black', size = 1.2) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4.5)  +
  plot_theme()


fragpipe.dia %>% head()

mq.dia
mq.dia.allPeptides <- read.table(mq.dia[1], sep = "\t", header = T)
mq.dia.ev <- read.table(mq.dia[3], sep = "\t", header = T)
head(mq.dia.allPeptides)

mq.dia.ev %>%
  distinct(Sequence, Raw.file) %>%
  mutate(Run = str_sub(Raw.file, -2, -1)) %>%
  ggplot(aes(x = Run)) +
  geom_bar(stat = 'count', fill = 'lightgrey', color = 'black', size = 1.2, width = 0.8) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3) +
  plot_theme()

mq.dia.ev$Run <- str_sub(mq.dia.ev$Raw.file, -2 , -1)
head(mq.dia.ev$Retention.time)
mq.dia.allPeptides$Run <- str_sub(mq.dia.allPeptides$Raw.file, -2 , -1)

ggplot(mq.dia.allPeptides, aes(x = Run, y = Retention.length)) + 
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(width = 0.3, alpha = 0.7) +
  plot_theme() + ylab("Peak width (s)")


head(mq.dia.allPeptides)
mq.dia.allPeptides %>%
  distinct(Raw.file, Sequence)

osw <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/20220714_libray_benchmark_osw/openswath/",
                  pattern = "tpplib(.*).tsv", full.names = T)
osw
osw <- fread(osw)
osw <- osw[peak_group_rank == 1 & decoy == 0]
osw %>%
  filter(peak_group_rank == 1 & decoy == 0) %>%
  distinct(Sequence, filename) %>%
  group_by(filename) %>%
  summarise(Peptides = n()) %>%
  mutate(Run = str_sub(filename, -10, -9)) %>%
  ggplot(aes(x = Run, y = Peptides)) +
  geom_col(width = 0.6, fill = 'darkgrey', colour = 'black', size = 1.2) +
  geom_text(aes(label = Peptides), vjust = -0.3) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12000)) +
  plot_theme()

head(osw)

osw %>%
  filter(peak_group_rank == 1 & decoy == 0) %>%
  distinct(Sequence, filename, .keep_all = T) %>%
  mutate(Run = str_sub(filename, -10, -9)) %>%
  mutate(PeakWidth = rightWidth - leftWidth) %>%
  ggplot(aes(x = Run, y = PeakWidth)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(alpha = 0.6) +
  plot_theme()



# load tpp library, analyzed with diann
dia.folder
tpp_diann <- list.files(dia.folder[3], recursive = T, pattern = "diann-output.tsv", full.names = T)
tpp_diann

tpp_diann <- lapply(tpp_diann, fread)
tpp_diann <- do.call('rbind', tpp_diann) %>% as.data.table()
str(tpp_diann)
tpp_diann$Run <- str_sub(basename(tpp_diann$File.Name), -7, -6)
tpp_diann %>%
  mutate(PeakWidth = (RT.Stop - RT.Start)*60) %>%
  ggplot(aes(x = Run , y = PeakWidth)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot()

dia.folder
fragpipe.dia <- list.files(dia.folder[6], recursive = T, # CHANGE IT ###
                              pattern = 'diann-output.tsv', full.names = T)
fragpipe.dia
fragpipe.dia <- lapply(fragpipe.dia, fread)
fragpipe.dia <- do.call('rbind', fragpipe.dia) %>% as.data.table()
fragpipe.dia$Run <- str_sub(basename(fragpipe.dia$File.Name), -7, -6)
fragpipe.dia$PeakWidth <- (fragpipe.dia$RT.Stop - fragpipe.dia$RT.Start)*60
ggplot(fragpipe.dia, aes(x = Run, y = PeakWidth)) +
  #geom_point(alpha= 0.8) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(alpha = 0.8, width = 0.3) +
  plot_theme() + ylab("Peak width (s)")

fragpipe.dia %>% arrange(PeakWidth) %>% tail()
  filter(Genes == "TF") %>% 
  select(RT, PeakWidth, Stripped.Sequence) %>% View()
  #filter(Stripped.Sequence %in% iRTs$Sequence)

fragpipe.xic <- list.files(dia.folder[6], recursive = T, 
                                pattern = "XIC", full.names = T)
fragpipe.xic

fragpipe.xic <- lapply(fragpipe.xic, fread)
head(fragpipe.xic[[1]])
fragpipe.xic <- do.call('rbind', fragpipe.xic) %>% as.data.table()
fragpipe.xic %>% filter(MS.Level == 2) %>%
  mutate(Run = str_sub(basename(File.Name), -7, -6)) %>%
  filter(Run == "01") %>%
  group_by(Precursor.Id, FragmentSeriesNumber)

test <- fragpipe.xic %>% 
  filter(MS.Level == 2 & File.Name == "D:\\Data\\Annie\\20211219_MStern_DIA\\20201219_16mz_staggered_45min_02.mzML") %>%
  filter(Precursor.Id == "ADVTPADFSEWSK2" & Intensities != 0)
  
test$SumInt <- rowSums(test[, 15:115])
test <- test[SumInt > 1]

test %>% 
  pivot_longer(!c("File.Name", "Precursor.Id", "Modified.Sequence",
                                "Stripped.Sequence", "Q.Value", "MS.Level", "Intensities", "Retention.Times",
                                "Theoretical.Mz", "Reference.Intensity", "FragmentType", "FragmentCharge", 
                                "FragmentSeriesNumber", "FragmentLossType"), names_to = "DataPoint") %>%
  mutate(DataPoint = as.numeric(DataPoint),
         FragmentSeriesNumber = as.character(FragmentSeriesNumber)) %>%
  ggplot(aes(x = DataPoint, y = value, color = FragmentSeriesNumber)) +
  geom_point() + geom_line()


# Quant -------------------------------------------------------------------
library(diann)
library(SWATH2stats)
library(aLFQ)

head(fragpipe.dia)

glimpse(fragpipe.dia)

fragpipe.dia %>%
  ggplot(aes(x = Run, y = log2(Precursor.Quantity))) +
  geom_jitter(alpha = 0.6, width = 0.1) +
  geom_boxplot(width = 0.3, alpha = 0.8) +
  plot_theme() +
  ylab(expression(log[2]*"(Precursor Intensity)"))

fragpipe.PrecursorInt <- subset(fragpipe.dia, 
                                select = c("Run", "Precursor.Quantity","Precursor.Id", 
                                           "Genes", "Protein.Group", "Modified.Sequence", 
                                           "Stripped.Sequence", "Precursor.Charge"))

head(fragpipe.PrecursorInt)

fragpipe.PrecursorInt %>% 
  pivot_wider(id_cols = c( "Precursor.Id", "Protein.Group", "Genes", 
                          "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge"),
              names_from = Run, values_from = "Precursor.Quantity" ) %>%
  ggplot(aes(x = log2(`02`), y = log2(`03`))) +
  geom_point() +
  stat_cor(cor.coef.name = "rho") +
  plot_theme() +
  xlab("rep 2") + ylab("rep 3") + ggtitle(expression(log[2]*"(Precursor Intensity)"))


fragpipe.PrecursorInt %>% 
  distinct(Modified.Sequence) %>% head()
fragpipe.dia$File.Name <- fragpipe.dia$Run

# this function takes the highest Int peptide precursor Intensity
diann.peptides <- diann_matrix(fragpipe.dia, id.header = "Stripped.Sequence")

head(diann.peptides)
diann.CV <- getCV(as.data.table(diann.peptides))
diann.CV$Sequence <- rownames(diann.peptides)

tpp_diann.peptides <- diann_matrix(tpp_diann, id.header = "Stripped.Sequence")
tpp_diann.CV <- getCV(as.data.table(tpp_diann.peptides))
tpp_diann.CV$Sequence <- rownames(tpp_diann.peptides)
head(tpp_diann.CV)

ggplot(tpp_diann.CV, aes(x = CV)) + geom_density()
median(tpp_diann.CV$CV, na.rm = T)

# Quant for OSW -----------------------------------------------------------

osw %>% distinct(FullPeptideName)
head(osw)
osw.PrecursorInt <- subset(osw, select = c("filename", "Sequence", "FullPeptideName", "Charge"))

osw <- reduce_OpenSWATH_output(osw)
osw.annotate <- osw %>% distinct(filename) %>% 
  mutate(BioReplicate = str_sub(basename(filename), -10,-9), 
         Condition = "16mz staggered") %>%
  as.data.table()
osw.annotate$Run <- osw.annotate$BioReplicate
setnames(osw.annotate,"filename", "Filename")
osw.annotate <- as.data.frame(osw.annotate)
osw <- as.data.frame(osw)
osw <- sample_annotation(osw, osw.annotate)

disagg <- disaggregate(osw)
alfq <- convert4aLFQ(disagg, check_transitions = F)

osw.peptides <- PeptideInference(alfq, transition_topx = 6, transition_strictness = 'strict',
                 transition_summary = 'sum', consensus_transitions = F)
osw.peptides

osw.peptides <- dcast(osw.peptides, peptide_sequence~run_id,
      drop = FALSE, fill = NaN, max, value.var = "peptide_intensity")

head(osw.peptides)
osw.CV <- getCV(osw.peptides[, 2:4])
osw.CV$Sequence <- osw.peptides$peptide_sequence
osw %>% 
  mutate(PeakWidth = rightWidth - leftWidth) %>% 
  filter(PeakWidth > 100) %>% 
  distinct(Sequence, ProteinName)


# Compare quant -----------------------------------------------------------

CVs <- merge(diann.CV, osw.CV, by = "Sequence", suffixes = c(".Fragpipe", ".OSW"))
CVs %>%
  ggplot(aes(x = log2(mean.OSW), y = log2(mean.Fragpipe))) +
  geom_point() +
  stat_cor(cor.coef.name = "rho") + plot_theme() +
  xlab("OpenSWATH") + ylab("Fragpipe (DIANN)") +
  ggtitle(expression(log[2]*"(Peptide Intensity)"))

CVs %>% 
  ggplot(aes(x = log2(mean.Fragpipe), y = CV.Fragpipe)) +
  geom_point()

CVs %>%
  pivot_longer(cols = c("CV.Fragpipe", "CV.OSW"), names_to = "Algorithm", values_to = "CV") %>%
  mutate(Algorithm = gsub("CV\\.", "", Algorithm)) %>%
  ggplot(aes(x = CV, fill = Algorithm)) +
  geom_density(alpha = 0.6, color = "lightgrey") +
  scale_fill_brewer(palette = "Set1") +
  # ggplot(aes(x = Algorithm, y = CV)) +
  # geom_jitter(width = 0.1, alpha = 0.6) +
  # geom_boxplot(alpha = 0.8, width = 0.3) +
  plot_theme() + theme(legend.position = c(0.8, 0.8))

CVs %>% 
  summarise(MedianDiann = median(CV.Fragpipe, na.rm = T),
            MedianOSW = median(CV.OSW, na.rm = T))


