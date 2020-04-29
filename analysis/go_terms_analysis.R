########## This script is to analyze the GO term analysis ##########

library(data.table)
library(tidyr)
library(ggplot2)
library(VennDiagram)

########### Load the data ############
# data was from list of genes shared between directEPS and postDRE urine
basefolder <- "D:/projects/pca_urine_spectral_lib/results/"

gProfiler <- list.files(basefolder, pattern = "gProfiler", full.names = T)

####### Read the first one ##########
GO_MF <- fread(gProfiler[1])
max(GO_MF$adjusted_p_value)
GO_MF <- GO_MF[adjusted_p_value <= 0.01]
GO_MF <- GO_MF[order(-negative_log10_of_adjusted_p_value)]
GO_MF$term_name <- factor(GO_MF$term_name, levels = GO_MF$term_name)
GO_MF_sub <- GO_MF[1:20]
GO_MF_sub$term_name <- factor(GO_MF_sub$term_name, levels = GO_MF_sub$term_name)

ggplot(GO_MF_sub, aes(x = term_name, y = negative_log10_of_adjusted_p_value, fill = log(negative_log10_of_adjusted_p_value))) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
  scale_x_discrete(limits=rev(levels(GO_MF_sub$term_name))) +
  scale_fill_viridis_c(direction = -1) + coord_flip() + labs(fill = "log10_pvalue") +
  theme_classic(base_size = 16)

GO_MF_cell_adhesion <- GO_MF[term_name == "cell adhesion molecule binding",] %>% as.character()
genes <- grep("[A-Z],", GO_MF_cell_adhesion)
genes <- colnames(GO_MF)[genes]
