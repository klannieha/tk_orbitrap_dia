#### This script is to get Gene lists for Enrichment Analysis ##########
library(data.table)

library(ggpubr)
library(ggrepel)
library(ggpmisc)

# load data from RData

################## Formate gene lists ###########3#######################
# Lists:
# 1. cISUP - DIA
# 2. cISUP - DDA
# 3. pISUP - DIA
# 4. pISUP - DDA
# 5. UP - DIA
# 6. UP - DDA

# cISUP - DIA

head(all.FC)
# using gProfiler:
dia.gPro <- all.FC %>% filter(!is.na(p.adj.value.DIA.cISUP)) %>%
  select(peptide_sequence, GeneName, protein_id.cISUP, `1.DIA.cISUP`, `2+.DIA.cISUP`, log2FC.DIA.cISUP,
         p.value.DIA.cISUP, p.adj.value.DIA.cISUP) %>%
    as.data.table()

head(dia.gPro)
# filter duplicated GeneNames
dia.gPro <- dia.gPro[order(p.adj.value.DIA.cISUP)]
dia.GeneList <- dia.gPro[!duplicated(GeneName)]$GeneName
dia.GeneList <- lapply(dia.GeneList, function(x) GeneFormat(x)) %>% unlist()
head(dia.GeneList)
write.table(dia.GeneList, 
            file = "D:/projects/pca_urine_spectral_lib/data/openswath/20220628_PathwayEnrichment/20220628_dia_cISUP_gProfilerGeneList.txt",
            sep = "\n", row.names = F, col.names = F, quote = F)
