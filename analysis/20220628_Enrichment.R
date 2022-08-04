#### This script is to get Gene lists for Enrichment Analysis ##########
library(data.table)

library(ggpubr)
library(ggrepel)
library(ggpmisc)
library(topGO)
library(org.Hs.eg.db)
library(BoutrosLab.plotting.general)

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
dia.GeneList <- dia.gPro[!duplicated(GeneName)]
dia.GeneList <- subset(dia.GeneList, select = c("GeneName", "p.value.DIA.cISUP"))
dia.GeneLst <- dia.GeneList$p.value.DIA.cISUP
names(dia.GeneLst) <- dia.GeneList$GeneName
dia.GeneList <- lapply(dia.GeneList, function(x) GeneFormat(x)) %>% unlist()
head(dia.gPro)
dia.GeneExp <- dia.gPro[order(-log2FC.DIA.cISUP)]
dia.GeneExp <- dia.GeneExp[!duplicated(GeneName)]
dia.GeneExp <- subset(dia.GeneExp, select = c("GeneName", "log2FC.DIA.cISUP"))

head(dia.GeneList)
head(dia.GeneExp)
write.table(dia.GeneList, 
            file = "D:/projects/pca_urine_spectral_lib/data/openswath/20220628_PathwayEnrichment/20220628_dia_cISUP_gProfilerGeneList.txt",
            sep = "\n", row.names = F, col.names = F, quote = F)
write.table(dia.GeneExp, 
            file = "D:/projects/pca_urine_spectral_lib/data/openswath/20220628_PathwayEnrichment/20220628_dia_cISUP_gProfilerGeneExp.txt",
            sep = "\t", row.names = F, col.names = F, quote = F)

################# First format gene list ###################
# required data:
# allGenes  named vector with p.values from DiffE and names from genes
# geneSel   function of gene selection of diffE genes
# annot function for annotation
# Gene2GO   named list of genes mapping to GO terms

# getting the named list
head(dia.GeneList)
dia.GeneList <- dia.gPro$p.adj.value.DIA.cISUP
names(dia.GeneList) <- dia.gPro$GeneName
#dia.GeneList <- dia.GeneList[1:10]

# assigning geneSel function
topDiffGenes <- function(score){
  return(score < 0.05) # changed from 0.01 FDR to 0.05 FDR for more interesting genes
}

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

#dia.GeneLst <- head(dia.GeneLst)
top10 <- dia.gPro[p.adj.value.DIA.cISUP < 0.01]$GeneName %>% unique()
  
# formatting named list of GO term and gene name mapping

map <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
gene2Go <- inverseList(map)
head(map)
head(gene2Go)
sampleGO <- new("topGOdata", ontology = "BP", allGenes = dia.GeneList, 
                geneSel = topDiffGenes,
                annot = annFUN.gene2GO, gene2GO = gene2Go)

##################### Start analysis of enrichment tests ###########

resultFisher <- runTest(sampleGO, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(sampleGO, algorithm = "classic", statistic = "ks")
resultKS
resultKS.elim <- runTest(sampleGO, algorithm = "elim", statistic = "ks")
resultKS.elim

res.dt <- GenTable(sampleGO,  classicFisher = resultFisher,
         classicKS = resultKS, elimKS = resultKS.elim,
         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)

# p-values
pV.classic <- score(resultFisher)
pV.adj.classic <- p.adjust(pV.classic, method = "BH")
pV.KS <- score(resultKS)
pV.adj.KS <- p.adjust(pV.KS, method = "BH")
pV.KSelim <- score(resultKS.elim)
pV.adj.KSelim <- p.adjust(pV.KSelim, method = "BH")
# summarise the statistics
gStat <- termStat(sampleGO, names(pV.classic))
gStat
gSize <- gStat$Annotated / max(gStat$Annotated) * 4
gSize
gCol <- colMap(gStat$Significant)

plot(pV.classic, pV.KSelim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)


res.dt
sel.go <- names(pV.classic)[pV.KSelim < pV.classic]
sel.go

cbind(termStat(sampleGO, sel.go),
      elim = pV.KSelim[sel.go],
      classic = pV.classic[sel.go])

# adjust for the p.values
res.dt$adj.p.value.classic <- pV.adj.classic[res.dt$GO.ID]
res.dt$adj.p.value.KS <- pV.adj.KS[res.dt$GO.ID]
res.dt$adj.p.value.elim <- pV.adj.KSelim[res.dt$GO.ID]

# plot
create.heatmap(
  x = res.dt[1:10, "adj.p.value.elim"],
  clustering.method = "none"
)
