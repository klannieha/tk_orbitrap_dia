############## Analysis and Comparison of Libraries #############
library(data.table)
library(ggplot2)
library(tidyr)
library(ggExtra)

source("D:/source/tk_orbitrap_dia/tk_orbitrap_dia/analysis/Utilities.R")
################### Load data #########################

basefolder <- "D:/projects/pca_urine_spectral_lib/"

nonlinear_lib <- list.files(paste0(basefolder, "data/irt"), pattern = "nonlinear_assaylib.tsv", full.names = T)
nonlinear_lib <- fread(nonlinear_lib)
# From riun E156

postDRE_nonlinear_lib <- list.files(paste0(basefolder, "data/irt"), pattern = "osw_assaylib.tsv", full.names = T) %>% fread()

pca_eps.msms <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "msms.txt", full.names = T)
pca_eps.evidence <- list.files(paste0(basefolder, "/data/pca_dda/pca_directEPS"), pattern = "evidence.txt", full.names = T)

irt <- list.files(paste0(basefolder, "/data/irt"), pattern="irt_assay", full.names = T)


pca_eps.msms <- fread(pca_eps.msms)
pca_eps.evidence <- fread(pca_eps.evidence)
pca_eps.proteinGroups <- fread(pca_eps.proteinGroups) 
pca_eps.proteinGroups <- subset(pca_eps.proteinGroups,select=c("Protein IDs", "Majority protein IDs", "Gene names"))
irt <- fread(irt)

irt <- subset(irt, select = c("ModifiedPeptideSequence", "NormalizedRetentionTime"))
setnames(irt, "NormalizedRetentionTime", "iRT")


irt <- irt[!duplicated(irt)]


runs <- pca_eps.msms$`Raw file` %>% unique()

ev <- subset(pca_eps.evidence, select=c("id", "Max intensity m/z 0"))
setnames(ev, "Max intensity m/z 0", "Intensity")
ms <- merge(pca_eps.msms, ev, by.x = "Evidence ID", by.y = "id")

spectral_lib <- list.files(paste0(basefolder, "data/library"), pattern = "nonlinear_lib.tsv", full.names = T)
spectral_lib <- fread(spectral_lib[1])
spectral_lib <- merge(spectral_lib, pca_eps.proteinGroups, by.x = "Proteins", by.y = "Majority protein IDs")


urine_lib <- list.files(paste0(basefolder, "data/library"), pattern = "urine_nonlinear.tsv", full.names = T)
urine_lib <- fread(urine_lib[1])

urine.proteinGroups <- list.files(paste0(basefolder, "/data/pca_dda/pca_postDRE_urine/txt_20200220_1-200_og"), 
                                  pattern = "proteinGroups.txt", full.names = T) 
urine.proteinGroups <- fread(urine.proteinGroups)

urine.proteinGroups <- subset(urine.proteinGroups,select=c("Protein IDs", "Majority protein IDs", "Gene names"))

urine_lib <- merge(urine_lib, urine.proteinGroups, by.x = "Proteins", by.y = "Majority protein IDs")

irtlib_count <- postDRE_nonlinear_lib[!duplicated(postDRE_nonlinear_lib,
                                                  by = c("ModifiedPeptideSequence", "PrecursorCharge"))]

nrow(irtlib_count)
# 17824

irtlib_count$PeptideSequence %>% unique() %>% length()
# 14255

irtlib_count$ProteinId %>% unique() %>% length()
# 3458

count_df <- data.frame(variable=c("Precursors", "Peptides", "ProteinGroups"),
                       counts=c(nrow(irtlib_count), 
                                                                                      length(unique(irtlib_count$PeptideSequence)),
                                                                                      length(unique(irtlib_count$ProteinId))),
                       cohort = "postDRE_urine")
count_df <- rbind(count_df, count)

##################### Analysis ############################

########check the spectra of library########

head(spectral_lib)

# Example peptide 1

tr1 <- spectral_lib[transition_group_id == "213776"]

# look at spectra

g <- ggplot(tr1, aes(x = ProductMz, ymin=0, ymax = LibraryIntensity)) + geom_linerange()

g + theme_bw() + ggtitle(tr1$FullUniModPeptideName %>% unique())

tr_ids <- spectral_lib$transition_group_id %>% unique()

# random sample 10

tr <- sample(tr_ids, 10)

for (transition in tr) {
  df <- spectral_lib[transition_group_id == transition]
  g <- ggplot(df, aes(x = ProductMz, ymin=0, ymax = LibraryIntensity)) + geom_linerange()
  g <- g + theme_bw() + ggtitle(transition)
  ggsave(paste0("library_spectra_", transition, ".png"))
}


######### Analyze score distribution ############

# Things to do
# Extract the information of peptide that were taken for library
# Compare the score and RT distribution 
# identify if they are found in more than 1 run

peptides <- spectral_lib$FullUniModPeptideName %>% unique()

ms$`Modified sequence` <- reformat_mods(ms$`Modified sequence`)

ms <- subset(ms, select = c("Raw file", "Sequence", "Modified sequence", "Proteins",
                            "Protein Names", "Charge", "Gene Names", "m/z", "Mass error [ppm]",
                            "Retention time", "PEP", "Score", "Delta score", "Peak coverage", "Intensity"))
setnames(ms, colnames(ms), c("raw", "sequence", "ModifiedSequence", "proteins", "ProteinNames", 
                             "Charge", "GeneNames", "mz", "ppmdiff", "RetentionTime", "pep", "score",
                             "deltaScore", "peakCoverage", "Intensity"))

ms <- ms[ModifiedSequence %in% peptides]
ms <- ms[order(pep)]
ms_unique <- ms[!duplicated(ms, by = c("Charge", "ModifiedSequence"))]

g <- ggplot(ms_unique, aes(x = score, y = RetentionTime, color = raw)) + 
  geom_point(size=1, alpha=0.5, show.legend = F) +
  scale_colour_viridis_d(option = "inferno")
  
ggMarginal(g, type="density")

g <- ggplot(ms_unique, aes(x = pep, y = score, color = raw)) +
  geom_point(size=1, alpha=0.5, show.legend = F) +
  scale_colour_viridis_d(option = "inferno") +
  theme_bw()
ggMarginal(g, type = "density")


g <- ggplot(ms_unique, aes(x = mz, y = ppmdiff, color = raw)) +
  geom_point(size=1, alpha=0.5, show.legend = F) +
  scale_color_viridis_d(option="inferno") + 
  theme_bw()
g <- g + geom_hline(yintercept = c(5, -5), show.legend = F, linetype = "dashed", color = "grey")
ggMarginal(g, type = "density", margins = "y", fill = sample(viridis(10),1))

#ms_unique$Intensity <- as.numeric(ms_unique$Intensity)
ms_unique$logInt <- log(ms_unique$Intensity)
ggplot(ms_unique, aes(x = RetentionTime, y = mz, color = Intensity)) + 
  geom_point(size=0.6) +
  theme_bw() +
  scale_color_gradientn(colors = viridis(10, option = "inferno"))


spectral_lib$UniProtID %>% unique() %>% head()

######### Compare libraries ##########################
spectral_lib <- spectral_lib[!duplicated(spectral_lib, by = c("ModifiedPeptideSequence", "Charge"))]
eps.count <- data.frame(variable = c("Precursors", "Peptides", "Proteins", "Genes"),
                        counts = c(nrow(spectral_lib), spectral_lib$PeptideSequence %>% unique() %>% length(),
                                   spectral_lib$`Protein IDs` %>% unique() %>% length(),
                                   spectral_lib$`Gene names` %>% unique() %>% length()),
                        cohort = "directEPS")
urine_lib <- urine_lib[!duplicated(urine_lib, by = c("ModifiedPeptideSequence", "Charge"))]

urine.count <- data.frame(variable = c("Precursors", "Peptides", "Proteins", "Genes"),
                          counts = c(nrow(urine_lib), urine_lib$PeptideSequence %>% unique() %>% length(),
                                     urine_lib$`Protein IDs` %>% unique() %>% length(),
                                     urine_lib$`Gene names` %>% unique() %>% length()),
                          cohort = "postDRE_urine")
counts <- rbind(eps.count, urine.count)
counts$variable <- factor(counts$variable, levels = c("Precursors", "Peptides", "Proteins", "Genes"))
counts$cohort <- factor(counts$cohort, levels = c("directEPS", "postDRE_urine"))
counts <- as.data.table(counts)
#counts[variable == "ProteinGroups"]$variable <- "Proteins"
ggplot(counts[variable %in% c("Peptides", "Proteins")], aes(x = variable, y = counts, fill = cohort, label = counts)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), vjust = -.9 ) +
  theme_classic(base_size = 16) + scale_fill_viridis_d(option = "plasma", end = 0.5, direction = -1) +
  xlab("")

peptides <- union(spectral_lib$PeptideSequence %>% unique(), urine_lib$PeptideSequence %>% unique()) %>% unique()
urine.peptides <- urine_lib$PeptideSequence %>% unique()
eps.peptides <- spectral_lib$PeptideSequence %>% unique()

pep_in_both <- intersect(urine.peptides, eps.peptides)
urine_only <- urine.peptides[- which(urine.peptides %in% pep_in_both)]
eps_only <- eps.peptides[-which(eps.peptides %in% pep_in_both)]

pep_in_both <- data.frame(Peptide = pep_in_both, Identified = "Both", stringsAsFactors = F)
urine_only <- data.frame(Peptide = urine_only, Identified = "postDRE_urine", stringsAsFactors = F)
eps_only <- data.frame(Peptide = eps_only, Identified = "directEPS", stringsAsFactors = F)

ID <- rbind(pep_in_both, urine_only, eps_only)


total_count <- data.frame(variable = rep(c("Peptide", "Protein", "Genes"), each = 3),
                          Identified = rep(c("Both", "directEPS", "postDRE_urine"), times = 3),
                          counts = numeric(9)
                          ) %>% as.data.table()

total_count[variable == "Peptide"]$counts <- table(ID$Identified)

proteins <- union(spectral_lib$`Protein IDs` %>% unique(), urine_lib$`Protein IDs` %>% unique()) %>% unique()
urine.proteins <- urine_lib$`Protein IDs` %>% unique()
eps.proteins <- spectral_lib$`Protein IDs` %>% unique()

protein_in_both <- intersect(urine.proteins, eps.proteins)
urine_only <- urine.proteins[- which(urine.proteins %in% protein_in_both)]
eps_only <- eps.proteins[-which(eps.proteins %in% protein_in_both)]

protein_in_both <- data.frame(Protein = protein_in_both, Identified = "Both", stringsAsFactors = F)
urine_only <- data.frame(Protein = urine_only, Identified = "postDRE_urine", stringsAsFactors = F)
eps_only <- data.frame(Protein = eps_only, Identified = "directEPS", stringsAsFactors = F)

ID <- rbind(protein_in_both, urine_only, eps_only)
table(ID$Identified)
protein_in_both.canonical <- lapply(protein_in_both, function(x){
  if(grepl(";", x)){
  n <- strsplit(x, ";")[[1]]
  n <- n[-grep("-", n)]
  n} else{ x }
}) %>% unlist()

write(protein_in_both.canonical, file = "D:/projects/pca_urine_spectral_lib/results/proteins_common.txt", sep = "\n")

total_count[variable == "Protein"]$counts <- table(ID$Identified)

genes <- union(spectral_lib$`Gene names` %>% unique(), urine_lib$UniProtID %>% unique()) %>% unique()
urine.genes <- urine_lib$`Gene names` %>% unique()
eps.genes <- spectral_lib$`Gene names` %>% unique()
"FN1" %in% gene_in_both
gene_in_both <- intersect(urine.genes, eps.genes)
urine_only <- urine.genes[- which(urine.genes %in% gene_in_both)]
eps_only <- eps.genes[-which(eps.genes %in% gene_in_both)]

gene_in_both <- data.frame(Gene = gene_in_both, Identified = "Both", stringsAsFactors = F)
urine_only <- data.frame(Gene = urine_only, Identified = "postDRE_urine", stringsAsFactors = F)
eps_only <- data.frame(Gene = eps_only, Identified = "directEPS", stringsAsFactors = F)

ID <- rbind(gene_in_both, urine_only, eps_only)
total_count[variable == "Genes"]$counts <- table(ID$Identified)
total_count$variable <- factor(total_count$variable, levels = c("Peptide", "Protein", "Genes"))
setnames(total_count, "Identified", "Detected")
ggplot(total_count, aes(x = variable, y = percent, label= percent, fill = Detected)) + 
  geom_bar(stat = "identity", position = position_stack()) + 
  geom_text(position= position_stack(), vjust = 1.8, color = "white") +
  scale_fill_viridis_d(end = 0.5, option = "plasma", direction = -1) + theme_classic(base_size = 16) + xlab("")

########### Venn Diagram #################

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
V = venn.diagram(
  x = list(urine = urine.genes, EPS =eps.genes),
  category.names = c('postDRE_urine', "directEPS"),
  filename = NULL,
  output=FALSE,
  # Output features
  #imagetype="png" ,
  #height = 480 , 
  #width = 480 , 
  #resolution = 300,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #rotation = 1,
  ext.text = FALSE
)

library(VennDiagram);
library(grid)

grid.newpage()
grid.draw(V)

cpcgene <- fread("D:/projects/pca_urine_spectral_lib/data/pca_dda/Ankt_Tissue_results/1-s2.0-S153561081930100X-mmc3.txt")
cpcgene_uni <- cpcgene$Gene %>% unique()
V = venn.diagram(
  x = list(urine = urine.genes, EPS =eps.genes, cpcgene_uni),
  category.names = c('postDRE_urine', "directEPS", "Tissue"),
  filename = NULL,
  output=FALSE,
  # Output features
  #imagetype="png" ,
  #height = 480 , 
  #width = 480 , 
  #resolution = 300,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 153),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  #rotation = 1,
  ext.text = FALSE
)
grid.newpage()
grid.draw(V)

all_genes_in_both <- lapply(gene_in_both, function(x){
  g <- strsplit(x, ";")[[1]]
  g
})
all_genes_in_both <- unlist(all_genes_in_both)
write(all_genes_in_both, file = "D:/projects/pca_urine_spectral_lib/results/genes_common.txt", sep = "\t")

########### protein canonical mapped ###############
# load the libraries that are mapped to canonical proteins
lib_urine <- fread("D:/projects/pca_urine_spectral_lib/data/library/urine_protein_mapped_20200402.tsv")
lib_eps <- fread("D:/projects/pca_urine_spectral_lib/data/library/directEPS_protein_mapped_20200402.tsv")

common_proteins <- intersect(lib_urine$ProteinName, lib_eps$ProteinName) %>% unique()
write(common_proteins, file = "D:/projects/pca_urine_spectral_lib/results/proteins_mapped_common.txt", sep = "\t")
