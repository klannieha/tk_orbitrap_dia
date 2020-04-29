########## Script for combining dEPS and uEPS ############
library(data.table)
library(ggplot2)
library(tidyr)
library(VennDiagram)
#library(hrbrthemes)
setwd("D:/projects/pca_urine_spectral_lib/")
#basefolder <- "D:/projects/pca_urine_spectral_lib/"
########### load data ####################
libraries <- list.files("data/library", pattern = "se_filtered_spectrallib.tsv", full.names = T)
dEPS <- libraries[1] %>% fread()
uEPS <- libraries[2] %>% fread()

########## Compare and look at overlap #############
library(RColorBrewer)
myCol <- brewer.pal(8, "Pastel2")
# number of precursors intersect
intersect(dEPS$FullUniModPeptideName, uEPS$FullUniModPeptideName) %>% length()
# 34091 modified peptide overlap

dEPS.unique <- dEPS[!duplicated(dEPS, by = c("FullUniModPeptideName", "PrecursorCharge"))]
uEPS.unique <- uEPS[!duplicated(uEPS, by =c("FullUniModPeptideName", "PrecursorCharge"))]
EPS.overlap <- merge(dEPS.unique, uEPS.unique, by = c("FullUniModPeptideName", "PrecursorCharge"))
# 43564 overlapping precursors
EPS.overlap <- subset(EPS.overlap, select=c("FullUniModPeptideName", "PrecursorCharge"))
EPS.overlap$precursor_id <- 1:nrow(EPS.overlap)

dEPS.spec1 <- dEPS[FullUniModPeptideName == EPS.overlap$FullUniModPeptideName[1] & PrecursorCharge == EPS.overlap$PrecursorCharge[1]]
uEPS.spec1 <- uEPS[FullUniModPeptideName == EPS.overlap$FullUniModPeptideName[1] & PrecursorCharge == EPS.overlap$PrecursorCharge[1]]
dEPS.spec1$cohort <- "dEPS"
uEPS.spec1$cohort <- "uEPS"
spec1 <- rbind(dEPS.spec1, uEPS.spec1)
ggplot(spec1, aes(x = ProductMz, ymin=0, ymax = LibraryIntensity)) + geom_linerange() + facet_grid(~cohort)
p <- ggplot(spec1, aes(x= ProductMz) ) +
  geom_linerange(data = spec1[cohort == "dEPS"], 
                 aes(x = ProductMz, ymin=0, ymax = LibraryIntensity, color = myCol[1]),
                 show.legend = F, size = 1) +
  geom_linerange(data = spec1[cohort == "uEPS"],
                 aes(x = ProductMz, ymin = 0, ymax = -LibraryIntensity, color = myCol[2]),
                 show.legend = F, size = 1) +
  geom_hline(yintercept = 0) + ylab("Intensity") + 
  geom_label(aes(x=900, y=max(LibraryIntensity)*0.7, label="dEPS"), color=myCol[1]) +
  geom_label(aes(x = 900, y = -max(LibraryIntensity)*0.7), label = "uEPS", color = myCol[2]) +
  geom_label(aes(x = 700, y = max(LibraryIntensity)), label = unique(spec1$FullUniModPeptideName)) +
  geom_label(aes(x = 950, y = max(LibraryIntensity)), label = paste0(unique(spec1$PrecursorCharge), "+")) +
  theme_minimal(base_size = 16)
p

dEPS.specs <- merge(dEPS, EPS.overlap, by = c("FullUniModPeptideName", "PrecursorCharge"))
dEPS.specs$cohort <- "dEPS"
uEPS.specs <- merge(uEPS, EPS.overlap, by = c("FullUniModPeptideName", "PrecursorCharge"))
uEPS.specs$cohort <- "uEPS"
EPS.specs <- rbind(dEPS.specs, uEPS.specs)


dEPS.spec1 <- NULL
uEPS.spec1 <- NULL


## Output comparisons
pdf("D:/projects/pca_urine_spectral_lib/results/2020Apr_pca_combine/spectral_comparison.pdf", onefile = T)

precursors <- EPS.overlap$precursor_id
for(p in precursors){
  sub <- EPS.specs[precursor_id == p]
  p <- ggplot(sub, aes(x= ProductMz) ) +
    geom_linerange(data = sub[cohort == "dEPS"], 
                   aes(x = ProductMz, ymin=0, ymax = LibraryIntensity, color = myCol[1]),
                   show.legend = F, size = 1) +
    geom_linerange(data = sub[cohort == "uEPS"],
                   aes(x = ProductMz, ymin = 0, ymax = -LibraryIntensity, color = myCol[2]),
                   show.legend = F, size = 1) +
    geom_hline(yintercept = 0) + ylab("Intensity") + 
    geom_label(aes(x=max(ProductMz)*0.9, y=max(LibraryIntensity)*0.7, label="dEPS"), color=myCol[1]) +
    geom_label(aes(x = max(ProductMz)*0.9, y = -max(LibraryIntensity)*0.7), label = "uEPS", color = myCol[2]) +
    geom_label(aes(x = max(ProductMz)*0.7, y = max(LibraryIntensity)), label = unique(spec1$FullUniModPeptideName)) +
    geom_label(aes(x = max(ProductMz)*0.7, y = max(LibraryIntensity)*0.85), label = paste0(unique(spec1$PrecursorCharge), "+")) +
    theme_minimal(base_size = 16)
  print(p)
}
dev.off()

############# Calculate spectral simularity by dot product ####################

############# Example ################
max(spec1$ProductMz)
min(spec1$ProductMz)

# create bins
bins <- seq(min(spec1$ProductMz), max(spec1$ProductMz), by = 1)

# assign bins
binned_spectra <- findInterval(spec1$ProductMz, bins, rightmost.closed = T)
spec1[, bins := binned_spectra]

# assign bin weight
bin_weights <- lapply(seq_along(bins), function(x){
  b <- x
  if(b %in% spec1$bins){
    w <- max(spec1[bins == b]$LibraryIntensity) / max(spec1$LibraryIntensity) * 100
  } else{
    w <- 1
  }
  w
}) %>% unlist()

# bin spectra
s1 <- lapply(seq_along(bins), function(x){
  b <- x
  spec <- spec1[cohort == "uEPS"]
  if(b %in% spec$bins){
    int <- spec[bins == b]$LibraryIntensity / max(spec$LibraryIntensity) * 100
  }else{
    int <- 0
  }
  int
}) %>% unlist()

s2 <- lapply(seq_along(bins), function(x){
  b <- x
  spec <- spec1[cohort == "dEPS"]
  if(b %in% spec$bins){
    int <- spec[bins == b]$LibraryIntensity / max(spec$LibraryIntensity) * 100
  }else{
    int <- 0
  }
  int
}) %>% unlist()

# s1 is a vector of 889 bins with library intensity
# s2 is a vector of 889 bins with library intensity

# Calculate dot product
s1 <- s1 * bin_weights

dp <- s1 %*% s2
absmag <- sqrt(sum(s1^2)) * sqrt(sum(s2^2))

dp <- dp/absmag

################# Create functions #######################

Createbins <- function(allSpecdf, pid){
  # Creates data.frame of bins for one peptide precursor
  specs <- allSpecdf[precursor_id == pid]
  bins <- seq(round(min(specs$ProductMz)), round(max(specs$ProductMz)), by = 0.005)
  binned_spectra <- findInterval(specs$ProductMz, bins, rightmost.closed = T)
  specs[, bin := binned_spectra]
  bin_weights <- lapply(seq_along(bins), function(x){
    b <- x
    if(b %in% specs$bin){
      w <- max(specs[bin == b]$LibraryIntensity) #/ max(specs$LibraryIntensity) * 100
    } else{
      w <- 1
    }
    w
  }) %>% unlist()
  bins <- data.frame(bins = bins, weights = bin_weights, stringsAsFactors = F)
  return(bins)
}

binSpectrum <- function(spectrumdf, bins){
  # spectrumdf dataframe/data.table  peptide-fragment info
  # bin  vector   list of bins
  # formatted with OpenMS
  binned_spectra <- findInterval(spectrumdf$ProductMz, bins, rightmost.closed = T)
  spectrumdf[, binned := binned_spectra]
  s <- lapply(seq_along(bins), function(x){
    b <- x
    if(b %in% spectrumdf$binned){
      int <- spectrumdf[binned == b]$LibraryIntensity# / max(spectrumdf$LibraryIntensity) * 100
    }else{
      int <- 0
    }
    int
  }) %>% unlist()
  return(s)
}

NormDotP <- function(spec1, spec2, weights){
  # spec1 & spec2 are vectors of equal length
  # weights  vector     bin weights
  spec1 <- spec1 * weights
  dp <- spec1 %*% spec2
  absmag <- sqrt(sum(spec1^2)) * sqrt(sum(spec2^2))
  dp <- dp/absmag
  return(dp[1])
}


test <- EPS.specs[precursor_id == 34316]
test.bins <- Createbins(EPS.specs, 34316)
test.dEPS <- test[cohort == "dEPS"]
test.uEPS <- test[cohort == "uEPS"]

test.uspec <- binSpectrum(test.uEPS, test.bins$bins)
test.dspec <- binSpectrum(test.dEPS, test.bins$bins)

test.dp <- NormDotP(test.uspec, test.dspec, test.bins$weights)

########### Run bulk ###############################

NormalizedDotProduct <- numeric(length(precursors))

pdf("D:/projects/pca_urine_spectral_lib/results/2020Apr_pca_combine/spectral_comparison_dotproduct.pdf", onefile = T)

precursors <- EPS.overlap$precursor_id
for(p in precursors){
  sub <- EPS.specs[precursor_id == p]
  sub.bins <- Createbins(EPS.specs, p)
  sub.dEPS <- sub[cohort == "dEPS"]
  sub.uEPS <- sub[cohort == "uEPS"]
  sub.uspec <- binSpectrum(sub.uEPS, sub.bins$bins)
  sub.dspec <- binSpectrum(sub.dEPS, sub.bins$bins)
  NormDP <- NormDotP(sub.uspec, sub.dspec, sub.bins$weights)
  NormalizedDotProduct[p] <- NormDP
#  p <- ggplot(sub, aes(x= ProductMz) ) +
#    geom_linerange(data = sub.dEPS, 
#                   aes(x = ProductMz, ymin=0, ymax = LibraryIntensity, color = myCol[1]),
#                   show.legend = F, size = 1) +
#    geom_linerange(data = sub.uEPS,
#                   aes(x = ProductMz, ymin = 0, ymax = -LibraryIntensity, color = myCol[2]),
#                   show.legend = F, size = 1) +
#    geom_hline(yintercept = 0) + ylab("Intensity") + 
#    geom_label(aes(x=max(ProductMz)*0.9, y=max(LibraryIntensity)*0.7, label="dEPS"), color=myCol[1]) +
#   geom_label(aes(x = max(ProductMz)*0.9, y = -max(LibraryIntensity)*0.7), label = "uEPS", color = myCol[2]) +
#    geom_label(aes(x = max(ProductMz)*0.7, y = max(LibraryIntensity)), label = unique(sub$FullUniModPeptideName)) +
#    geom_label(aes(x = max(ProductMz)*0.1, y = max(LibraryIntensity)*0.85), label = paste0(unique(spec1$PrecursorCharge), "+")) +
#    theme_minimal(base_size = 16) +
#    geom_label(aes(x = max(ProductMz)*0.9, y = -max(LibraryIntensity)*0.8), label = NormDP, color = myCol[3], size = 2)
#  print(p)
}
dev.off()

