#=======
# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
library(VennDiagram)
library(Matrix)
library(eulerr)
wd <- "/Users/annieha/OneDrive - University of Toronto/Attachments/"
src_wd <- paste0(wd, "source/tk_orbitrap_dia/analysis")
source(paste0(wd, "projects/pca_uEPS1_spectral_lib/src/Amanda_source.R"))

#################### Functions for Use ###################

# Median normalization
medNorm <- function(x){
  med <- apply(x, 2, median, na.rm=T)
  av.med <- mean(med)
  norm.factors.med <- av.med/med 
  y <- as.data.frame(t(t(x) * norm.factors.med))
  row.names(y) <- row.names(x)
  return(y)
}

getCV <- function(x){
  CV = apply(x, 1, function(y) 100*sd(y, na.rm = T)/mean(y, na.rm = T))
  mean = apply(x, 1, function(y) mean(y, na.rm = T))
  x$CV <- CV
  x$mean <- mean
  return(x)
}
# Reformat tpp library
GeneFormat <- function(x){
  if(grepl("Subgroup", x)){
    p <- unlist(strsplit(x, "\\|"))
    p <- p[grepl("^[A-Z][0-9](.*)", p)]
    p <- gsub("-.*", "", p)
    p <- unique(p)
    p <- sort(p)
    p <- paste0(p, collapse=";")
  } else{
    p <- strsplit(x, ";")[[1]]
    p <- gsub("-.*", "", p)
    p <- unique(p)
    p <- sort(p)
    p <- paste0(p, collapse=";")
  }
  return(p)
}
a <- "Q71U36;Q71U36-2;P68363;P68363-2;Q13748;Q13748-2;Q6PEY2;P68366;P68366-2"
a <- "Subgroup_0_2/DECOY_sp|Q00325|MPCP_HUMAN/DECOY_sp|Q00325"
a <- "Subgroup_1_5/sp|Q56UQ5|TPT1L_HUMAN/sp|P13693"
GeneFormat(a)  


# filtering of data

AllPeptides_filter <- function(Peptidesdf){
  #filters out reverse and potential contaminents
  df <- Peptidesdf
  df <- df[Reverse != "+" | `Potential contaminant` != "+"]
  return(df)
}


peptideCountOverObservations <- function(peptidesdf){
  select_cols <- colnames(peptidesdf)[grep("Identification(.*)", colnames(peptidesdf))]
  df <- subset(peptidesdf, select = select_cols)
  counts <- data.frame(Peptide = peptidesdf$Sequence, Found = numeric(length(peptidesdf$Sequence)),
                       Matched = numeric(length(peptidesdf$Sequence)))
  s <- lapply(1:nrow(df), function(x){
    r <- df[x,]
    found <- r == "By MS/MS"
    found <- sum(found)
    matched <- r == "By matching"
    matched <- sum(matched)
    c(found, matched)
  })
  s <- do.call('rbind', s) %>% as.data.frame(stringsAsFactors=F) %>% as.data.table()
  colnames(s) <- c("Found", "Matched")
  counts$Found <- s$Found
  counts$Matched <- s$Matched
  return(counts)
}


CommonPeptides <- function(msms_lst){
  # list of dataframe from multiple runs
  runs <- 1:length(msms_lst)
  lst <- lapply(runs, function(x){
    run <- msms_lst[[x]]
    run <- run[order(PEP)]
    sequence <- run['Modified sequence']
    sequence <- unique(sequence)
    sequence <- reformat_mods(sequence)
    sequence})
  common <- Reduce(intersect, lst)
  return(common)
}
  
# DIA runs
# Peptide detection / reproducibility:
Peptide_Detection <- function(df, annotation, n){
  # df - the annotated pyprophet
  # n  - number of replicates
  # annotation - df of experiments, filenames, condition and replicates
  methods <- unique(annotation$Condition)
  detection <- lapply(methods, function(x){
    r <- df[Condition == x]
    r <- r[!duplicated(r, by=c("Sequence", "Replicate"))]
    r <- dcast(r, Sequence~Replicate, value.var = "Intensity")
    count <- r %>% is.na() %>% rowSums()
    r$Detection <- count
    r$Method <- x
    r
  })
  detection <- do.call('rbind', detection) %>% as.data.table()
  detection$Detection <- n - detection$Detection
  return(detection) 
}



# Median normalization
medNorm <- function(x){
  med <- apply(x, 2, median, na.rm=T)
  av.med <- mean(med)
  norm.factors.med <- av.med/med 
  y <- as.data.frame(t(t(x) * norm.factors.med))
  row.names(y) <- row.names(x)
  return(y)
}

# Need a function to merge list of dataframes and get mean or sd or na


CommonPeptideStats <- function(msms_lst, run_ids, column,  mean = T, sd = T, na = T){
  "
  msms_lst	list of dataframes from msms files first column contains peptide sequence
  run_ids 	list of run ids to replace for column names
  column	string of the column that is interested in analyzing e.g. RT or Intensity
  "
  df <- lapply(1:length(run_id), function(x){
    rep <- subset(msms_lst[[x]], select=c("PeptideSequence", column))
    setnames(rep, column, run_ids[x])
    rep})
  df <- Reduce(function(x, y) merge(x,y, by = "PeptideSequence", all=TRUE, df))
  if(mean){
    mean <- lapply(1:nrow(df), function(x){
      m <- df[x, 2:ncol(df)] %>% as.numeric() %>% mean(na.rm=T)
      m}) %>% unlist()
    df$mean <- mean
  }else if(sd){
    sd <- lapply(1:nrow(df), function(x){
      s <- df[x, 2:ncol(df)] %>% as.numeric() %>% sd(na.rm=T)
      s}) %>% unlist()
    df$sd <- sd
  }else if(na){
    na <- lapply(1:nrow(df), function(x){
      n <- df[x, 2:ncol(df)] %>% is.na() %>% sum()
      n}) %>% unlist()
    df$na <- na}
  return(df)
}

reduce_OpenSWATH_output <- function(data, column.names=NULL){
  if(is.null(column.names)){
    column.names <- c('ProteinName', 'FullPeptideName', 'Sequence', 'Charge', 'aggr_Fragment_Annotation', 'aggr_Peak_Area', 'filename', 'm_score', 'decoy', "Intensity", "RT", "run_id", "transition_group_id")
  }
  if(length(column.names) > length(column.names[column.names %in% colnames(data)])){
    col.names.missing <- column.names[!column.names %in% colnames(data)]
    warning("These columns are missing from the data:", paste(unlist(col.names.missing), collapse=", "))
    
  }
  # Keep only required columns for MSStats and mapDIA
  if(length(column.names) == length(column.names[column.names %in% colnames(data)])){
    data.filtered <- data[, ..column.names]
    return(data.filtered)
  }
}

# Need to right a function to automate using pairwise comparison for finding a reference run


filter_pyprophet <- function(pyprophet_out, fdr = 0.01, Decoy = 0, peak_group = 1){
  # Filter the pyprophet output by m_score, decoy and peak groups
  tmp <- copy(pyprophet_out)
  tmp <- tmp[decoy == Decoy]
  tmp <- tmp[m_score <= fdr]
  tmp <- tmp[peak_group_rank == peak_group]
  tmp <- tmp[!duplicated(tmp, by=c("filename", "FullPeptideName", "Charge"))]
  return(tmp)
}

annotate_pyprophet <- function(pyprophet_out){
  df <- data.frame(filename = unique(pyprophet_out$filename),
                   stringsAsFactors = F)
  method <- lapply(basename(df$filename), function(x){
    n <- gsub("\\.(.*)", "", x)
    n <- strsplit(n, "_")[[1]]
    m <- paste0(n[grepl("[A-Za-z]", n)], collapse = " ")
    r <- n[length(n)]
    c(m, r)
  })
  method <- do.call('rbind', method)
  df$Condition <- method[,1]
  df$Replicate <- method[,2]
  return(df)
}

save_svg <- function(plot, file ,show.plot = F, width = 6, height = 8, unit = "in"){
  if(show.plot){
    plot
  }
  ggsave(plot, filename = file, width = width, height = height, unit = unit,
         device = 'svg')
}

if(!exists("runs_annotation")){
  uEPS.names <- fread(paste0(wd, "projects/pca_uEPS1_spectral_lib/data/pca_dda/uEPS1_names.txt"))
  uEPS.names$Sample <- str_extract(uEPS.names$Old_name, "S[0-9]+")
  uEPS.names$SampleID <- str_extract(uEPS.names$New_name, "UP[0-9]+")
  runs_annotation <- fread(paste0(wd, "projects/pca_uEPS1_spectral_lib/data/clinical/annotation.tsv"))
  runs_annotation <- merge(uEPS.names, runs_annotation, by.x = "Sample",
                           by.y = "Condition", all.x = T)
  runs_annotation$Batch <- gsub("[^0-9]", "", runs_annotation$Batch)
}

#MQ data extraction 
readMQ <- function(filename){
  #' This funciontion is to read MQ dataframes
  #' filename    character path
  #' Read only 1 file at a time
  print(paste("Reading file:", filename, sep = " "))
  dataframe <- read.table(filename, header = T, sep = "\t")
  return(dataframe)
}

getMQsubset <- function(dataframe, proteinGroup = T, peptides = F, column.prefix = "Intensity\\."){
  #' This function is to extract the subset columns of the MQ output
  #' Converts the wide dataframe to long format
  #' dataframe  Inputfile usually peptides.txt or proteinGroups.txt
  newName <- gsub("\\.(.*)", "", column.prefix)
  if(proteinGroup){
    print("Reading proteinGroups dataframe")
    df.Int <- subset(dataframe, 
                     select = c("Gene.names", "Protein.IDs", "Protein.names",
                                colnames(dataframe)[grepl(column.prefix,
                                                          colnames(dataframe))]))
    df.Int <- melt(df.Int, 
                   id.vars = c("Gene.names", "Protein.IDs", "Protein.names"),
                   variable.name = "Sample", value.name =  newName)
  }else if(peptides){
    print("Reading peptides dataframe")
    df.Int <- subset(dataframe,
                     select = c("Sequence",
                                colnames(dataframe)[grepl(column.prefix,
                                                          colnames(dataframe))]))
    df.Int <- melt(df.Int, id.vars = "Sequence", variable.name = "Sample", 
                   value.name = newName)
  }
  df.Int$Sample <- gsub(column.prefix, "", df.Int$Sample)
  return(df.Int)
}

# read iRTs
wd <- "/Users/annieha/OneDrive - University of Toronto/Attachments/"
# 
iRTs <- list.files(paste0(wd, "/projects/pca_uEPS1_spectral_lib/data"), pattern = "irt", full.names = T)
iRTs
iRTs <- readxl::read_xls(iRTs, sheet = 2)

iRTs <- iRTs %>%
 filter(rank == 1) %>%
 select(`Q1 monoisotopic`, `precursor charge`, `nominal sequence`, `sequence id` ) %>%
 distinct_all()
# 
colnames(iRTs) <- c("PrecursorMz", "PrecursorCharge", "Sequence", "PeptideID")

# perform wilcoxon / Mann-Whitney U test per protein
Htest <- function(protein_id,variable, category1, category2, key, df, test = 'utest', intensity_col){
  #' protein_id   chr   protein accession of interest
  #' category1    chr   group 1 category e.g. ISUP 1
  #' category2    chr   group 2 category e.g. ISUP 2
  #' key          chr   group variable
  #' df        dataframe  long dataframe containing at least columns with protein_id and key
  #' test         chr   statistical tests to run e.g. 'utest', 'ttest'
  #
  p <- protein_id
  #print(p)
  g1 <- df %>% filter({{key}} == category1 & {{variable}} == p) %>% pull({{intensity_col}}) %>% as.numeric()
  g2 <- df %>% filter({{key}} == category2 & {{variable}} == p) %>% pull({{intensity_col}}) %>% as.numeric()
  if(test == 'utest'){
    Ht <- wilcox.test(g1, g2, paired = F)
  }else if(test == 'ttest'){
    Ht <- t.test(g1, g2, paired = F, alternative = "two.sided")
  }
  return(Ht)
}

# test
#Htest(p, "1", ">1", ClinicISUP, pg.de.impute, "utest", logTransform = T)

# calculate log2FC

calculatePeptidelogFC <- function(data, colname, group1, group2, level){
  #' data  long dataframe with Log2Intensity
  #' variable column name
  #' group1 group 1 from the column
  #' group2 group2 from the column
  #' level    mean or median
  FC <- data %>%
    group_by({{colname}}, Sequence) %>%
    summarise(Intensity = ifelse(level == "mean", mean(Log2Intensity, na.rm = T), median(Log2Intensity, na.rm = T)), 
              Observation = n()) %>% 
    pivot_wider(id_cols = "Sequence", names_from = {{colname}}, values_from = c("Intensity", "Observation")) %>%# head()
    mutate(FC = {{group2}} - {{group1}}) %>% 
    as.data.table()
  return(FC)
}

#calculatePeptidelogFC(peptides.top3.lnorm, Grade, `1`, `2+`)

calculateProteinlogFC <- function(data, colname, group1, group2, intensity_col, level = "mean", paired = F){
  #' data  long dataframe with Log2Intensity
  #' variable column name
  #' group1 group 1 from the column
  #' group2 group2 from the column
  stat.res <- data %>% 
    select({{intensity_col}}, {{colname}}, Gene.names) %>% 
    mutate(int = {{intensity_col}}, 
           group = {{colname}}) %>% 
    group_by(Gene.names) %>% 
    summarise(pvalue = wilcox.test(int ~ group, paired = paired, exact = F)$p.value) %>% 
    as.data.table()
  
  FC <- data %>%
    group_by({{colname}}, Gene.names) %>%
    summarise(Intensity = ifelse(level == "mean", mean({{intensity_col}}, na.rm = T), median({{intensity_col}}, na.rm = T))) %>% 
    pivot_wider(id_cols = "Gene.names", names_from = {{colname}}, values_from = "Intensity") %>%# head()
    mutate(FC = {{group2}} - {{group1}}) %>% 
    as.data.table()
  
  FC <- merge(FC, stat.res, by = "Gene.names")
  FC$FDR <- p.adjust(FC$pvalue, method = "fdr")
  
  return(FC)
}


plot_PeptideVolcano <- function(data, label = F, pvalue_column,padj_column, FC_column, sig_level = 0.05, legends = F){
  #' data    dataframe with FC and p values
  #' label   boolean   include labels or not
  #' add other legends if needed
  df <- copy(data)
  
  df <- df %>% mutate(Sig = ifelse({{padj_column}} < sig_level, TRUE, FALSE)) %>%
    mutate(Enriched = case_when((Sig & {{FC_column}} > 0) ~ "UP",
                                (Sig & {{FC_column}} < 0) ~ "DOWN"))
  p <- df %>%
      mutate(`-log10(p-value)` = -log10({{pvalue_column}})) %>%
      ggplot(aes(x = {{FC_column}}, y = `-log10(p-value)`)) +
      geom_point(aes(color = Enriched), show.legend = legends) +
      scale_color_manual(breaks = c("UP", "DOWN"), values = c("red", "blue")) +
      plot_theme()

  if(label){
    df <- df %>% mutate(LABEL = ifelse(Sig, paste(peptide_sequence, GeneName, sep = " * "), NA))
    p <- p + geom_label_repel(aes(label = LABEL), max.overlaps = 30)
  }
  p
  return(p)
}
 

plot_PGVolcano <- function(data, label = F, label_col, pvalue_column,padj_column, 
                           FC_column, sig_level = 0.05, legends = F, colours = c("UP" = "red", "DOWN" = "blue")){
  #' data    dataframe with FC and p values
  #' label   boolean   include labels or not
  #' add other legends if needed
  df <- copy(data)
  df <- df %>% mutate(Sig = ifelse({{padj_column}} < sig_level, TRUE, FALSE)) %>%
    mutate(Enriched = case_when((Sig & {{FC_column}} > 0) ~ "UP",
                                (Sig & {{FC_column}} < 0) ~ "DOWN")) %>% 
    as.data.table()
  p <- df %>%
    ggplot(aes(x = {{FC_column}}, y = -log10({{pvalue_column}}))) +
    geom_point(aes(color = Enriched), show.legend = legends) +
    scale_color_manual(breaks = c("UP", "DOWN"), values = colours) +
    theme_ez()
  
  if(label==TRUE){
    setnames(df, label_col, "Gene.names")
    #df <- df %>% mutate(LABEL = ifelse(Sig == TRUE, Gene.names, NA))
    p <- p + geom_label_repel(data = df[Sig == TRUE], aes(label = Gene.names), max.overlaps = Inf)
  }
  p
  return(p)
}


spot.size.function <- function(x) { 0.1 + (1.5* abs(x)); }
spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1];
  colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2];
  #colours[sign(x) == 1] <- default.colours(1, palette.type = "spiral.sunrise");
  return(colours);
}
spot.single.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x > 0] <- default.colours(1, palette.type = "spiral.sunrise");
  return(colours);
}
# create dot map using BPG
plot_dotmap <- function(dataset, effectSize_col, background_col, xlabs, yaxis.lab, bgBins = seq(0, 0.1, 0.01)){
  #' This function is to create dotmap using BPG creat.dotmap function
  #' effectSize_col   vector of chr   column names of the effect size
  #' dataset          dataframe       df to plot
  #' backgroun_col    vector of chr   column names of the p values
  #' xlabs            vector of chr   column labels
  #' yaxis.lab        vector of chr   row names on dot map
  #' bgBins           vector of num   breaks for background colors

  dm <- create.dotmap(subset(dataset, select = effectSize_col),
              yaxis.cex = 1.5,
              xaxis.cex = 1.5,
              na.pch = 1,
              na.spot.size = 0,
              yaxis.lab = yaxis.lab,
              xaxis.lab = xlabs, 
              spot.size.function = spot.size.function,
              #spot.colour.function = spot.single.colour.function,
              spot.colour.function = spot.colour.function,
              key = list(
                space = "right",
                points = list(
                  cex = spot.size.function(seq(-2, 2, 0.25)),
                  col = spot.colour.function(seq(-2, 2, 0.25)),
                  pch = 19
                ),
                text = list(
                  lab = as.character(seq(-2, 2, 0.25)),
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
              bg.data = subset(dataset, select = background_col),
              # add a colourkey
              colourkey = TRUE,
              # set colour scheme for background data
              colour.scheme = c("black", "white"),
              # make bg colour scheme a discrete colour scheme, with breaks at these places
              at = bgBins)

  return(dm)
}

# reformat fragpipe output

reformat_fragpipe <- function(lstData){
  #' input is a list of dataset
  #' reformat the list of fragpipe output to simplier dataset
  tmp <- do.call('rbind', lstData)
  tmp <- as.data.table(tmp)
  tmp <- dplyr::select(tmp, select = -c("Predicted.IM", "IM", "iIM", "Predicted.iIM"))
  tmp <- tmp[order(-Precursor.Normalised)]
  tmp <- unique(tmp, by = c("File.Name", "Precursor.Id"))
  setnames(tmp, c("File.Name", "Stripped.Sequence", "Precursor.Normalised"),
           c("filename", "PeptideSequence", "Intensity"))
  return(tmp)
}

# plot venn for different libraries


plot_LibVenn <- function(peptides.lst, fills = c("#FAE5A1", "#ACE8E9", "#F7BEBE", "palevioletred")){
  #' this is color scheme made for just uEPS dEPS combEPS and sEV
  fit <- euler(peptides.lst)
  par(oma = c(5,5,5,5))
  plt <- plot(fit,
     #fills = c("#FAE5A1", "lightcoral"),
     col = 'white',
     fills = fills,
     edges = T,
     quantities = list(fontsize = 14), adjust_labels = T, lwd = 2, 
     legend = list(fontsize = 12, side = "bottom", 
                   nrow = 1, ncol = length(peptides.lst)), rotation = 1)
  return(plt)
}

lib_colors <- list(colors = c("#FAE5A1", "#ACE8E9", "#F7BEBE", "palevioletred"), 
                   breaks = c("uEPS1", "dEPS1", "sEV", "lEV"))


# Map fragpipe protein groups from library files
fetchQuery <- function(lst){
  #'lst vecotr list of ids to be fetched
  url <- "https://www.uniprot.org/uploadlists/"
  params = list(
    from = "ACC+ID",
    to = "GENENAME",
    format = "tsv",
    query = paste0(lst, collapse = " ")
  )
  r <- httr::POST(url, body = params, encode = "form")
  #cat(httr::content(r))
  return(r)
} 


create_proteinGroup <- function(peptides){
  #' this function is to map shared peptides to proteinGroups
  #' peptides    dataframe     peptide.tsv from the dda search from fragpipe
  #' 
  ProteinGroups <- subset(peptides,
                          select = c("Protein ID", "Gene", "Mapped Genes", "Mapped Proteins"))
  
  ProteinGroups <- unique(ProteinGroups)
  ProteinGroups$ProteinGroups <- ProteinGroups$`Protein ID`
  
  for(i in seq_along(ProteinGroups$`Protein ID`)){
    entry <- ProteinGroups$`Mapped Proteins`[i]
    if(grepl("[A-Z]", entry)){
      if(grepl(",", entry)){
        pg <- str_split(entry, ", ")[[1]]
        pgs <- lapply(pg, function(x) sub('^...([^\\|]+).*', '\\1', x))
        pg <- paste0(pgs, collapse = ";")
      }else{
        pg <- sub('^...([^\\|]+).*', '\\1', entry)
        pg <- paste(ProteinGroups$`Protein ID`[i], pg , sep = ";")}
      
      ProteinGroups$ProteinGroups[i] <- pg
    }
  }
  
  ProteinGroups$Genes <- ProteinGroups$Gene
  ProteinGroups[grepl("[A-Z]", `Mapped Genes`)]$Genes <- paste(ProteinGroups[grepl("[A-Z]", `Mapped Genes`)]$Gene,
                                                               ProteinGroups[grepl("[A-Z]", `Mapped Genes`)]$`Mapped Genes`, sep = ";")
  
  ProteinGroups$Genes <- gsub("\\,", ";", ProteinGroups$Genes)
 
  ProteinGroups <- merge(ProteinGroups, peptides,
                               by = c("Protein ID", "Gene", "Mapped Genes", "Mapped Proteins"))
  
  ProteinGroups <- subset(ProteinGroups, select = c("Peptide", "ProteinGroups", "Genes", "Protein ID"))
  ProteinGroups <- unique(ProteinGroups)
  return(ProteinGroups)
}

formatFragpipePeptideIntensity <- function(peptidesdf, file = NULL){
  #' This functions is to reformat the diann output to peptides files
  #' peptidesdf    datatable    frappipe output of the DIANN search
  
  df <- peptidesdf[Global.Q.Value < 0.01]
  df <- subset(df, select = c("SampleID", "Run", "PeptideSequence", "Protein.Ids", "ProteinGroups",
                              "Genes", "Precursor.Id", "Intensity"))
  df <- unique(df)
  df.Intensity <- dcast(df, ProteinGroups + Genes + PeptideSequence ~ SampleID, 
                        value.var = "Intensity", fun.aggregate = sum)
  df.Intensity <- as.data.table(df.Intensity)
  if(!is.null(file)){
    write.table(df.Intensity, file = file, quote = F, sep = "\t", row.names = F)
  }
  return(df.Intensity)
}

theme_ez <- function(){
  theme_classic() %+replace%
    theme(
      axis.text = element_text(size = 16, face = "plain"),
      axis.text.x = element_text(vjust = -0.5, face = 'plain'),
      axis.title = element_text(size = 20, face = "plain"),
      axis.title.x = element_text(vjust = -0.6, face = 'plain'),
      panel.border = element_rect(colour = "black", linewidth = 0.8, fill = "transparent"),
      legend.position = "bottom", legend.direction = "horizontal", 
      legend.title = element_text(size = 16, face ="bold"), legend.text = element_text(size = 14))
}

theme_library_counts <- function(){
  theme_ez() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.9))
}

plot_decile <- function(df, intensity_col, key){
  #' df   dataframe   long dataframe with intensity and samples
  #' intensity_col    column_name   column name of log2 intensity
  g <- df %>% group_by({{key}}) %>% 
    summarise(Observations = n(),
              MedianInt = median({{intensity_col}})) %>% 
    mutate(Percent = Observations*100/max(Observations)) %>% 
    mutate(Decile = cut(MedianInt, breaks = quantile(MedianInt, probs = seq(0, 1, length = 11), type = 5),
                        labels = 1:10, include.lowest = T))  %>% 
    ggplot(aes(x = Percent, y = MedianInt)) +
    geom_point(alpha = 0.6, aes(color = Decile)) +
    scale_color_brewer(palette = "Spectral", direction = -1) +
    scale_x_reverse() +
    theme_ez() + xlab("Observations (%)") + ylab(expression(log[2]~"Intensity"))
  return(g)
}


imputeLQuantile <- function(data, shift=1.8, size=0.2) {
  set.seed(8)
  data[data==0]<-NA
  mdf <- melt(data, na.rm=T)
  n<-round(dim(mdf)[1])
  sdev <- sd(mdf$value)
  m <- mean(mdf$value) - (sdev * shift)
  d<-rnorm(n*size, mean=m, sd=sdev*size)
  data[is.na(data)] <- d[round(runif(sum(is.na(data)),1,n*size))]
  return(data)
}

imputeLowerGuassian <- function(data, shift = 2, width = 0.3){
  #' data    wide dataframe   GeneNames as rows, SampleIDs as columns
  #' This function is for imputation across samples per gene
  set.seed(123)
  data[data==0] <- NA
  tmpdf <- melt(data) # get the values of only the Non-zero measurements
  if(sum(!is.na(tmpdf$value))==1){
    sdev <- 0.3
  }else{
    sdev <- sd(tmpdf$value, na.rm = T) * width
  }
  newMean <- mean(tmpdf$value, na.rm = T) - shift*sdev
  missing <- sum(is.na(tmpdf$value))
  idx <- which(is.na(tmpdf$value))
  lowGaussian <- rnorm(missing, mean = newMean, sd = sdev)
  #lowGaussian <- sample(lowGaussian) # shuffle vector
  data[is.na(data)] <- sample(lowGaussian, missing, replace = T)
  return(data)
}

plot_UpSet <- function(combination_list, top_axis = seq(0, 40000, 10000), top_axis_label = seq(0, 40000, 10000), top_yaxis_limits = c(0, 42000), bottom_axis = seq(0, 65000, 20000)){
  ta <- HeatmapAnnotation(
    " Peptides \nIntersection" = anno_barplot(comb_size(combination_list), height = unit(8, "cm"), 
                                              axis_param = list(gp = gpar(fontsize = 16),
                                                                at = top_axis), 
                                                                ylim = top_yaxis_limits,# ! CHANGE AXES
                                              labels = top_axis_label),
    annotation_name_side = "left", annotation_name_rot = 0, annotation_name_align = T,
    annotation_name_gp = gpar(fontsize = 24))
  la <- rowAnnotation(
    "Total peptides" = anno_barplot(set_size(combination_list),
                                    axis_param = list(direction = 'reverse', 
                                                      gp = gpar(fontsize = 16),
                                                      at = bottom_axis),
                                    width = unit(3, "cm"),
                                    border = F),
    annotation_name_side = 'bottom', annotation_name_rot = 0, 
    annotation_name_gp = gpar(fontsize = 20)
  )
  return(UpSet(combination_list, top_annotation = ta, left_annotation = la ))
}

# formatting
getPeptides_diann <- function(diann_df){
  #' function to subset peptide columns for the diann dataframe
  #' only peptide information is subsetted
  #' 
  df <- subset(diann_df, select = c("File.Name","Run", "Protein.Ids", "Genes", 
                              "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
                              "Precursor.Normalised"))
  setnames(df, c("Run","Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Normalised"),
           c("RunID","PeptideSequence", "ModifiedPeptideSequence", "Charge", "Intensity"))
  df.int <- dcast(df, PeptideSequence ~ RunID, fun.aggregate = sum, fill = NaN, 
                  value.var = "Intensity")
  df.int <- melt(df.int, id.vars = "PeptideSequence", value.name = "Intensity", variable.name = "RunID")
  df.var <- unique(df, by = c("File.Name", "RunID", "PeptideSequence", "Genes", "Protein.Ids"))
  df.int <- merge(df.int, df.var, by = c("RunID", 'PeptideSequence'))
  return(df)
}

plot_total_pg_barplot <- function(dataset, group, fill, 
                                  yaxis_limits = c(0, 5010),
                                  ylab = "No. of proteins \ndetected (DIA)"){
  #' dataset  df    dataframe of protein grouped data in long format
  #' group    col_name  column name of the group comparison
  #' fill     vector of chr    a named vector of colours for the corresponding groups
  #' yaxis_limits   vector of int   a vector of the yaxis limits
  #' ylab     chr   yaxis label
  #' xlab     chr   xaxis label
  df <- copy(dataset)
  plt <- df %>% group_by({{group}}) %>% 
    distinct(Gene.names, Uniprot.IDs) %>% 
    summarise(Total = n()) %>% 
    ggplot(aes(x = {{group}}, y = Total)) +
    geom_col(aes(fill = {{group}}), width = 0.6, colour = 'black', size = 0.3) +
    geom_text(aes(label = Total), vjust = 1, size = 5) +
    scale_y_continuous(expand = c(0,0), limits = yaxis_limits, labels = scales::comma) +
    scale_fill_manual(values = fill) +
    theme_ez() + ylab(ylab)
  return(plt)
}

#plot_total_pg_barplot(uEPS.proteinGroup, Dataset, fill = c(dataset_colours, "uEPS5" = "goldenrod1"))

plot_sample_pg_boxplot <- function(dataset, group, colours, 
                                  yaxis_limits = c(0, 4100),
                                  ylab = "No. of proteins \ndetected (DIA)"){
  #' dataset  df    dataframe of protein grouped data in long format
  #' group    col_name  column name of the group comparison
  #' fill     vector of chr    a named vector of colours for the corresponding groups
  #' yaxis_limits   vector of int   a vector of the yaxis limits
  #' ylab     chr   yaxis label
  #' xlab     chr   xaxis label
  df <- copy(dataset)
  plt <- df %>% group_by({{group}}, Run) %>% 
    distinct(Gene.names, Uniprot.IDs) %>% 
    summarise(Count = n()) %>% 
    ggplot(aes(x = {{group}}, y = Count)) +
    geom_quasirandom(aes(colour = {{group}})) +
    geom_boxplot(fill = "transparent", outlier.shape = F, width = 0.6) +
    stat_compare_means(method = "wilcox.test", paired = F,label.y.npc = 0.9, label.x = 1.1, size = 5) +
    scale_y_continuous(expand = c(0,0), limits = yaxis_limits, labels = scales::comma) +
    scale_colour_manual(values = colours) +
    theme_ez() + ylab(ylab)
  return(plt)
}
#plot_sample_pg_boxplot(uEPS.proteinGroup, Dataset, colours = c(dataset_colours, "uEPS5" = "goldenrod1"))

plot_sample_pg_lineplot <- function(dataset, group,matching_col, colours, 
                                   yaxis_limits = c(0, 4100),
                                   ylab = "No. of proteins \ndetected (DIA)"){
  #' dataset  df    dataframe of protein grouped data in long format
  #' group    col_name  column name of the group comparison
  #' colours     vector of chr    a named vector of colours for the corresponding groups
  #' yaxis_limits   vector of int   a vector of the yaxis limits
  #' ylab     chr   yaxis label
  #' xlab     chr   xaxis label
  df <- copy(dataset)
  plt <- df %>% group_by({{group}}, Run, {{matching_col}}) %>% 
    distinct(Gene.names, Uniprot.IDs) %>% 
    summarise(Count = n()) %>% 
    ggplot(aes(x = {{group}}, y = Count)) +
    geom_point(aes(colour = {{group}})) +
    geom_line(aes(group = {{matching_col}}), colour = "grey") +
    geom_boxplot(fill = "transparent", outlier.shape = NA, width = 0.6) +
    stat_compare_means(method = "wilcox.test", paired = F,label.y.npc = 0.9, label.x = 1.1, size = 5) +
    scale_y_continuous(expand = c(0,0), limits = yaxis_limits, labels = scales::comma) +
    scale_colour_manual(values = colours) +
    theme_ez() + ylab(ylab)
  return(plt)
}
# plot_sample_pg_lineplot(uEPS.proteinGroup, 
#                         Dataset,matching_col = pkey,
#                         colours = c(dataset_colours, "uEPS5" = "goldenrod1"))


plot_sample_pg_correlation <- function(dataset,dataset_summary, 
                                       group, group_category = c("inter", "intra"),
                                       colours = c("darkcyan", "grey"),
                                       xlabel = "Dataset" ){
  #' dataset  df    dataframe of protein correlation data in long format
  #' dataset_summary df   dataframe of the summary of correlation
  #' group    col_name  column name of the group comparison
  #' group_category chr   vector of the categories of the group
  #' colours     vector of chr    a named vector of colours for the corresponding groups
  #' xlabel     chr   xaxis label
  df.corr <- copy(dataset)
  sum.corr <- copy(dataset_summary)

  plt <- df.corr %>% 
    ggplot(aes(x = {{group}}, y = rho)) +
    geom_quasirandom(aes(colour = {{group}}), width = 0.3, show.legend = F) +
    geom_boxplot(fill = "transparent", outlier.shape = NA, width = 0.6) +
    geom_hline(data = sum.corr, aes(yintercept = MedianRho, colour = {{group}}), lty = 2, show.legend = F) +
    geom_text(data = sum.corr, aes(label = round(MedianRho, 2), y = MedianRho, colour = {{group}}),
              x = 2.5, vjust = -0.2, show.legend = F) +
    stat_compare_means(method = "wilcox.test", paired = F, label.y = 1, label.x = 1.2, size = 5) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
    scale_colour_manual(breaks = {{group_category}}, values = colours) +
    ylab(expression(rho)) + xlab(xlabel) +
    theme_ez() 
  return(plt)
}

#plot_sample_pg_correlation(uEPS.corr,dataset_summary = uEPS.corr.sum,group = Dataset)

groupISUP <- function(dataset, GG1, GG2, GGS){
  new_column <- case_when(dataset[[GG1]] == 3 & dataset[[GG2]] == 3 ~ "1",
                          dataset[[GGS]] == 6 ~ "1",
                          dataset[[GG1]] == 3 & dataset[[GG2]] == 4 ~ "2",
                          dataset[[GG1]] == 4 & dataset[[GG2]] == 3 ~ "3",
                          dataset[[GG1]] == 4 & dataset[[GG2]] == 4 ~ "4",
                          dataset[[GGS]] == 8 ~ "4",
                          dataset[[GGS]] > 8 ~ "5",
                          .default = NA)
  return(new_column)
}

treatment_cleanup <- function(procName){
  #' procName   vector of chr   vector of treatments
  tx <- copy(procName)
  tx <- gsub("Active Surveillance|Watchful Waiting", "Active Surveillance & Watchful Waiting", tx)
  tx <- gsub("DVP|RP", "Radical Prostatectomy", tx)
  tx <- gsub("Cryoablation(.*)", "Cryoablation", tx)
  tx <- gsub("AA|LHRH|HORM", "Hormone Therapy", tx)
  tx[grepl("Beam|Brachytherapy", tx)] <- "Radiation therapy"
  tx <- gsub("CHEMO", "Chemotherapy", tx)
  tx <- gsub("ND", "No data", tx)
  tx[is.na(tx)] <- "No data"
  tx[grepl("OTHER", tx)] <- "Others"
  return(tx)
}


isup_colours <- c("1" = "cornsilk",
                  "2" =  "yellow",
                  "3" =  "orange",
                  "4"= "maroon3",
                  "5" = "red", "No data" = "white")

psa_colours <- c("<4" = "peachpuff", "4-10" = "lightsalmon", "10-20" = "coral",
                 "20-50" = "orangered", ">50" = "orangered4", "No data" = "white")

psa_categorical_colours <- c("0 - 9.9" = "#FEE6CE", 
                             "10 - 19.9" = "#FDAE6B", 
                             ">= 20" = "#E6550D",
                             "No data" = "white")

df_colours <- c("uEPS1" = "#188B91", "uEPS5" = "#FB8072", "Overlap" = "slategray3")

race_colours <- c("White" = "#4B4B4B", 
                  "Black" = "#ED6C4C", 
                  "Asian" = "#4491F1", 
                  "South Asian" = "#F2A93B",
                  "Hispanic" = "#CBCCCC",
                  "No data" = "white")

tx_colours <- c("Active Surveillance & Watchful Waiting" = "khaki",
                "Radiation therapy" = "tan1",
                "Radical Prostatectomy" = "firebrick1",
                "Cryoablation" = "pink1",
                "Hormone Therapy" = "orchid3", 
                "Chemotherapy" = "turquoise4",
                "Others"= "darkslateblue",
                "No data" = "white")

group_colours <- c("Cancer" = "brown", "Non-cancer" = "aquamarine3")
ct_colours <- c("T0" = "#FFF7DF",
                "T1" = "#DDF68B",
                "T2" = "#6DC46E",
                "T3" = "#2F6D60", 
                "No data" = "white")

cM_colours <- c("MX" = "#E6E6E6",
                "M0" = "grey31", 
                "M1a" = "#F7BEBE",
                "M1" = "darkred",
                "No data" = "white")

cN_colours <- c("NX" = "#E6E6E6",
                "N0" = "grey31",
                "N1" = "#45B4BB",
                "N2" = "#188B91", 
                "No data" = "white")



