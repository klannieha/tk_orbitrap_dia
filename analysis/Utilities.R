# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
#library(VennDiagram)
library(Matrix)
library(httr)

#################### Functions for Use ###################

binColumn <- function(data, column_name, binwidth){
  # Purpose:
  #       Bin the data to shrink the data size
  # parameters:
  #       data          data.table  full maxquant output library
  #       column_name   str         columns to be binned
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
}

SwathWindows <- function(data, start, end, WinNum){
  mq <- data[ `m/z` > start &  `m/z` < end]
  mq <- mq[order( `m/z`)]
  mq <- mq[!duplicated(mq, by=c("ModifiedPeptideSequence", "Charge"))]
  mq <- binColumn(mq, "m/z", 1)
  average_ions <- nrow(mq)/WinNum
  swathdf <- data.frame(Window = seq(1, WinNum), Start = numeric(length = WinNum),
                        End = numeric(length = WinNum), width = numeric(length = WinNum), stringsAsFactors = F)
  id_start <- 1
  id_end <- round(average_ions)
  for (i in seq(1, WinNum)) {
    prec <- mq[id_start: id_end]
    swathdf$Start[i] <- prec$`m/z`[1]
    swathdf$End[i] <- prec$`m/z`[nrow(prec)]
    id_start <- id_start + nrow(prec)
    id_end <- id_end + nrow(prec)
  }
  swathdf$End[WinNum] <- end
  swathdf$Start <- round(swathdf$Start)
  swathdf$End <- round(swathdf$End)
  swathdf$width <- swathdf$End - swathdf$Start
  return(swathdf)
}



# function for uniform modification tags
reformat_mods <- function(col){
  modifiedSequence <- col
  modifiedSequence <- gsub("_", "", modifiedSequence)
  modifiedSequence <- gsub("(ox)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("(ph)","Phospho", modifiedSequence)
  modifiedSequence <- gsub("C","C(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("(ac)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:4\\)","(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:1\\)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:35\\)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:21\\)","Phospho", modifiedSequence)
  
  column <- modifiedSequence
  return(column)
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
    if(length(p) > 1){
      p <- paste0(p, collapse=";")
    }
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
    r <- r[!duplicated(r, by=c("Sequence", "BioReplicate"))]
    r <- dcast(r, Sequence~BioReplicate, value.var = "Intensity")
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

getCV <- function(x){
  CV = apply(x, 1, function(y) 100*sd(y, na.rm = T)/mean(y, na.rm = T))
  mean = apply(x, 1, function(y) mean(y, na.rm = T))
  x$CV <- CV
  x$mean <- mean
  return(x)
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
    column.names <- c('ProteinName', 'FullPeptideName', 'Sequence', 'Charge', 'aggr_Fragment_Annotation', 
                      'aggr_Peak_Area', 'filename', 'm_score', 'decoy', "Intensity", "RT", "run_id", 
                      "transition_group_id", "d_score", "peak_group_rank")
  }
  if(length(column.names) > length(column.names[column.names %in% colnames(data)])){
    col.names.missing <- column.names[!column.names %in% colnames(data)]
    warning("These columns are missing from the data:", paste(unlist(col.names.missing), collapse=", "))
    
  }
  # Keep only required columns for MSStats and mapDIA
  if(length(column.names) == length(column.names[column.names %in% colnames(data)])){
    data.filtered <- data[,column.names]
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
    m <- paste0(n[2:4], collapse = " ")
    r <- n[length(n)]
    c(m, r)
  })
  method <- do.call('rbind', method)
  df$Condition <- method[,1]
  df$Replicate <- method[,2]
  return(df)
}

annotate_MStern <- function(pyprophet_out){
  df <- data.frame(filename = unique(pyprophet_out$filename),
                   stringsAsFactors = F)
  df$Run <- str_extract(df$filename, "R[0-9]+")
  df$Condition <- str_extract(df$filename, "S[0-9]+")
  df$BioReplicate <- str_sub(df$filename, start = -10, end = -9)
  return(df)
}

#library(VennDiagram)
library(RColorBrewer)

create_venn <- function(data, category_names, category_colours, file = NULL){
  if((length(data) != length(category_names)) | (length(data) != length(category_colours))){
    print("Arguments are not of the same length")
  }else if((length(data) == length(category_names)) & length(data) == length(category_colours)){
    v1 <- venn.diagram(
      x = data,
      scaled=T,
      area.vector = T,
      category.names = category_names,
      fill=category_colours,
      col = category_colours,
      filename = file,
      height = 480 , 
      width = 480 , 
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      cex = 2,
      fontfamily = "sans",
      cat.cex = 2,
      cat.default.pos = "outer",
      cat.fontfamily = "sans",
      cat.col = category_colours
    )
  }
  grid.newpage()
  grid.draw(v1)
}

fetchQuery <- function(lst){
  #'lst vecotr list of ids to be fetched
  url <- "https://www.uniprot.org/uploadlists/"
  
  params = list(
    from = "ACC+ID",
    to = "GENENAME",
    format = "tab",
    query = paste0(lst, collapse = " ")
  )
  r <- httr::POST(url, body = params, encode = "form")
  #cat(httr::content(r))
  return(r)
} 

map2Gene <- function(protein_id){
  #' protein_id vector list of protein_ids to be mapped to Gene names

  cleanp <- protein_id[!grepl(";", protein_id)] # IDs that are unique, not groups
  fetch_cleanp <- fetchQuery(cleanp)
  strsp <- strsplit(content(fetch_cleanp), "\n")[[1]]
  strsp2 <- strsplit(strsp, "\t")
  protein_df <- lapply(2:length(strsp), function(x){
    from <- strsp2[[x]][1]
    to <- strsp2[[x]][2]
    c(from, to)
  })
  
  protein_df <- do.call('rbind', protein_df) %>% as.data.table() %>% as.data.table()
  colnames(protein_df) <- c("ProteinName", "GeneName")
  protein_df <- protein_df[!duplicated(ProteinName)]

  # got df with unique protein names and genes
  # this takes too long
  # change it into canonical and protein groups
  groupP <- protein_id[grepl(";", protein_id)] # get protein group IDs
  # test: 2503 IDs
  groupP.df <- data.frame(Group = groupP, Canonical = character(length = length(groupP)), stringsAsFactors = F)
  groupP.df$Canonical <- unlist(lapply(groupP, GeneFormat))
  
  groupP_noIso <- groupP.df$Canonical[!grepl(";", groupP.df$Canonical)] %>% unique()
  # test: 1358
  fetch_noIso <- fetchQuery(groupP_noIso)
  s <- strsplit(content(fetch_noIso), "\n")[[1]]
  s2 <- strsplit(s, "\t")
  s2 <- do.call('rbind', s2[-1])
  groupP_noIso <- data.frame(ProteinName = s2[,1], GeneName = s2[, 2],
                             stringsAsFactors = F)
  groupP_noIso <- as.data.table(groupP_noIso)
  groupP_noIso <- groupP_noIso[!duplicated(ProteinName)]
  #groupP.df <- merge(groupP.df, groupP_noIso, by.x = "Canonical", by.y = "ProteinName")
  
  groupP <- groupP.df$Canonical[grepl(";", groupP.df$Canonical)] %>% unique()
  groupP.str <- unique(unlist(strsplit(groupP, ";")))
  fetch_pg <- fetchQuery(groupP.str)
  s <- strsplit(content(fetch_pg), "\n")[[1]]
  s2 <- strsplit(s, "\t")
  s2 <- do.call('rbind', s2[-1]) %>% as.data.frame()
  
  groupP_Iso <- data.frame(ProteinName = groupP, stringsAsFactors = F)
  GeneNames <- lapply(1:nrow(groupP_Iso), function(x){
    id <- which(s2[, 1] %in% strsplit(groupP_Iso$ProteinName[x], ";")[[1]])
    g <- paste0(s2[id,2], collapse = ";")
    g
  })
  GeneNames <- unlist(GeneNames)
  groupP_Iso$GeneName <- GeneNames
  groupP_Iso <- distinct_all(groupP_Iso)
  groupP_allIso <- rbind(groupP_noIso, groupP_Iso)
  groupP.df <- merge(groupP.df, groupP_allIso, by.x = "Canonical", by.y = "ProteinName", all.x = T)
  groupP.df$Canonical <- NULL
  setnames(groupP.df, "Group", "ProteinName")

  proteinGroups <- rbind(groupP.df, protein_df)
  return(proteinGroups)
}

# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
#library(VennDiagram)
#library(Matrix)

#################### Functions for Use ###################

binColumn <- function(data, column_name, binwidth){
  # Purpose:
  #       Bin the data to shrink the data size
  # parameters:
  #       data          data.table  full maxquant output library
  #       column_name   str         columns to be binned
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
}

SwathWindows <- function(data, start, end, WinNum){
  mq <- data[ `m/z` > start &  `m/z` < end]
  mq <- mq[order( `m/z`)]
  mq <- mq[!duplicated(mq, by=c("ModifiedPeptideSequence", "Charge"))]
  mq <- binColumn(mq, "m/z", 1)
  average_ions <- nrow(mq)/WinNum
  swathdf <- data.frame(Window = seq(1, WinNum), Start = numeric(length = WinNum),
                        End = numeric(length = WinNum), width = numeric(length = WinNum), stringsAsFactors = F)
  id_start <- 1
  id_end <- round(average_ions)
  for (i in seq(1, WinNum)) {
    prec <- mq[id_start: id_end]
    swathdf$Start[i] <- prec$`m/z`[1]
    swathdf$End[i] <- prec$`m/z`[nrow(prec)]
    id_start <- id_start + nrow(prec)
    id_end <- id_end + nrow(prec)
  }
  swathdf$End[WinNum] <- end
  swathdf$Start <- round(swathdf$Start)
  swathdf$End <- round(swathdf$End)
  swathdf$width <- swathdf$End - swathdf$Start
  return(swathdf)
}



# function for uniform modification tags
reformat_mods <- function(col){
  modifiedSequence <- col
  modifiedSequence <- gsub("_", "", modifiedSequence)
  modifiedSequence <- gsub("(ox)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("(ph)","Phospho", modifiedSequence)
  modifiedSequence <- gsub("C","C(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("(ac)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:4\\)","(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:1\\)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:35\\)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:21\\)","Phospho", modifiedSequence)
  
  column <- modifiedSequence
  return(column)
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
    data.filtered <- data[,..column.names]
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
    m <- paste0(n[2:4], collapse = " ")
    r <- n[length(n)]
    c(m, r)
  })
  method <- do.call('rbind', method)
  df$Condition <- method[,1]
  df$Replicate <- method[,2]
  return(df)
}


# All functions

library(data.table)
library(tidyverse)
library(stringr)
library(dbplyr)
library(reshape2)
library(VennDiagram)
library(Matrix)
library(httr)
library(BoutrosLab.plotting.general)

#################### Functions for Use ###################

binColumn <- function(data, column_name, binwidth){
  # Purpose:
  #       Bin the data to shrink the data size
  # parameters:
  #       data          data.table  full maxquant output library
  #       column_name   str         columns to be binned
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
}

SwathWindows <- function(data, start, end, WinNum){
  mq <- data[ `m/z` > start &  `m/z` < end]
  mq <- mq[order( `m/z`)]
  mq <- mq[!duplicated(mq, by=c("ModifiedPeptideSequence", "Charge"))]
  mq <- binColumn(mq, "m/z", 1)
  average_ions <- nrow(mq)/WinNum
  swathdf <- data.frame(Window = seq(1, WinNum), Start = numeric(length = WinNum),
                        End = numeric(length = WinNum), width = numeric(length = WinNum), stringsAsFactors = F)
  id_start <- 1
  id_end <- round(average_ions)
  for (i in seq(1, WinNum)) {
    prec <- mq[id_start: id_end]
    swathdf$Start[i] <- prec$`m/z`[1]
    swathdf$End[i] <- prec$`m/z`[nrow(prec)]
    id_start <- id_start + nrow(prec)
    id_end <- id_end + nrow(prec)
  }
  swathdf$End[WinNum] <- end
  swathdf$Start <- round(swathdf$Start)
  swathdf$End <- round(swathdf$End)
  swathdf$width <- swathdf$End - swathdf$Start
  return(swathdf)
}



# function for uniform modification tags
reformat_mods <- function(col){
  modifiedSequence <- col
  modifiedSequence <- gsub("_", "", modifiedSequence)
  modifiedSequence <- gsub("(ox)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("(ph)","Phospho", modifiedSequence)
  modifiedSequence <- gsub("C","C(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("(ac)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:4\\)","(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:1\\)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:35\\)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:21\\)","Phospho", modifiedSequence)
  
  column <- modifiedSequence
  return(column)
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
    if(length(p) > 1){
      p <- paste0(p, collapse=";")
    }
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
    r <- r[!duplicated(r, by=c("Sequence", "BioReplicate"))]
    r <- dcast(r, Sequence~BioReplicate, value.var = "Intensity")
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

getCV <- function(x){
  CV = apply(x, 1, function(y) 100*sd(y, na.rm = T)/mean(y, na.rm = T))
  mean = apply(x, 1, function(y) mean(y, na.rm = T))
  x$CV <- CV
  x$mean <- mean
  return(x)
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
    column.names <- c('ProteinName', 'FullPeptideName', 'Sequence', 'Charge', 'aggr_Fragment_Annotation', 
                      'aggr_Peak_Area', 'filename', 'm_score', 'decoy', "Intensity", "RT", "run_id", 
                      "transition_group_id", "d_score", "peak_group_rank")
  }
  if(length(column.names) > length(column.names[column.names %in% colnames(data)])){
    col.names.missing <- column.names[!column.names %in% colnames(data)]
    warning("These columns are missing from the data:", paste(unlist(col.names.missing), collapse=", "))
    
  }
  # Keep only required columns for MSStats and mapDIA
  if(length(column.names) == length(column.names[column.names %in% colnames(data)])){
    data.filtered <- data[,column.names]
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
    m <- paste0(n[2:4], collapse = " ")
    r <- n[length(n)]
    c(m, r)
  })
  method <- do.call('rbind', method)
  df$Condition <- method[,1]
  df$Replicate <- method[,2]
  return(df)
}

annotate_MStern <- function(pyprophet_out){
  df <- data.frame(filename = unique(pyprophet_out$filename),
                   stringsAsFactors = F)
  df$Run <- str_extract(df$filename, "R[0-9]+")
  df$Condition <- str_extract(df$filename, "S[0-9]+")
  df$BioReplicate <- str_sub(df$filename, start = -10, end = -9)
  return(df)
}

library(VennDiagram)
library(RColorBrewer)

create_venn <- function(data, category_names, category_colours, file = NULL){
  if((length(data) != length(category_names)) | (length(data) != length(category_colours))){
    print("Arguments are not of the same length")
  }else if((length(data) == length(category_names)) & length(data) == length(category_colours)){
    v1 <- venn.diagram(
      x = data,
      scaled=F,
      euler.d = T,
      rotation = 1, 
      category.names = category_names,
      fill=category_colours,
      col = category_colours,
      filename = file,
      height = 480 , 
      width = 480 , 
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      cex = 2,
      fontfamily = "sans",
      #cat.dist = c(0.05, 0.05, 0.02),
      #cat.pos = c(-20, 0, 20),
      cat.cex = 2,
      #ext.pos = 2,
      cat.default.pos = "outer",
      cat.fontfamily = "sans",
      cat.col = category_colours
    )
  }
  grid.newpage()
  grid.draw(v1)
}


isJobReady <- function(jobId) {
  pollingInterval = 5
  nTries = 20
  for (i in 1:nTries) {
    url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return (FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}

getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
    url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    url <- gsub("/results/", "/results/stream/", redirectURL)
  }
}



files = list(
  ids = ids,
  from = "UniProtKB_AC-ID",
  to = "Gene_Name"
)
r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
submission <- content(r, as = "parsed")

if (isJobReady(submission[["jobId"]])) {
  url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
  r <- GET(url = url, accept_json())
  details <- content(r, as = "parsed")
  url <- getResultsURL(details[["redirectURL"]])
  # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
  url <- paste(url, "?format=tsv", sep = "")
  r <- GET(url = url, accept_json())
  resultsTable = read.table(text = content(r), sep = "\t", header=TRUE)
  print(resultsTable)
}



fetchQuery <- function(lst){
  #'lst vecotr list of ids to be fetched
  
  params = list(
    from = "UniProtKB_AC-ID",
    to = "Gene_Name",
    ids = paste0(lst, collapse = ", ")
  )
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = params, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")
  
  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = content(r), sep = "\t", header=TRUE)
    #print(resultsTable)
  } else{
    
  }
  
  return(resultsTable)
} 


map2Gene <- function(protein_id){
  #' protein_id vector list of protein_ids to be mapped to Gene names
  
  cleanp <- protein_id[!grepl(";", protein_id)] # IDs that are unique, not groups
  fetch_cleanp <- fetchQuery(cleanp) # results is a table already
  
  protein_df <- fetch_cleanp %>% as.data.table()
  setnames(protein_df, c("From", "To"), c("ProteinName", "GeneName"))

  
  protein_df <- protein_df[!duplicated(ProteinName)]
  
  # got df with unique protein names and genes
  # this takes too long
  # change it into canonical and protein groups
  groupP <- protein_id[grepl(";", protein_id)] # get protein group IDs
  # test: 2503 IDs
  groupP.df <- data.frame(Group = groupP, Canonical = character(length = length(groupP)), stringsAsFactors = F)
  groupP.df$Canonical <- unlist(lapply(groupP, GeneFormat))
  
  groupP_noIso <- unique(groupP.df$Canonical[!grepl(";", groupP.df$Canonical)])
  # test: 1358
  fetch_noIso <- fetchQuery(groupP_noIso) %>%  as.data.table()# again, shows tables as a result
  
  setnames(fetch_noIso, c("From", "To"), c("ProteinName", "GeneName"))
  
  fetch_noIso <- fetch_noIso[!duplicated(ProteinName)]
  
  groupP <- groupP.df$Canonical[grepl(";", groupP.df$Canonical)] %>% unique()
  groupP.str <- unique(unlist(strsplit(groupP, ";")))
  fetch_pg <- fetchQuery(groupP.str) %>% as.data.table()
  
  fetch_pg <- setnames(fetch_pg, c("From", "To"), c("ProteinName", "GeneName"))
  
  groupP_Iso <- data.frame(ProteinName = groupP, stringsAsFactors = F)
  GeneNames <- lapply(1:nrow(groupP_Iso), function(x){
    id <- which(fetch_pg$ProteinName %in% strsplit(groupP_Iso$ProteinName[x], ";")[[1]])
    g <- paste0(fetch_pg$GeneName[id], collapse = ";")
    g
  })
  GeneNames <- unlist(GeneNames)
  groupP_Iso$GeneName <- GeneNames
  groupP_Iso <- distinct_all(groupP_Iso)
  groupP_allIso <- rbind(groupP_noIso, groupP_Iso, fetch_noIso)
  groupP.df <- merge(groupP.df, groupP_allIso, by.x = "Canonical", by.y = "ProteinName", all.x = T)
  groupP.df$Canonical <- NULL
  setnames(groupP.df, "Group", "ProteinName")
  
  proteinGroups <- rbind(groupP.df, protein_df)
  return(proteinGroups)
}
#=======
# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
#library(VennDiagram)
#library(Matrix)

#################### Functions for Use ###################

binColumn <- function(data, column_name, binwidth){
  # Purpose:
  #       Bin the data to shrink the data size
  # parameters:
  #       data          data.table  full maxquant output library
  #       column_name   str         columns to be binned
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
}

SwathWindows <- function(data, start, end, WinNum){
  mq <- data[ `m/z` > start &  `m/z` < end]
  mq <- mq[order( `m/z`)]
  mq <- mq[!duplicated(mq, by=c("ModifiedPeptideSequence", "Charge"))]
  mq <- binColumn(mq, "m/z", 1)
  average_ions <- nrow(mq)/WinNum
  swathdf <- data.frame(Window = seq(1, WinNum), Start = numeric(length = WinNum),
                        End = numeric(length = WinNum), width = numeric(length = WinNum), stringsAsFactors = F)
  id_start <- 1
  id_end <- round(average_ions)
  for (i in seq(1, WinNum)) {
    prec <- mq[id_start: id_end]
    swathdf$Start[i] <- prec$`m/z`[1]
    swathdf$End[i] <- prec$`m/z`[nrow(prec)]
    id_start <- id_start + nrow(prec)
    id_end <- id_end + nrow(prec)
  }
  swathdf$End[WinNum] <- end
  swathdf$Start <- round(swathdf$Start)
  swathdf$End <- round(swathdf$End)
  swathdf$width <- swathdf$End - swathdf$Start
  return(swathdf)
}



# function for uniform modification tags
reformat_mods <- function(col){
  modifiedSequence <- col
  modifiedSequence <- gsub("_", "", modifiedSequence)
  modifiedSequence <- gsub("(ox)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("(ph)","Phospho", modifiedSequence)
  modifiedSequence <- gsub("C","C(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("(ac)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:4\\)","(Carbamidomethyl)", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:1\\)","Acetylation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:35\\)","Oxidation", modifiedSequence)
  modifiedSequence <- gsub("\\(UniMod:21\\)","Phospho", modifiedSequence)
  
  column <- modifiedSequence
  return(column)
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

uEPS.names <- fread("D:/projects/pca_urine_spectral_lib/data/pca_dda/pca_postDRE_urine/uEPS1_names.txt")
uEPS.names
uEPS.names$Sample.Old <- str_extract(uEPS.names$Old_name, "S[0-9]+")
uEPS.names$Sample.New <- str_extract(uEPS.names$New_name, "UP[0-9]+")
runs_annotation <- fread("D:/projects/pca_urine_spectral_lib/data/openswath/20210201_MStern/annotation.tsv")
runs_annotation %>% View()
runs_annotation <- merge(uEPS.names, runs_annotation, by.x = "Sample.Old",
                         by.y = "Condition")

runs_annotation$Batch <- gsub("[^0-9]", "", runs_annotation$Batch)


#------ MQ data extraction ----------------------
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

iRTs <- fread("D:/projects/pca_urine_spectral_lib/data/irt/irt_transitions.tsv")
iRTs
iRTs <- iRTs %>%
  filter(rank == 1) %>%
  select(`Q1 monoisotopic`, `precursor charge`, `nominal sequence`, `sequence id` ) %>%
  distinct_all()

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
Htest(p, "1", ">1", ClinicISUP, pg.de.impute, "utest", logTransform = T)

# calculate log2FC

calculatePeptidelogFC <- function(data, colname, group1, group2){
  #' data  long dataframe with Log2Intensity
  #' variable column name
  #' group1 group 1 from the column
  #' group2 group2 from the column
  FC <- data %>%
    group_by({{colname}}, Sequence) %>%
    summarise(Mean = mean(Log2Intensity, na.rm = T), 
              Observation = n()) %>% 
    pivot_wider(id_cols = "Sequence", names_from = {{colname}}, values_from = c("Mean", "Observation")) %>%# head()
    mutate(FC = {{group2}} - {{group1}}) %>% 
    as.data.table()
  return(FC)
}

calculatePeptidelogFC(peptides.top3.lnorm, Grade, `1`, `2+`)

calculateProteinlogFC <- function(data, colname, group1, group2){
  #' data  long dataframe with Log2Intensity
  #' variable column name
  #' group1 group 1 from the column
  #' group2 group2 from the column
  FC <- data %>%
    group_by({{colname}}, ProteinName) %>%
    summarise(Mean = mean(Log2Intensity, na.rm = T)) %>% 
    pivot_wider(id_cols = "ProteinName", names_from = {{colname}}, values_from = "Mean") %>%# head()
    mutate(FC = {{group2}} - {{group1}}) %>% 
    as.data.table()
  return(FC)
}


plot_volcano <- function(data, label = F, pvalue_column,padj_column, FC_column, sig_level = 0.05, legends = F){
  #' data    dataframe with FC and p values
  #' label   boolean   include labels or not
  #' add other legends if needed
  df <- copy(data)
  
  df <- df %>% mutate(Sig = ifelse({{padj_column}} < sig_level, TRUE, FALSE)) %>%
    mutate(Enriched = case_when((Sig & {{FC_column}} > 0) ~ "UP",
                                (Sig & {{FC_column}} < 0) ~ "DOWN"))
  df <- df %>% mutate(LABEL = ifelse(Sig, paste(peptide_sequence, GeneName, sep = " * "), NA))
  p <- df %>%
      mutate(`-log10(p-value)` = -log10({{pvalue_column}})) %>%
      ggplot(aes(x = {{FC_column}}, y = `-log10(p-value)`)) +
      geom_point(aes(color = Enriched), show.legend = legends) +
      scale_color_manual(breaks = c("UP", "DOWN"), values = c("red", "blue")) +
      plot_theme()

  if(label){
    p <- p + geom_label_repel(aes(label = LABEL), max.overlaps = 30)
  }
  p
  return(p)
}
 

spot.size.function <- function(x) { 0.1 + (2 * abs(x)); }
spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1];
  colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2];
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


plot_LibVenn <- function(peptides.lst){
  #' this is color scheme made for just uEPS dEPS combEPS and sEV
  fit <- euler(peptides.lst)
  plt <- plot(fit,
     #fills = c("#FAE5A1", "lightcoral"),
     col = 'white',
     fills = c("#FAE5A1", "#ACE8E9", "#F7BEBE", "palevioletred"),
     edges = T,
     quantities = list(fontsize = 18), adjust_labels = T, lwd = 2, 
     legend = list(fontsize = 24, side = "bottom", nrow = 1, ncol = 4), rotation = 1)
  return(plt)

}

lib_colors <- list(colors = c("#FAE5A1", "#ACE8E9", "lightgrey", "#F7BEBE"), breaks = c("uEPS1", "dEPS1", "combEPS", "sEV"))


# Map fragpipe protein groups from library files

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

formatFragpipePeptideIntensity <- function(peptidesdf){
  #' This functions is to reformat the diann output to peptides files
  #' peptidesdf    datatable    frappipe output of the DIANN search
  
  df <- peptidesdf[Global.Q.Value < 0.01]
  df <- subset(df, select = c("SampleID", "Run", "PeptideSequence", "Protein.Ids", "ProteinGroups",
                              "Genes", "Precursor.Id", "Intensity"))
  df <- unique(df)
  df.Intensity <- dcast(df, ProteinGroups + Genes + PeptideSequence ~ SampleID, 
                        value.var = "Intensity", fun.aggregate = sum)
  df.Intensity <- as.data.table(df.Intensity)
  return(df.Intensity)
}


