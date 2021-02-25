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


