# All functions

library(data.table)
library(stringr)
library(dbplyr)
library(reshape2)
library(VennDiagram)
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

library(VennDiagram)
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

create_venn <- function(nKeys, BrewerPalette, KeyNames, data){
  #' nKeys: integer number of categories
  #' BrewerPalette: string choice of brewer palette
  #' KeyNames: vector of strings  category names
  #' data: list of data for venn diagram
  
  if(nKeys < 3){
    catfill <- brewer.pal(3, BrewerPalette)
    catfill <- sample(catfill, 2)
  }else{
    catfill <- brewer.pal(nKeys, BrewerPalette)
  }
  v1 <- venn.diagram(
    x = data,
    scaled=T,
    area.vector = T,
    category.names = KeyNames,
    fill=catfill,
    col = catfill,
    filename = NULL,
    output = FALSE ,
    imagetype="png" ,
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
    cat.col = catfill
  )
  grid.newpage()
  grid.draw(v1)
  return(v1)
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
=======
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


