# Load background libraries for quality assessments
library(data.table)

############################ Paths ######################
path <- 'D:/projects/pca_urine_spectral_lib/data/irt'

lib.files <- list.files(path, pattern = "", full.names = T)
lib.files

#################### Load files #######################

# Biognosys iRTs
irt <- fread(lib.files[9])

# CiRTs
cirt <- fread(lib.files[3])

# MQ nonlinear iRT file
uirt <- fread("14")

################## Format files ######################

irt
irt <- irt[!duplicated(irt, by = c("transition_group_id", "ModifiedPeptideSequence"))]
irt

