#!/bin/python

import os
import numpy as np
import pandas as pd


# assign directories
dir = "/home/annieha/pca_urine_spectral_lib/data/library/library_bootstrap"
all_dir = os.listdir(os.path.join(dir, "bootstrap_splib"))

# filter for directoris to work with: sample 2

sample2 = list(filter(lambda i: "sample180_" in i, all_dir))
# get full path

for i in range(0, len(sample2)):
    sample2[i] = os.path.join(dir, "bootstrap_splib", sample2[i])

# read peptide.tsv file is sufficient (library.tsv very big)
# only need the number of unique peptides
# do not waste memory storing df

lib_peptideLength=[]
lib_names=[]
lib_protein=[]

for j in range(0, len(sample2)):
    df = pd.read_csv(os.path.join(sample2[j], "peptide.tsv"), sep = "\t")
    proteins = pd.read_csv(os.path.join(sample2[j], "protein.tsv"), sep = "\t")
    lib_names.append(os.path.basename(sample2[j]))
    lib_peptideLength.append(len(df.Peptide.unique()))
    lib_protein.append(len(proteins.Protein.unique()))
    print("Bootstrap: ", lib_names[j])
    print("Number of peptides: ", lib_peptideLength[j])
    print("Number of proteins: ", lib_protein[j])

lib_numbers = {'Bootstrap': lib_names, 'Number of peptides': lib_peptideLength, 'Number of proteins': lib_protein}

lib_df = pd.DataFrame(lib_numbers)

lib_df.to_csv(os.path.join(dir, "analysis", "library_numbers_sample180.tsv"), sep = "\t", index = False)

