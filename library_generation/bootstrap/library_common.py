#!/bin/python

import os
import numpy as np
import pandas as pd
from itertools import chain

# assign directories
dir = "/home/annieha/pca_urine_spectral_lib/data/library/library_bootstrap"
all_dir = os.listdir(os.path.join(dir, "bootstrap_splib"))


# start loop with each increments

increments = [i for i in range(2, 20, 2)]
for i in range(20, 199, 20):
    increments.append(i)
increments.append(199)

all_peptide = []
all_protein = []
intersection_peptide = []
intersection_protein = []

for j in range(0, len(increments)):
    sampleSize = str(increments[j])
    sampleSize = "sample"+sampleSize+"_"
    # filter per increment sample size
    samplelst = list(filter(lambda i : sampleSize in i, all_dir))
    # get full path
    for k in range(0, len(samplelst)):
        samplelst[k] = os.path.join(dir, "bootstrap_splib", samplelst[k])
    lib_peptide = []
    lib_protein = []
    for x in range(0, len(samplelst)):
        df = pd.read_csv(os.path.join(samplelst[x], "peptide.tsv"), sep = "\t")
        df_protein = pd.read_csv(os.path.join(samplelst[x], "protein.tsv"), sep = "\t")
        lib_peptide.append(df.Peptide.unique())
        lib_protein.append(df_protein.Protein.unique())
    total_peptide = len(set(chain(*lib_peptide)))
    common_peptide = len(set.intersection(*map(set, lib_peptide)))
    total_protein = len(set(chain(*lib_protein)))
    common_protein = len(set.intersection(*map(set, lib_protein)))
    all_peptide.append(total_peptide)
    intersection_peptide.append(common_peptide)
    all_protein.append(total_protein)
    intersection_protein.append(common_protein)


lib_numbers = {'Sample size': increments, 'Number of total peptides': all_peptide, 'Number of common peptides': intersection_peptide, "Number of total proteins": all_protein, "Number of common proteins": intersection_protein}

lib_df = pd.DataFrame(lib_numbers)

lib_df.to_csv(os.path.join(dir, "analysis", "library_common_entries.tsv"), sep = "\t", index = False)


