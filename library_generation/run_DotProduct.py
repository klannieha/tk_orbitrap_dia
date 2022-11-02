#!/usr/env/ python

############## This script is to run the similarity score ##############
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re
from spectral_angle import *
import seaborn as sns
import multiprocessing as mp
from multiprocessing import Manager
from itertools import repeat
import ctypes

################ Performance #####################
def spitTime(func, params):
    wrapped = wrapper(func, params)
    time = timeit.timeit(wrapped, number = 1)
    print(time)
    return(time)


#def calculateDotPro(precursor_id, bin_width=0.01, data):
#    print("Analyzing dot product for precursor # ", precursor_id)
#    specs = data[data['precursor_id'] == precursor_id]
#    print("Creating bins")
#    bins = Createbins(specs, precursor_id, bin_width)
#    print("Subsetting spectra by cohort")
#    binned_spec1 = specs[specs['cohort'] == "uEPS"]
#    binned_spec2 = specs[specs['cohort'] == "dEPS"]
#    print("Binning spec 1")
#    binned_spec1 = BinSpectrum(binned_spec1, bins.bins)
#    binned_spec1 = mp.Array('d', binned_spec1)
#    print("Binning spec2")
#    binned_spec2 = BinSpectrum(binned_spec2, bins.bins)
#    binned_spec2 = mp.Array('d', binned_spec2)
#    print("calculating cos simularity")
#    bin_weights = mp.Array('d', bins.weight)
#    DotProduct = NormDotP(binned_spec1, binned_spec2, bins.weight)
#    print("DotProduct of peptide #  ", precursor_id, "is ", DotProduct)
#    out_dp = pd.DataFrame({"peptide" : [specs['FullUniModPeptideName'].values[1]], "charge": [specs["PrecursorCharge"].values[1]], 'Score':[DotProduct[0]]}, columns = ['peptide', 'charge', 'Score'])
#    if precursor_id == 0:
#        print("Started calculating dot products...")
#    return(out_dp)


################ load data #####################

base = "/project/6002011/annieha/pca_urine_spectral_lib/data/library/FragpipeLibrary"
libraries = os.listdir(base)
libraries = [x for x in libraries if re.search(".tsv", x)]

dEPS = pd.read_csv(os.path.join(base, libraries[0]), sep = "\t")
sEV = pd.read_csv(os.path.join(base, libraries[1]), sep = "\t")
lEV = pd.read_csv(os.path.join(base, libraries[2]), sep = "\t")

uEPS = pd.read_csv(os.path.join(base,libraries[3]), sep = "\t")

# intersection peptides
#peptxt = os.path.join(base, 'peptideIntersect.txt')
#peptides = peptxt.read().splitlines()
precursors = pd.read_csv(os.path.join(base, "precursorIntersect.txt"), sep = "\t")

precursors = precursors.assign(precursor_id = range(0, precursors.shape[0]))

#uEPS_unique = uEPS.loc[:, ['FullUniModPeptideName', 'PrecursorCharge']]
#uEPS_unique = uEPS_unique.drop_duplicates()
#print("PostDRE urine unique precursors: ", uEPS_unique.shape[0])

#dEPS_unique = dEPS.loc[:, ['FullUniModPeptideName', 'PrecursorCharge']]
#dEPS_unique = dEPS_unique.drop_duplicates()
#print("DirectEPS unique precursors: ", dEPS_unique.shape[0])

#EPS_overlap = pd.merge(uEPS_unique, dEPS_unique, on = ["FullUniModPeptideName", "PrecursorCharge"])
#EPS_overlap = EPS_overlap.assign(precursor_id = range(0, EPS_overlap.shape[0]))
#print("EPS overlapping precursors: ", EPS_overlap.shape[0]) #43564

# Filter for the libraries with the overlapping precursors

uEPS_specs = pd.merge(uEPS, precursors, on = ["ModifiedPeptideSequence", "PrecursorCharge"])
uEPS_specs = uEPS_specs.assign(cohort = "uEPS")
dEPS_specs = pd.merge(dEPS, precursors, on = ["ModifiedPeptideSequence", "PrecursorCharge"])
dEPS_specs = dEPS_specs.assign(cohort = "dEPS")
sEV_specs = pd.merge(sEV, precursors, on = ["ModifiedPeptideSequence", "PrecursorCharge"])
sEV_specs = sEV_specs.assign(cohort = "sEV")

lEV_specs = pd.merge(lEV, precursors, on = ["ModifiedPeptideSequence", "PrecursorCharge"])
lEV_specs = lEV_specs.assign(cohort = "lEV")

#EPS_specs = uEPS_specs.append(dEPS_specs)
#EPS_specs = EPS_specs.sort_values(by = "precursor_id")

all_specs = pd.concat([uEPS_specs, dEPS_specs, sEV_specs, lEV_specs], ignore_index = True)
all_specs = all_specs.drop(columns=['Library'])
#all_specs = all_specs.rename(columns={'cohort':'Library'})

all_specs = all_specs.sort_values(by = 'precursor_id')

##################### test ###########################

s = all_specs[all_specs['precursor_id'] == 107]

bins = Createbins(s, 107, 0.005)
binned_specu = s[s['cohort'] == "uEPS"]
binned_specd = s[s['cohort'] == "dEPS"]
binned_specu = BinSpectrum(binned_specu, bins.bins)
binned_specd = BinSpectrum(binned_specd, bins.bins)

#dp = NormDotP(binned_specu, binned_specd, bins.weight, 107)


################ Combine pipeline to one #####################
import time
# first create the labeled peptides in memory
precursorIDs = all_specs['precursor_id'].unique()


binned_specu = dict()
binned_specd = dict()
binned_specsEV = dict()
binned_speclEV = dict()
binned_weights = dict()
bins = dict()
print("Creating bins.....")
start = time.time()
for p in range(30000, len(precursorIDs)):
#    print("creating bins...")
    specs = all_specs[all_specs['precursor_id'] == p]
    b = Createbins(specs, p, 0.005)
    b1 = specs[specs['cohort'] == 'uEPS']
    b2 = specs[specs['cohort'] == 'dEPS']
    b3 = specs[specs['cohort'] == 'sEV']
    b4 = specs[specs['cohort'] == 'lEV']
    b1 = BinSpectrum(b1, b.bins)
#    binned_specu.append(b1)
    binned_specu[p] = b1
    b2 = BinSpectrum(b2, b.bins.to_numpy())
#    binned_specd.append(b2)
    binned_specd[p] = b2
    b3 = BinSpectrum(b3, b.bins)
    binned_specsEV[p] = b3
    b4 = BinSpectrum(b4, b.bins)
    binned_speclEV[p] = b4
#    binned_weights.append(bins.weight)
    binned_weights[p] = b.weight.to_numpy()
#    bins.append(bins.bins)
    bins[p] = b.bins.to_numpy()

end = time.time()
print("Created bins in", end-start, "s.")
print("Done processing bins.")



#threads = os.cpu_count()
threads = 4
results = pd.DataFrame()
results = results.assign(precursor_id = precursorIDs[30000:len(precursorIDs)])
results = pd.merge(results, precursors.iloc[30000:len(precursorIDs)], on="precursor_id")
#args = zip(binned_specu, binned_specd, binned_weights, precursorIDs[0:10].tolist())

#args = [(binned_specu[i], binned_specd[i], binned_weights[i], precursorIDs[i]) for i in range(0, 200)] 
res = []
res_sEV = []
res_lEV = []

for i in range(30000:len(precursorIDs)):
    dp_dEPS = NormDotP(binned_specu, binned_specd, binned_weights, i)
    dp_sEV = NormDotP(binned_specu, binned_specsEV, binned_weights, i)
    dp_lEV = NormDotP(binned_specu, binned_speclEV, binned_weights, i)
    res.append(dp_dEPS[0])
    res_sEV.append(dp_sEV[0])
    res_lEV.append(dp_lEV[0])

res = np.hstack(res)
res_sEV = np.hstack(res_sEV)
res_lEV = np.hstack(res_lEV)

results = results.assign(dEPS_score = res, sEV_score = res_sEV, lEV_score = res_lEV)

#if __name__ == '__main__':
#    with mp.Pool(processes=4) as pool:
#        cal = pool.starmap(NormDotP, args)
#        dp = [res for res in cal]
        
#results = results.assign("Score", dp)

#with mp.get_context("spawn").Pool(processes=4) as pool:
#    cal = pool.map(calculateDotPro, precursors[0:4])




#rs = []


#results_df = pd.DataFrame(columns = ['peptide', 'charge', 'Score'])

#for i in range(0, 40):
#    out = calculateDotPro(precursor_id = i)
#    results = results.append(out, ignore_index = True)
#    if i == 0:
#        print("Start analysis with first peptide")
#    elif i == max(precursors)*0.3:
#        print("Analyzed 30% of spectra")
#    elif i == max(precursors)*0.5:
#        print("Analyzed 50% of spectra")
#    elif i == max(precursors)*0.8:
#        print("Analyzed 80* of spectra")
#    elif i == max(precursors):
#        print("Analyzed all of spectra")
#
#base_results = "/project/6002011/annieha/pca_urine_spectral_lib/results/"

results.to_csv(os.path.join(base, "20221013_spectral_similarity_fragpipeLibs_batch7.tsv"), sep = "\t", index = False)


#if __name__ == '__main__':
#    with mp.Pool(processes=4) as pool:
#        cal = pool.map(calculateDotPro, precursors[0:4])
#        print("Done")

