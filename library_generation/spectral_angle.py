########## Python script for computing spectral similarity #########
# This script is used for improving speed from R computing

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mpu
import math

import sys
import os

from sklearn.metrics.pairwise import cosine_similarity

############## Functions ##########################

def Createbins(allSpecdf, pid, bin_width):
    '''
    allSpecdf	dataframe/table	
    pid 	precursor id
    bin_with    float
    Returns dataframe  bins based on the ProductMz range of the precursor spectrum.
    '''
    specs = allSpecdf.loc[allSpecdf['precursor_id'] == pid]
    FragMz = specs['ProductMz'].to_numpy()
    upper = math.ceil(max(FragMz)) + bin_width # upper boundary
    lower = math.floor(min(FragMz)) # lower boundary
    bins = np.arange(lower, upper, bin_width)
    specs = specs.assign(bin_id=np.digitize(FragMz, bins, right = False))
    # assign bin weights
    bin_weights = np.zeros(len(bins), dtype=float)
    bin_weights[:] = 1
    for i in specs.bin_id.unique():
        w = specs.loc[specs['bin_id'] == i, "LibraryIntensity"].to_numpy() #/max(specs['LibraryIntensity'].to_numpy())
        if len(w) > 1:
            bin_weights[i] = max(w)
        else:
            bin_weights[i] = w
    bins = pd.DataFrame({'bins':bins, 'weight':bin_weights})
    return(bins)

def BinSpectrum(spectrumdf, bins):
    '''
    spectrumdf	dataframe	peptide-fragment info (1 peptide)
    bin		np array	list of bins
    Returns binned spectra as numpy array.
    '''
    binned_spectra = np.zeros(len(bins))
    spectrumdf = spectrumdf.assign(binned = np.digitize(spectrumdf['ProductMz'].to_numpy(), bins))
    binned_spectra[spectrumdf.binned.to_numpy()] = spectrumdf['LibraryIntensity'].to_numpy() #/ max(spectrumdf['LibraryIntensity'].to_numpy())
    return(binned_spectra) 

def NormDotP(spec1, spec2, weights):
    """
    spec1	binned spectrum 	peptide-fragment spectrum
    spec2	binned spectrum 	peptide-fragment spectrum
    Returns the cosine similarity of the same peptide between 2 spectra.
    """
    spec1 = np.multiply(spec1, weights)
    dp = cosine_similarity([spec1], [spec2])
    #cos_sim = dot(spec1, spec2)/(norm(spec1)*norm(spec2))
    #dp = 1 - spatial.distance.cosine(spec1, spec2)
    return(dp)

def spitTime(func, params):
    wrapped = wrapper(func, params)
    time = timeit.timeit(wrapped, number = 1)
    print(time)
    return(time)


################ Load data ##########################

# input data should include all spectral information from OpenMS

# Called from another script

