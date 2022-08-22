#!/usr/env python

import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import msproteomicstoolslib.data_structures as db
import SqlDataAccess
import os
import sys
from pyopenms import *

# pseudo-code:
# read diann XIC data only take 1 file at a time (no combined file)
diann = sys.argv[1]
filename = sys.argv[2]
xics = pd.read_csv(diann, sep = "\t")
xics = xics.reset_index(drop = True)
# rename all the column names
xics = xics.rename(columns = {"Precursor.Id": "PrecursorId", "File.Name": "filename", "Retention.Times": "RetentionTimes", "Theoretical.Mz":"TheoreticalMz", "MS.Level": "MSLevel", "Modified.Sequence":"ModifiedSequence", "Stripped.Sequence":"Sequence"})
runs = xics.filename.unique()

# separate xics by file
xics = xics.loc[xics.filename == runs[0]] 
# data columns start from file1.columns[14:600]

# annotate Fragment ions/transitions
xics = xics.astype({"FragmentSeriesNumber":'str'})
xics['PrecursorCharge'] = xics.PrecursorId.str[-1].astype("int")
conditions = [xics["MSLevel"] == 1, xics["MSLevel"] == 2]
values = [xics["PrecursorId"] + "_Prec", xics["PrecursorId"] + "_" + xics['FragmentSeriesNumber']+xics['FragmentType']]
xics['FragmentAnnotation'] = np.select(conditions, values)

# store as mzML
exp = MSExperiment()

# use pyopenms MSChromatogram()
# assign peaks
# chromatogram.set_peaks([rt,i])

# set loop
# for each precursor Id:
#   set peaks
#   set Precursor

# store as mzML in the end
#exp = MSExperiment()


# get values from each row

precursors = xics.PrecursorId.unique()
transitions = xics.FragmentAnnotation.unique()

trRT = xics.loc[xics.RetentionTimes == 1]
trInt = xics.loc[xics.Intensities == 1]
trID = xics[["filename", "TheoreticalMz", "MSLevel", "ModifiedSequence", "Sequence", "PrecursorId","PrecursorCharge", "FragmentType", "FragmentSeriesNumber", "FragmentAnnotation", "FragmentCharge"]].drop_duplicates()

# start loop

for i in range(0, len(precursors)):
    p = precursors[i]
    tr = trID.loc[trID.PrecursorId == p]
    pep = tr.ModifiedSequence.unique()[0]
    mz = tr.loc[trID.MSLevel == 1].TheoreticalMz.values[0]
    charge = tr.PrecursorCharge.unique()[0]
    tr = tr.reset_index(drop = True)
    tr = tr.FragmentAnnotation
    Prec = Precursor()
    Prec.setMZ(mz)
    Prec.setMetaValue("peptide_sequence", pep)
    Prec.setCharge(int(charge))
    for j in range(0, len(tr)):
        tID = tr[j]
        chrom = MSChromatogram()
        rt = trRT.loc[trRT.FragmentAnnotation == tID]
        rt = rt.iloc[:, 14:615].values[0]
        rt = rt * 60 # convert minutes to seconds
        Int = trInt.loc[trInt.FragmentAnnotation == tID]
        Int = Int.iloc[:, 14:615].values[0]
        chrom.set_peaks([rt, Int])
        chrom.setNativeID(tID)
        prod = Product()
        prod.setMZ(trID.loc[trID.FragmentAnnotation == tID].TheoreticalMz.unique()[0])
        chrom.setProduct(prod)
        chrom.setPrecursor(Prec)
        exp.addChromatogram(chrom)
    

# end loop
# save file
MzMLFile().store(filename, exp)

# plotting is optional

