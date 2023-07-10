#!/usr/env python

#import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import msproteomicstoolslib.data_structures as db
import SqlDataAccess
import os
import sys
#from pyopenms import *

# pseudo-code:
# read diann XIC data only take 1 file at a time (no combined file)
diann = sys.argv[1]
dirname = os.path.dirname(diann)
xics = sys.argv[2]
filename = sys.argv[3]
diann = pd.read_csv(diann, sep = "\t")
xics = pd.read_csv(xics, sep = "\t")
xics = xics.reset_index(drop = True)
plotting = sys.argv[4]
plottingIDs = sys.argv[5]

# read peptides for the plots
plottingPeps = open(plottingIDs, "r").read().split('\n')

# rename all the column names
diann = diann.rename(columns = {"File.Name":"Filename", "Precursor.Id":"PrecursorId", "RT.Start":"RTstart", "RT.Stop":"RTend"})
# subset diann output for only necessary columns
diann = diann.loc[:,["Filename", "PrecursorId", "RTstart", "RTend", "RT"]]

xics = xics.rename(columns = {"Precursor.Id": "PrecursorId", "File.Name": "filename", "Retention.Times": "RetentionTimes", "Theoretical.Mz":"TheoreticalMz", "MS.Level": "MSLevel", "Modified.Sequence":"ModifiedSequence", "Stripped.Sequence":"Sequence"})
runs = xics.filename.unique()
len(runs)
# separate xics by file
#xics = xics.loc[xics.filename == runs[0]] 
# data columns start from file1.columns[14:600]

# annotate Fragment ions/transitions
xics = xics.astype({"FragmentSeriesNumber":'str'})
xics['PrecursorCharge'] = xics.PrecursorId.str[-1].astype("int")
conditions = [xics["MSLevel"] == 1, xics["MSLevel"] == 2]
values = [xics["PrecursorId"] + "_Prec", xics["PrecursorId"] + "_" + xics['FragmentSeriesNumber']+xics['FragmentType']]
xics['FragmentAnnotation'] = np.select(conditions, values)

precursors = xics.PrecursorId.unique()
transitions = xics.FragmentAnnotation.unique()


# for each precursor
# get precursorID, get RT start and end from diann
# filter the list to the rt start and end
# find out the number of intensties that are non-zero`
# get values from each row
PointsPerPeak = pd.DataFrame(columns = ('PrecursorId', "RTStart", "RTStop", "FragmentAnnotation", "Charge", "Points"))
trRT = xics.loc[xics.RetentionTimes == 1]
trInt = xics.loc[xics.Intensities == 1]
trID = xics[["filename", "TheoreticalMz", "MSLevel", "ModifiedSequence", "Sequence", "PrecursorId","PrecursorCharge", "FragmentType", "FragmentSeriesNumber", "FragmentAnnotation", "FragmentCharge"]].drop_duplicates()

# start loop
if plotting:
    for i in range(0, len(plottingPeps)):
        p = plottingPeps[i]
        print(p)
        tr = trID.loc[trID.PrecursorId == p]
        pep = tr.ModifiedSequence.unique()[0]
        mz = tr.loc[trID.MSLevel == 1].TheoreticalMz.values[0]
        charge = tr.PrecursorCharge.unique()[0]
        trRank = trInt.loc[trInt.PrecursorId == p]
        #trRank["SumIntensity"] = trRank.iloc[:, 14:615].sum(axis="columns").values
        #trRank['Rank'] = trRank['SumIntensity'].rank(method='dense', ascending=False)
        #trRank["Top6"] = (trRank.Rank <= 6)
        #trRank = trRank.drop_duplicates(subset=["FragmentAnnotation", "Rank", "Top6"], keep=False)
        tr = tr.loc[tr.MSLevel == 2]
        tr = tr.reset_index(drop = True)
        tr = tr.FragmentAnnotation
        if p in diann.PrecursorId.values:
            rtstart = diann.loc[diann.PrecursorId == p].RTstart.values[0]
            rtend = diann.loc[diann.PrecursorId == p].RTend.values[0]
            print("plotting chromatogram for precursor %s" %p)
            fig, ax = plt.subplots(figsize=(6,5))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            #plt.rcParams["font.size"] = "20"
            for k in range(0, len(tr)):
                tID = tr[k]
                #top6 = trRank.loc[trRank.FragmentAnnotation == tID].Top6.values[0]
                #c = 'r' if top6 else 'k'
                rt = trRT.loc[trRT.FragmentAnnotation == tID]
                rt = rt.iloc[:, 14:615].values[0]
                rt0 = (np.abs(rt - rtstart)).argmin() # find the closest number
                rt1 = (np.abs(rt - rtend)).argmin()
                Int = trInt.loc[trInt.FragmentAnnotation == tID]
                Int = Int.iloc[:, 14:615].values[0]
                Int = Int[rt0-2:rt1+2]
                ax.plot(rt[rt0-2:rt1+2], Int, label=tID)
            ax.set_ylabel("Intensity")
            ax.set_xlabel("Retention time (min)")
            fig.show()
            fig.savefig(os.path.join(dirname,"%s_XIC.pdf" %p), dpi = 600)
            plt.close(fig)

for i in range(0, len(precursors)):
    p = precursors[i]
    tr = trID.loc[trID.PrecursorId == p]
    pep = tr.ModifiedSequence.unique()[0]
    mz = tr.loc[trID.MSLevel == 1].TheoreticalMz.values[0]
    charge = tr.PrecursorCharge.unique()[0]
    tr = tr.reset_index(drop = True)
    tr = tr.FragmentAnnotation
    if p in diann.PrecursorId.values:
        rtstart = diann.loc[diann.PrecursorId == p].RTstart.values[0]
        rtend = diann.loc[diann.PrecursorId == p].RTend.values[0]
        for j in range(0, len(tr)):
            tID = tr[j]
            rt = trRT.loc[trRT.FragmentAnnotation == tID]
            rt = rt.iloc[:, 14:615].values[0]
            rt0 = (np.abs(rt - rtstart)).argmin() # find the closest number
            rt1 = (np.abs(rt - rtend)).argmin()
            Int = trInt.loc[trInt.FragmentAnnotation == tID]
            Int = Int.iloc[:, 14:615].values[0]
            Int = Int[rt0:rt1]
            points = sum(Int > 0)
            df = {'PrecursorId':p,  'RTStart':rtstart, 'RTStop':rtend, 'FragmentAnnotation':tID, 'Charge':charge, 'Points':points}
            PointsPerPeak = PointsPerPeak.append(df, ignore_index = True)


PointsPerPeak.to_csv(filename, sep = "\t", index=False)    

# end loop
# save file

# plotting is optional

