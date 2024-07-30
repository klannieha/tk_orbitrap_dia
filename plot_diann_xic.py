#!/usr/env python

#import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
#from pyopenms import *
import re

# pseudo-code:
# read diann XIC data only take 1 file at a time (no combined file)
xics = sys.argv[1]
dirname = os.path.dirname(xics)
xics = pd.read_csv(xics, sep = "\t")
xics = xics.reset_index(drop = True)
#plottingIDs = sys.argv[2]

# read peptides for the plots
#if plotting:
#    plottingPeps = open(plottingIDs, "r").read().split(', ')

#plottingPeps = [re.sub("\n", "", x) for x in plottingPeps]
#plottingPeps = [re.sub(" ", "", x) for x in plottingPeps]

xics = xics.rename(columns = {"Precursor.Id": "PrecursorId", "File.Name": "filename", "Retention.Times": "RetentionTimes", "Theoretical.Mz":"TheoreticalMz", "MS.Level": "MSLevel", "Modified.Sequence":"ModifiedSequence", "Stripped.Sequence":"Sequence"})
xics["basename"] = xics['filename'].apply(os.path.basename)
runs = xics.basename.unique()
len(runs)
precursors = xics.PrecursorId.unique()

# add annotation to the precursors and transitions
xics = xics.astype({"FragmentSeriesNumber":'str'})
xics['PrecursorCharge'] = xics.PrecursorId.str[-1].astype("int")
conditions = [xics["MSLevel"] == 1, xics["MSLevel"] == 2]
values = [xics["PrecursorId"] + "_Prec", xics["PrecursorId"] + "_" + xics['FragmentSeriesNumber']+xics['FragmentType']]
xics['FragmentAnnotation'] = np.select(conditions, values)

# separate xics by file
for f in range(0, len(runs)):
    print(runs[f])
    xics_sub = xics.loc[xics.basename == runs[f]]
    transitions = xics_sub.FragmentAnnotation.unique()
    # for each precursor: extract RT, INT, map together
    for i in range(0, len(precursors)):
        p = precursors[i]
        print(p)
        prec = xics_sub.loc[xics_sub.PrecursorId == p]
        # melt wide RT and int format to long and merge
        s = re.compile("\d+")
        cols = [j for j in prec.columns if s.match(j)]
        rt = prec.loc[prec.RetentionTimes == 1]
        values = prec.loc[prec.Intensities == 1]
        rt = rt.melt(id_vars = ["PrecursorId", "FragmentAnnotation"], value_vars = cols, value_name = "RT", var_name = "ID")
        values = values.melt(id_vars = ["PrecursorId", "FragmentAnnotation"], value_vars = cols, value_name = "intensity", var_name = "ID")
        df = pd.merge(rt, values, on = ["PrecursorId", "FragmentAnnotation", "ID"])
        ids = prec[["filename", "TheoreticalMz", "MSLevel", "ModifiedSequence", "Sequence", "PrecursorId", "PrecursorCharge", "FragmentType", "FragmentSeriesNumber", "FragmentAnnotation", "FragmentCharge"]].drop_duplicates()
        transitions = ids.FragmentAnnotation.unique()
        # assign ranges
        xmin = df['RT'].agg('min') - 2
        xmax = df['intensity'].agg('max') + 2
        # start plot
        print("plotting chromatogram for precursor %s" %p)
        fig, ax = plt.subplots(figsize=(6,5))
        plt.ticklabel_format(style="sci", axis = "y", scilimits = (0,0))
        # start loop for transition
        for k in range(0, len(transitions)):
            tr = transitions[k]
            x = df.loc[df.FragmentAnnotation == tr]
            x = x.RT.values
            y = df.loc[df.FragmentAnnotation == tr]
            y = y.intensity.values
            ax.plot(x, y, label = tr)
        ax.set_ylabel("Intensity")
        ax.set_xlabel("RT (min)")
        fig.show()
        fig.savefig(os.path.join(dirname, "%s_%s_XIC.pdf" % (runs[f], p)), dpi = 360)
        plt.close(fig)

# end loop
# save file


