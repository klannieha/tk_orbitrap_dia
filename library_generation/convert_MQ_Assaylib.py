#!/usr/env python

################## Convert MQ format to Assay lib format #########
# this script is to reformat the merged EV and MSMS files from MQ
# into fragment ion spectral library

# input is file merged with msms and ev with NormalizedRetentionTime

import sys
import os
import matplotlib
import pandas as pd
import statsmodels.formula.api as sm
import numpy as np

ms = sys.argv[1]
out = sys.argv[2]

ms = pd.read_table(ms, sep = "\t")
ms = ms.loc[ms.groupby(['ModifiedPeptideSequence', 'Charge'])['PEP'].idxmin()]

ms = ms.loc[:, ["id", "m/z", "Masses", "Charge", "NormalizedRetentionTime", "Intensities", "PeptideSequence", "ModifiedPeptideSequence", "Proteins", "UniProtID"]]

masses = ms['Masses'].str.split(';', expand = True).stack().str.strip().reset_index(level = 1, drop = True)
intensities = ms['Intensities'].str.split(';', expand = True).stack().str.strip().reset_index(level = 1, drop=True)
df = pd.concat([masses, intensities], axis = 1, keys = ['Masses', 'Intensities'])
mslib = ms.drop(['Masses', 'Intensities'], axis = 1)
mslib = mslib.join(df).reset_index(drop = True)

mslib.columns = ["transition_group_id","PrecursorMz","PrecursorCharge","NormalizedRetentionTime", "PeptideSequence","FullUniModPeptideName","ProteinName", "UniProtID","ProductMz", "LibraryIntensity"]

mslib['decoy'] = 0

mslib['transition_name'] = ['_'.join(str(i) for i in z) for z in zip(mslib.transition_group_id,mslib.PrecursorMz,mslib.ProductMz)]
mslib=mslib.dropna(subset=['LibraryIntensity'])
mslib=mslib.drop_duplicates()


# write output
mslib.to_csv("%s.tsv" % out ,sep = "\t", index=False)
