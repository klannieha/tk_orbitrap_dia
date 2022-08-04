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
# read diann XIC data
# store as mzML
exp = MSExperiment()

# use pyopenms MSChromatogram()
chromatogram = MSChromatogram()

# set precursor
# rename all the column names
xics = xics.rename(columns = {("Precursor.Id": "PrecursorId", "File.Name": "filename", "Retention.Times": "RetentionTimes", "Theoretical.Mz":"TheoreticalMz", "MS.Level": "MSLevel", "Modified.Sequence":"ModifiedSequence", "Stripped.Sequence":"Sequence"})

# separate xics by file

# get values from each row

precursors = xics.PrecursorId.unique()

# assign peaks
# chromatogram.set_peaks([rt,i])

# set loop
# for each precursor Id:
#   set peaks
#   set Precursor

# store as mzML in the end
#exp = MSExperiment()
exp.addChromatogram(chromatogram)

# end loop
# save file
MzMLFile().store(filename, exp)

# plotting is optional

