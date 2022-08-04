# This script makes tryptic peptides only
# source python36 from bin
from pyopenms import *
import pandas as pd
import numpy as np
import sys
import os
from itertools import zip_longest

# load the uniprot database with all the protein accessions
fastaFile = sys.argv[1]
outFile = sys.argv[2]

# format output table
dicts = {}

proteins = list() # make protein objects in the list
# entries here are FASTAEntry objects

# check length of the proteins
FASTAFile().load(fastaFile, proteins) # load the protein accessions to openms objects
len(proteins)

# Get sequence with only fully digested tryptic peptides
dig = ProteaseDigestion()
dig.setEnzyme('Lys-C')
dig.setMissedCleavages(0)

results = []

# Get a list of sequence objects from 1 protein ID

for p in proteins:
    ProteinID = p.identifier
    s = AASequence.fromString(p.sequence)
    res = []
    dig.digest(s, res, 10, 100)
    res_str = [x.toString() for x in res]
    dicts[ProteinID] = res_str

df = pd.DataFrame.from_dict(dicts, orient = 'index').T.melt(var_name="ProteinName", value_name = "Sequence").dropna(subset=['Sequence'])

df.to_csv(outFile, sep = "\t", index = False)

