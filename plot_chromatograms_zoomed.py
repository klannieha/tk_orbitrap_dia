#!/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import chrom_functions as ch
import numpy as np
import os
import sys
from scipy import signal
from scipy.ndimage import gaussian_filter1d

#plt.rcParams["figure.figsize"] = (8,4)
peptides = sys.argv[1]
peptides = open(peptides, "r").read().split('\n')

sqmass = sys.argv[2]

pyprophet = sys.argv[3]
pyprophet = pd.read_csv(pyprophet, sep = "\t")
#sqmass='/project/6011811/frankmak/tims/Mar2018_helaBenchmark/data/openswath/20201020_NewiRT_Analysis/2019_May200ng/'
#IMNONE_sqmass= sqmass + '/20190507_TIMS1_FlMe_SA_HeLa_diaPASEF_py3_1_A1_1_119_0_IMNONE.sqMass'
#IM6_sqmass= sqmass + '/20190507_TIMS1_FlMe_SA_HeLa_diaPASEF_py3_1_A1_1_119_0_IM6.sqMass'
#IM9_sqmass= sqmass + '/20190507_TIMS1_FlMe_SA_HeLa_diaPASEF_py3_1_A1_1_119_0_IM9.sqMass'

suffix = sys.argv[4]
#suffix = "IM6"
filename = sys.argv[5]
#pep = peptides[0]
#chrom1 = ch.getChromValues(pep, sqmass)

# plot
# test
#fig, ax = plt.subplots()

#for i in range(0, len(chrom1[1])):
#    if "Precursor" not in chrom1[1][i]:
#        ax.plot(chrom1[3][i], chrom1[4][i], label=chrom1[1][i]) 
#ax.legend()

#fig.savefig("%s_test.png"%pep)
#plt.close(fig)

for p in peptides:
    print("Plotting Chromatogram for peptide %s" %p)
    fig, ax = plt.subplots(figsize=(12, 4))
    chrom0 = ch.getChromValues(p, sqmass)
    chrom_id = {"native_id" : chrom0[1]}
    print(chrom_id)
    tp = pyprophet.loc[(pyprophet['Sequence'] == p) & (pyprophet['filename'] == filename), ['leftWidth', 'rightWidth']].values.tolist()[0]
    if len(chrom_id['native_id']) != 0:
        chrom_id = pd.DataFrame(chrom_id)
        chrom_id = chrom_id[-chrom_id['native_id'].str.contains("Precursor")]
        chrom_id = chrom_id.sort_values('native_id')
        idx = list(chrom_id.index)
        allINt = []
        for i in idx:
            for j in chrom0[4][i]:
                allINt.append(j)
        #allINt = [j for i in chrom0[4] for j in i]
        #print(allINt)
        maxInt = max(allINt)
        minInt = 0
        #minInt = min(allINt)
        print(maxInt)
        for i in idx:
            x = chrom0[3][i]
            y = chrom0[4][i]
            print(y)
            #filter = gaussian_filter1d(y, sigma = 9, mode = "constant")
            #if max(y) > 0:
            #ref = [(Int/maxInt)*100 for Int in y]
            scaled = [((Int - minInt)/(maxInt - minInt))*100 for Int in y]
            filter = gaussian_filter1d(y, sigma = .5)
            x = np.array(x)
            scaled = np.array(x)
            midx = np.where((x > tp[0]) & (x < tp[1]))
            x_zoomed = x[midx]
            ax.plot(x_zoomed, scaled[midx], label=chrom0[1][i])
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_title(p)
        ax.set_ylim([0, 110])
        ax.set_ylabel("Relative Intensity (%)")
        ax.set_xlabel("Retention time (s)")
        fig.show()
        fig.savefig("%s_%s_zoomed.png" % (p,suffix), dpi=600)
        plt.close(fig)
    else:
        print("Peptide not found.")
