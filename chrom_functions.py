#!/usr/env python

import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import msproteomicstoolslib.data_structures as db
import SqlDataAccess
import os
import sys

def testNativeID(nativeid, filename):
    conn = sqlite3.connect(filename)
    query = '''SELECT * FROM CHROMATOGRAM \
    INNER JOIN DATA ON DATA.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    INNER JOIN PRECURSOR ON PRECURSOR.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    where NATIVE_ID="{0}"'''
    q = conn.execute(query.format(nativeid))
    tmp = q.fetchall()
    return(tmp)

def getChromValues(peptide, filename):
    chrom = SqlDataAccess.SqlDataAccess(filename)
    query = '''SELECT * FROM CHROMATOGRAM \
    INNER JOIN DATA ON DATA.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    INNER JOIN PRECURSOR ON PRECURSOR.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    where PEPTIDE_SEQUENCE="{0}"'''
    conn = sqlite3.connect(filename)
    q = conn.execute(query.format(peptide))
    tmp = q.fetchall()
    idx = []
    nativeid = []
    for i in range(0, len(tmp)):
        idx.append(tmp[i][0])
        nativeid.append(tmp[i][2])
    idx = list(set(idx))
    nativeid = list(set(nativeid))
    rt = []
    intensity = []
    for t in nativeid:
        selected = chrom.getDataForChromatogramFromTransitionNativeId(t)
        rt.append(selected[0])
        intensity.append(selected[1])
    return(idx,nativeid, peptide, rt, intensity)

def subset_df(filename, mzml):
    f = pd.read_csv(filename, sep='\t')
    f = f[f.filename == mzml]
    f = f[f.decoy == 0]
    f = f[f.peak_group_rank == 1]
    f = f[f.m_score <= 0.01]
    f = f.loc[:,  ['transition_group_id', 'run_id', 'filename','RT', 'id', 'Sequence', \
                   'FullPeptideName', 'Charge', 'mz', 'Intensity','ProteinName', 'decoy','leftWidth', 'rightWidth',\
                   'peak_apices_sum', 'd_score', 'm_score', 'peak_group_rank']]
    return(f)

def get_peakwidth(sequence, dataframe):
    ''' sequence -> list of peptide sequence '''
    t = dataframe.loc[dataframe.Sequence.isin(sequence),["leftWidth", "rightWidth", "RT", "Sequence", \
                                                         "m_score"]]
    df = pd.DataFrame()
    for p in sequence:
        tmp = t.loc[t.Sequence == p]
        tmp = tmp.dropna(subset=["leftWidth", "rightWidth", "RT"])
        if tmp.shape[0] > 1:
            df = df.append(tmp.ix[tmp.m_score.idxmin()], ignore_index=True)
        else:
            df = df.append(tmp, ignore_index=True)
#    for i in range(0, len(df)):
#        if np.isnan(df[i][0][1]):
#            rt = df[i][0][2]
#            lw = df[i][0][0]
#            df[i][0][1] = rt - lw + rt
    return(df)

def get_sn(sequence, df, chrom):
    p = sequence
    tp = df.loc[df['Sequence'] == p, ['leftWidth', 'rightWidth']].values.tolist()[0]
    out = pd.DataFrame(columns = ('Sequence', 'ChromID', 'NativeID','MaxInt', 'AvgNoise', 'sn','lw','rw'))
    for i in range(0, len(chrom[0])):
        tmp_rt = np.array(chrom[3][i])
        #print("length of rt: %d" % len(tmp_rt))
        tmp_int = np.array(chrom[4][i])
        #print("length of int: %d" %len(tmp_int))
        if (max(tmp_rt) > tp[1]) & (min(tmp_rt) < tp[0]):
            idx = np.where((tmp_rt < tp[0]) | (tmp_rt > tp[1])) # indices where the rt is outside of the peakwidth
           # print("length of idx: %d" %len(idx))
            n_int = tmp_int[idx] #noise intensities
            midx = np.where((tmp_rt > tp[0]) & (tmp_rt < tp[1])) # max peak intensity
            #print("length of max idx: %d" %len(midx))
            m = max(tmp_int[midx])
            sn = m/n_int.mean()
            #print("Peptide: %s   Transition ID: %s " %(p, chrom[1][i]))
            #print("Max Intensity: %f    Mean Noise: %f  " %(m,  n_int.mean()))
            #print( "S/N: %f" % sn)
            out.loc[i] = [p, chrom[0][i], chrom[1][i], m, n_int.mean(), sn, tp[0], tp[1]]
    return(out)


def plotChrom(sequence, sqmass, legend=False, savefig=False):
    print("Plotting Chromatogram for peptide %s" %sequence)
    fig, ax = plt.subplots()
    chrom0 = getChromValues(sequence, sqmass)
    chrom_id = {"native_id" : chrom0[1]}
    chrom_id = pd.DataFrame(chrom_id)
    chrom_id = chrom_id[-chrom_id['native_id'].str.contains("Precursor")]
    chrom_id = chrom_id.sort_values('native_id')
    idx = list(chrom_id.index)
    allINt = [j for i in chrom0[4] for j in i]
    maxInt = max(allINt)
    print(maxInt)
    for i in idx:
        x = chrom0[3][i]
        y = chrom0[4][i]
        #filter = gaussian_filter1d(y, sigma = 9, mode = "constant")
        #if max(y) > 0:
        ref = [(Int/maxInt)*100 for Int in y]
        filter = gaussian_filter1d(y, sigma = .5)
        ax.plot(x, ref, label=chrom0[1][i])
#    box = ax.get_position()
#    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title(sequence)
    ax.set_ylim([0, 110])
    ax.set_ylabel("Relative Intensity (%)")
    ax.set_xlabel("Retention time (s)")
    fig.show()
    #fig.savefig("%s_%s.png" % (p,suffix), dpi=600)
    return(fig) 


def plotChromWithinPeakWidth(sequence, sqmass, df,legend=False, savefig=False):
    p = sequence
    tp = df.loc[df['Sequence'] == p, ['leftWidth', 'rightWidth']].values.tolist()[0]
    print("Plotting Chromatogram for peptide %s" %sequence)
    fig, ax = plt.subplots()
    chrom0 = getChromValues(sequence, sqmass)
    chrom_id = {"native_id" : chrom0[1]}
    chrom_id = pd.DataFrame(chrom_id)
    chrom_id = chrom_id[-chrom_id['native_id'].str.contains("Precursor")]
    chrom_id = chrom_id.sort_values('native_id')
    idx = list(chrom_id.index)
    allINt = [j for i in chrom0[4] for j in i]
    maxInt = max(allINt)
    print(maxInt)
    for i in idx:
        x = chrom0[3][i]
        y = chrom0[4][i]
        peak = np.where((x > tp[0]) & (x < tp[1]))
        x = x[peak]
        y = y[peak]
        filter = gaussian_filter1d(y, sigma = .5)
        ax.plot(x, y, label=chrom0[1][i])
#    box = ax.get_position()
#    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title(sequence)
    #ax.set_ylim([0, 110])
    ax.set_ylabel("Intensity")
    ax.set_xlabel("Retention time (s)")
    fig.show()
    #fig.savefig("%s_%s.png" % (p,suffix), dpi=600)
    return(fig) 




