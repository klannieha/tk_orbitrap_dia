#!/usr/env/ python

############ THIS IS A CUSTOM SCRIPT ##################################
# Purpose:
# 	With input files of a pyprophet output tsv with respect to its 
#	sqMass file, calculate number of non-zero points across the 
#	detected peptide/precursor/transition. Output a tsv file.
#	This script is to be used as a commandline tool.

# Required arguments:
#	Input tsv: pyprophet file
#	Input filename: filename of the mzml file
#	Input sqMass: the corresponding sqMass file containing the XICs
#	Output path: output path for tsv file

#######################################################################


import sqlite3
import numpy as np
import pandas as pd
import matplotlib as plt
import msproteomicstoolslib.data_structures as db
import SqlDataAccess
import os
import sys
import time

filename=sys.argv[1]
pyprophet=sys.argv[2]
sqMass=sys.argv[3]
output=sys.argv[4]

# Subset the data
stats = pd.read_csv(pyprophet, sep = "\t")
stats = stats.loc[:, ['transition_group_id', 'filename', 'RT', 'mz', 'Sequence', 'decoy','id','FullPeptideName', 'Charge', 'm_score', 'leftWidth', 'rightWidth', 'peak_group_rank', 'd_score']]

stats = stats[stats.decoy == 0] # remove decoys
stats = stats[stats.m_score <= 0.01] # filter by 1% FDR
stats = stats[stats.filename == filename]
# Rank by m_score (fdr)
stats = stats.loc[stats.groupby(["FullPeptideName", "Charge"])['m_score'].idxmin()]
stats = stats[stats.peak_group_rank == 1]
# Remove duplicates
stats = stats.drop_duplicates(subset=['FullPeptideName', 'Charge'])


# Modified Functions from DIAPASEF scripts:
def get_PeakWidth(df):
    '''
    Given the dataset, calculate the peakwidth for each peptide.
    Input: pyprophet pandas dataframe
    Output: dataframe of peakwidth per peptide.
    '''
    pw = df.loc[:, ['Sequence', 'FullPeptideName', 'Charge', 'leftWidth', 'rightWidth']]
    width = df.rightWidth - df.leftWidth
    pw['PeakWidth'] = width
    return(pw)

def readSQL(sql):
    conn = sqlite3.connect(sql)
    return(conn)

def getChromData(sql):
    chrom = SqlDataAccess.SqlDataAccess(sql)
    return(chrom)

def getChromValues(peptide, charge, chrom, conn):
    start = time.time()
    query = '''SELECT * FROM CHROMATOGRAM \
    INNER JOIN DATA ON DATA.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    INNER JOIN PRECURSOR ON PRECURSOR.CHROMATOGRAM_ID = CHROMATOGRAM.ID \
    where PEPTIDE_SEQUENCE="{0}" AND CHARGE="{1}"'''
    q = conn.execute(query.format(peptide, charge))
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
    selected = chrom.getDataForChromatograms(idx)
    # selected: nested list
	# list:
	#	- chrom id
	#		- rt
	#		- int
   # for t in nativeid:
   #     selected = chrom.getDataForChromatogramFromTransitionNativeId(t)
   #     rt.append(selected[0])
   #     intensity.append(selected[1])
    end = time.time()
    print("Time elapsed: ", end-start)
    return(idx,nativeid, peptide, selected)
	
def getPointsPerChrom(pw, peptide, chrom):
    p = peptide
    tp = pw.loc[pw['Sequence'] == p, ['leftWidth', 'rightWidth']].values.tolist()[0]
    charge = pw.loc[pw['Sequence'] == p, "Charge"].values.tolist()[0]
    out = pd.DataFrame(columns = ('Sequence','lw','rw', 'ChromID', 'Points', "Charge"))
    chromID = chrom[0]
    chromINFO = chrom[3]
    for chrID in range(0, len(chrom[3])):
        tmp_rt = np.array(chromINFO[chrID][0])
        tmp_int = np.array(chromINFO[chrID][1])
        if (max(tmp_rt) > tp[1]) | (min(tmp_rt) < tp[0]): # tp[0] is the lw, tp[1] is rw
            idx = np.where((tmp_rt > tp[0]) & (tmp_rt < tp[1]))
            n_int = tmp_int[idx]
            n_point = sum(n_int > 0)
            out.loc[chrID] = [p, tp[0], tp[1], chromID[chrID], n_point, charge]
    return(out)           
#    for i in range(0, len(chrom[0])):
#        tmp_rt = np.array(chrom[3][i])
#        #print("length of rt: %d" % len(tmp_rt))
#        tmp_int = np.array(chrom[4][i])
#        #print("length of int: %d" %len(tmp_int))
#        if (max(tmp_rt) > tp[1]) & (min(tmp_rt) < tp[0]):
#            idx = np.where((tmp_rt < tp[0]) | (tmp_rt > tp[1]))
#           # print("length of idx: %d" %len(idx))
#            n_int = tmp_int[idx]
#            midx = np.where((tmp_rt > tp[0]) & (tmp_rt < tp[1]))
#            #print("length of max idx: %d" %len(midx))
#            n_points = sum(n_int > 0)
#            out.loc[i] = [p, chrom[0][i], chrom[1][i], n_points, tp[0], tp[1]]
#    return(out)

# Use functions to calculate for number of peaks per peptide:
chrom = getChromData(sqMass)
conn = readSQL(sqMass)

Peptide_peakWidths = get_PeakWidth(stats)

peptides =  stats.Sequence.tolist()

#Peptide_chroms = [getChromValues(pep, chrom, conn) for pep in peptides]

PointsPerPeak = pd.DataFrame()

for p in range(0, len(peptides)):
    pep = Peptide_peakWidths.Sequence.tolist()[p]
    charge = Peptide_peakWidths.Charge.tolist()[p]
    print(pep, charge)
    chromValue = getChromValues(pep, charge, chrom, conn)
    n_Points = getPointsPerChrom(Peptide_peakWidths, pep, chromValue)
    PointsPerPeak = PointsPerPeak.append(n_Points, ignore_index=True)
    


PointsPerPeak.to_csv(output, sep = "\t", index=False)


