#!/usr/env python

import diapysef as dp
import pandas as pd
import numpy as np
import os
import sys

dir="/home/klannie/Desktop/dia_pipeline/data/DIA"

mqdir=os.path.join(dir, "TXT")

evidence=pd.read_csv(os.path.join(mqdir, "evidence.txt"), sep = "\t")
msms = pd.read_csv(os.path.join(mqdir, "msms.txt"), sep = "\t")

irt = pd.read_csv("/home/klannie/Desktop/dia_pipeline/data/libraries/irt/CiRT_ALL.tsv", sep = "\t")

mqout = dp.pasef_to_tsv(evidence, msms, irt_file = irt, ion_mobility = False, im_alignment = False, rt_alignment = "nonlinear", pdfout="rtcalibration.pdf")

