#!/bin/python
# This script is to generate input FragPipe manifest files with
# bootstrapped filenames for library generation

from numpy import *
from pandas import *
import random
import os

# get filenames

files = os.listdir("/scratch/annieha/tklab/data/dda/urine")

# paste full paths
filenames = []
for i in range(0, len(files)):
    filenames.append(os.path.join("/scratch/annieha/tklab/data/dda/urine", files[i]))

##### Start bootstrapping ##################################
nums = [i for i in range(len(files))] # 199 files = 199 integers to pick from

random.seed(123)

# Bootstrap for 50 times per increment

increments = [i for i in range(2, 20, 2)]
for i in range(20, len(files), 20):
    increments.append(i)

for x in range(0, len(increments)):
    size = increments[x]
    print("Writing bootstrap files for size of :", size)
    iterations = 0
    while iterations < 25:
        idx = random.sample(nums, size)
        f = open(f"./manifests/manifest_sample{size}_bootstrap{iterations}.txt", "w")
        for y in range(0, len(idx)):
            quote = filenames[idx[y]] + "\t" + "DDA" + "\n"
            f.write(quote)
        f.close()
        iterations = iterations + 1



