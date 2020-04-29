#!/bin/sh

#SBATCH --account=def-tkisling
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=2
#SBATCH --time=0-07:59

basefolder=/project/6002011/annieha/pca_urine_spectral_lib/src/
cd $basefolder

source /home/annieha/bin/python36/bin/activate

python -u run_DotProduct.py > log.txt



