#!/bin/sh

#SBATCH --account=def-tkisling
#SBATCH --mem-per-cpu=16000M
#SBATCH --cpus-per-task=4
#SBATCH --time=0-00:59

dir=`pwd`
cd $dir

source /home/annieha/bin/python36/bin/activate

python -u run_DotProduct.py > log.txt



