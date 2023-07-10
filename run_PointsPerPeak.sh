#!/bin/sh
#SBATCH --account=def-tkisling
#SBATCH --time=0-07:59
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M

source /home/annieha/bin/python36/bin/activate

dir=/project/6002011/annieha/pca_urine_spectral_lib/data/fragpipe/2020_MethodOpt_XICs
src=/home/annieha/source/tk_orbitrap_dia

files=$(find $dir -type d -name "*lowAbundant*")

for f in $files; do
   echo $f
   diann=$(find $f -type f -name "diann-output.tsv")
   echo $diann
   xics=$(find $f -type f -name "diann-output.XIC.tsv")
   peps=$(find $f -type f -name "*peptidesforXICs.txt")
   #echo $xics
   runcommand="python $src/getPointsPerPeakDIANN.py $diann $xics $f/points.tsv True $peps"
   echo $runcommand
   $runcommand
done
