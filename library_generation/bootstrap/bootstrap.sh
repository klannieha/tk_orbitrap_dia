#!/bin/sh
#SBATCH --account=def-tkisling


# assign working directories

dir=/home/annieha/bin/fragpipe
wdir=`pwd`
output=$wdir/bootstrap_splib

# get all manifest files
manifiles=$(find $wdir/manifests/ -maxdepth 1 -type f -name "manifest_sample199_bootstrap*.txt")


for f in $manifiles; do
    echo "Manifest file: $f"
    output_name=$(basename -- "$f")
    output_name="${output_name##*manifest_}"
    output_name="${output_name%.*}"
    echo "Output directory: $output/$output_name"
    mkdir $output/$output_name
    sbatch $dir/run_fragpipe.sh $wdir/DIA_Speclib.workflow $f $output/$output_name annieha 8 8
    sleep 1s
done
