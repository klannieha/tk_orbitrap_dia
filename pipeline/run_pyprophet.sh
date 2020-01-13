#!/bin/sh

source /home/klannie/Desktop/codes/pyven_37/bin/activate

basefolder=/home/klannie/Desktop/dia_pipeline

data=$basefolder/results/ankit_HEKqc

input=$(find $data/ -maxdepth 1 -type f -name "*.osw")

for f in $input; do
	filename=$(basename -- "$f")
	filename="${filename%.*}"
	pyprophet score --classifier=XGBoost --level=ms2 --in=$f
	pyprophet peptide --context=run-specific --in=$f
	pyprophet protein --context=run-specific --in=$f
	pyprophet export --in=$f --out=${filename}_pyprophet_export.tsv --format=legacy_merged --ipf=disable
	sleep 1s
done

