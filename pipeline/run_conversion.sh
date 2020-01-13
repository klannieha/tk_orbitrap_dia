#!/bin/sh

dir=/home/klannie/Desktop/dia_pipeline/data/DIA/raw
outdir=/home/klannie/Desktop/dia_pipeline/data/DIA/mzML
input=$(find $dir -maxdepth 1 -type f -name "*dynamic*.raw")
src=/home/klannie/Desktop/dia_pipeline/src

for f in $input; do
	filename=$(basename -- "$f")
	outfile="${filename%.*}"
	echo $filename
	. $src/conversion.sh $dir $filename ${outfile}.mzML
done

mv $dir/*.mzML $outdir
