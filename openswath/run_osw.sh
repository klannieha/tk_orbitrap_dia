#!/bin/sh

basefolder=/home/klannie/Desktop/dia_pipeline
src=$basefolder/src
data=$basefolder/data/DIA/mzML
tr_lib=$basefolder/data/libraries/HEK_qc_mqout_SCX_10only_assaylib_targetdecoy.pqp
#tr_lib=$basefolder/data/libraries/Ankit_HEK_qc_targetdecoy_assaylib.pqp
irt=$basefolder/data/libraries/irt/iRTassays.TraML
dir=`pwd`
input=$(find $data -maxdepth 1 -type f -name "*.mzML")

for f in $input;do
	filename=$(basename -- "$f")
	filename="${filename%.*}"
	log=${filename}.log
	#script $log
	. $src/run_osw_Orbitrap_HEK.sh $f $tr_lib $irt $dir/$filename False 4 klannie \
		"-RTNormalization:alignmentMethod linear -RTNormalization:outlierMethod none -Scoring:TransitionGroupPicker:compute_peak_quality false -ms1_isotopes 3"
        mv debug_ppmdiff.txt ${filename}_debug_ppmdiff.txt	
	sleep 1s
	#exit
done

