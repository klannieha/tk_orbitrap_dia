#!/bin/sh

# This script is meant to facilitate th euse of the OpenSwathWorkflow executable


buildpath=/home/klannie/Desktop/codes/OpenMS-build/bin
OPENMS_DATA_PATH=/home/klannie/Desktop/codes/OpenMS/share/OpenMS/


#############################################
# Commandline options
#############################################

inputfile=$1
library=$2
irt_lib=$3
output=$4
doChrom=$5
threads=$6
user=$7
additonalparam=$8

echo $inputfile
echo $library
echo $irt_lib
echo $output
echo $user
echo `pwd`

# Usage

if [ "$user" == "" ]; then
	echo "Not all input arguments were specified, aborting... "
	echo ""
	echo "Complete workflow to run OpenSwath on DIA data"
		echo "--------------------------------------------" ;
		echo "Executes the OpenSwathWorkflow build"
		echo ""
		echo "USAGE INFORMATION:";
		echo "";
		echo "run_osw_Orbitrap.sh inputfile library irt_lib output doChrom threads user advanced";
		echo "";
		echo "inputfile Input fils separated by blank (valid formats: 'mzML', 'mzXML', 'sqMass')"
		echo "library Spectral library containing at least 6 transitions per peptides (valid formats: 'tsv', 'traML', 'pqp')"
		echo "irt_lib iRT spectral library containing retention time and precursors (valid formats: 'tsv', 'traML', 'pqp')"
		echo "output OSW output file (PyProphet compatible SQLite file) (valid formats: 'tsv' (only precursor scoring), 'osw')"
		echo "doChrom OSW output file of chromatographic traces (valid inputs: True, False)"
		echo "threads Number of cores to use"
		echo "user Name of user account"
		echo "advanced Additional parameters in the format \"-parameter value\". Can be multiple parameters not specified in this script."
		echo ""
		exit 1

fi

#########################################
# Helper functions
#########################################

err() {
	echo "$1...exciting";
	exit 1; 
}

ckFile() {
	if [ ! -e "$1"]; then
		err "$2 File '$1' not found";
	fi
}

ckFileSz() {
	ckFile $1 $2
	SZ=`ls -l $1 | awk '{print $5}'`;
	if [ "$SZ" == "0"]; then
		err "$2 file '$1' is zero length";
	else
		echo "$2 file '$1' checked";
	fi
}

ckFileW() {
	if [ ! -w "$1"]; then
		err "$2 File '$1' not writeable";
	fi
}


###########################################
# Check input parameters
###########################################

# Exit if a file the user provided doesn't exist
echo ""
ckFileSz $inputfile "Input file"
ckFileSz $library "Transition file"
ckFileSz $irt_lib "irt transition file"
touch $output
ckFileW $output "Output files"
rm $output
echo "The output will be stored in $output"
echo ""

############################################
# Parameter setup
###########################################

threads=4

dir=`pwd`
echo "Current directory: $dir "

TMPDIR=/home/klannie/tmp

echo " Temporary directory used: $TMPDIR "

OUTCHROM=""
if [[ "$doChrom" == "True" ]]
then
	OUTCHROM="-out_chrom ${output}.sqMass"
	echo "Chromatogram output is enabled"
fi

ADVANCED="-use_ms1_traces -rt_extraction_window 600 -extra_rt_extraction_window 100 -mz_correction_function quadratic_regression_delta_ppm" # !! changeit!!

EXTRACTION_WINDOW="-mz_extraction_window 20 -mz_extraction_window_unit ppm -mz_extraction_window_ms1 20 -mz_extraction_window_ms1_unit ppm -irt_mz_extraction_window_unit ppm -irt_mz_extraction_window 20"

#CALIBRATION="-Calibration:debug_mz_file ${output}_debug_ppmdiff.txt"

SCORING="-Scoring:TransitionGroupPicker:min_peak_width 5.0 -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length 11 -Scoring:stop_report_after_feature 5 -Scoring:TransitionGroupPicker:recalculate_peaks true -Scoring:Scores:use_mi_score -Scoring:Scores:use_ms1_mi"

ADDITIONAL_PARAM=$additionalparam

#############################################
# Run OpenSwath
############################################

# Change output format if needed #

runcommand="$buildpath/OpenSwathWorkflow \
	-in "$inputfile" \
	-tr $library \
	-tr_irt $irt_lib \
	$ADDITIONAL_PARAM \
	$SCORING $ADVANCED \
	-RTNormalization:estimateBestPeptides \
	-Debugging:irt_mzml ${output}_debug_irt.mzML \
	$EXTRACTION_WINDOW -batchSize 250 \
	-readOptions cache -out_osw ${output}.osw $OUTCHROM \
        -tempDirectory $TMPDIR \
	-debug 0 \
	-threads $threads \
	-force"
echo "The run command was: "
echo $runcommand

# Execute OpenSwath
$runcommand
ls -ltrh $TMPDIR


##############################################
# Finish up
##############################################
echo "done running"
rm $TMPDIR/*.mzML.cache* 

rm $TMPDIR/*.mzML


