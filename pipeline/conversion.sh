#!/bin/sh

echo `pwd`
dir=$1
in=$2
outfile=$3
ADDITIONAL=$4

if [ -f "$in" ]; then
	echo "$in exists"
else
	echo "$in does not exist"
fi

echo "mzML file will be written into $outfile"

pwiz=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

COMPRESS_OPTIONS="-z --numpressLinear"

runcommand="docker run -it --rm -e WINEDEBUG=-all -v $dir:/data $pwiz wine msconvert $in --outfile=$outfile $COMPRESS_OPTIONS $ADDITIONAL"

echo $runcommand

$runcommand
