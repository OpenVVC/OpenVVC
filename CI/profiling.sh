#!/bin/bash

if [ "$#" -ne 2 ];  then
    echo "Illegal number of parameters"
    echo
    echo "USE ./profiling.sh <bitstreams folder> <decoder> [<output folder>]"
    exit
fi

STREAM=$1
DECODER=$2
OUTDIR=$3

for file in $STREAM/*.266
do
  OUTDIR=$(dirname $file)
  mkdir -p $OUTDIR
  echo $file
  yuv=${file%.266}.yuv
  outfile=$OUTDIR/$(basename ${file%.266})
  valgrind --tool=cachegrind  --cachegrind-out-file=$outfile.cache $DECODER -i $file -o /dev/null
  rm -f $yuv
done

exit 0
