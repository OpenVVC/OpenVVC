#!/bin/bash

if [ "$#" -ne 2 ];  then
    echo "Illegal number of parameters"
    echo
    echo "USE ./extractProfileInfo.sh <profiling folder> <output file>"
    exit
fi

PROFILE_FLDR=$1
OUTFILE=$2
for QP in 22 27 32 37;
do
  rm -f $OUTFILE
  for i in $PROFILE_FLDR/*_${QP}_*.cache;
  do
    callgrind_annotate $i --inclusive=yes --threshold=100 > ${i%.cache}.annote;
    filetot=($(grep -E "([0-9,]+) *PROGRAM TOTALS" ${i%.cache}.annote | sed 's/,//g'))
    echo ${filetot[0]} ${filetot[-1]} >> $OUTFILE
    grep -E "^ *([0-9]+,)*[0-9]+.*\?:[^?]*$" ${i%.cache}.annote | while read -r line ;
    do
      match=($(echo $line | sed 's/[,?:]//g'))
      echo ${match[0]} ${match[-1]} >> $OUTFILE
    done
  done

  sort -k 2 $OUTFILE | awk '{a[$2] += $1} END{for (i in a) print a[i], i}' | sort -nr -o $OUTFILE$QP
  cat $OUTFILE$QP
done
