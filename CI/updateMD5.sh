#!/usr/bin/env bash

if [ "$#" -ne 2 ];  then
    echo "Illegal number of parameters"
    echo
    echo "USE ./updateMD5.sh <bitstreams folder> <decoder>"
    echo
    exit
fi

STREAM=$1
DECODER=$2

for file in $1/*.266
do
  yuv=${file%.266}.yuv
  $DECODER -b $file -o rec_orig.yuv --SEIFGSFilename=$yuv --OutputDecodedSEIMessagesFilename=sei.txt

  # $DECODER -b $file -o $yuv
  md5sum $yuv > ${file%.266}.md5
  rm -f $yuv
  rm -f rec_orig.yuv
done
