#!/usr/bin/env bash

STREAM=$1
DECODER=$2

for file in $1/*.266
do
  yuv=${file%.266}.yuv
  $DECODER -i $file -o $yuv
  md5sum $yuv > ${file%.266}.md5
  rm -f $yuv
done
