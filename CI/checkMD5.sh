#!/bin/bash

STREAM=$1
DECODER=$2

for file in $1/*.266
do
  yuv=${file%.266}.yuv
  $DECODER -i $file -o $yuv
  md5sum $yuv | diff ${file%.266}.md5 -
  error=$?
  rm -f $yuv
  if [ $error -ne 0 ]
  then
    exit $error
  fi
done

exit 0
