#!/bin/bash

STREAM=$1
DECODER=$2

for file in $STREAM/*.266
do
  echo $file
  yuv=${file%.266}.yuv
  $DECODER -i $file -o $yuv
  MD5=$(md5sum $yuv | grep -o '[0-9,a-f]*\ ')
  MD5fc=$(cat ${file%.266}.md5 | grep -o '[0-9,a-f]*\ ')
  echo $MD5
  echo $MD5fc
  if [ $MD5 == $MD5fc ]
  then
    error=0
  else
    error=1
  fi
  rm -f $yuv
  if [ $error -ne 0 ]
  then
    exit $error
  fi
done

exit 0
