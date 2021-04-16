#!/bin/bash

if [[ "$#" < 2 ]];  then
    echo "Illegal number of parameters"
    echo
    echo "USE ./checkMD5.sh <bitstreams folder> <decoder> [<url>]"
    exit
fi

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
STREAM=$1
DECODER=$2
URL=$3
error=0

mkdir -p $STREAM

if [[ "$#" == 3 ]];  then
  STREAMLIST=$(curl --silent $URL | grep -o -E '"([[:alnum:]]+_)+([[:alnum:]]+.266)"' | sed 's/\"/\ /g')
fi

for file in $STREAMLIST
do
  if [[ ! -f "$STREAM/$file" ]]; then
    echo Downloading $file ...
    curl --silent $URL/$file --output $STREAM/$file
    echo Downloading ${file%.266}.md5 ...
    curl --silent $URL/${file%.266}.md5 --output $STREAM/${file%.266}.md5
  fi
done

rm -f failed.txt
for file in $STREAM/*.266
do
  yuv=${file%.266}.yuv
  $DECODER -i $file -o $yuv 2> ${file%.266}.log
  MD5=$(md5sum $yuv | grep -o '[0-9,a-f]*\ ')
  MD5fc=$(cat ${file%.266}.md5 | grep -o '[0-9,a-f]*\ ')
  if [[ $MD5 == $MD5fc ]]; then
    echo -e $GREEN$(basename ${file%.266})
    echo -e Computed MD5:'\t'$MD5 $NC
    rm -f ${file%.266}.log
  else
    basename ${file%.266} >> failed.txt
    echo -e $RED$(basename ${file%.266})
    echo -e Computed MD5:'\t'$MD5
    echo -e Reference MD5:'\t'$MD5fc $NC
    cat ${file%.266}.log
    ((error=error+1))
  fi
  echo
  rm -f $yuv
done

exit $error
