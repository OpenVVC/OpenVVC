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

ext_list="266
          bin"

append(){
    var=${1}
    shift
    eval "${var}=\"\$${var} ${*}\""
}

mkdir -p $STREAM
if [[ "$#" == 3 ]];  then
  STREAMLIST=$(curl --silent $URL | grep -o -E '"([[:alnum:]]+_)+([[:alnum:]]+.(266|bin))"' | sed 's/\"/\ /g')
fi

for file in $STREAMLIST
do
  if [[ ! -f "$STREAM/$file" ]]; then
    echo Downloading $file ...
    curl --silent $URL/$file --output $STREAM/$file
    echo Downloading ${file%.266}.md5 ...
    curl --silent $URL/${file%.(266|bin)}.md5 --output $STREAM/${file%.(266|bin)}.md5
  fi
done

# Construct list of files based on extension rules
for ext in ${ext_list}; do
  append file_list $(find ${STREAM} -name "*.${ext}")
done

rm -f failed.txt
for file in ${file_list}
do
  name=$(basename ${file} | sed -e "s/\.bin$//g" | sed -e "s/\.266$//g")
  src_dir=$(dirname ${file})
  yuv_file="${src_dir}/${name}.yuv"
  log_file="${src_dir}/${name}.log"
  md5_file="${src_dir}/${name}.md5"
  $DECODER -i "${file}" -o ${yuv_file} 2> ${log_file}
  MD5=$(md5sum ${yuv_file} | grep -o '[0-9,a-f]*\ ')
  MD5fc=$(cat ${md5_file} | grep -o '[0-9,a-f]*\ ')
  if [[ $MD5 == $MD5fc ]]; then
    echo -e $GREEN${name}
    echo -e Computed MD5:'\t'$MD5 $NC
    rm -f ${log_file}.log
  else
    echo ${name} >> failed.txt
    echo -e $RED${name}
    echo -e Computed MD5:'\t'$MD5
    echo -e Reference MD5:'\t'$MD5fc $NC
    cat ${log_file}
    ((error=error+1))
  fi
  echo
  rm -f ${yuv_file}
done

exit $error
