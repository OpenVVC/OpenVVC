#!/bin/bash

print_usage(){
    cat <<EOF
USAGE ./checkMD5.sh <bitstreams folder> <decoder> [<url>]
EOF
    exit 0
}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# Default values
ext_list="266
          bin
          h266
          vvc"


log_error(){
    echo -e "$RED${*}$NC"
}

die(){
    echo -e "$RED${*}$NC"
    exit 1
}

if [[ "$#" < 2 ]];  then
    log_error "Illegal number of parameters"
    print_usage
    die
fi

STREAM=$1
DECODER=$2
URL=$3

append(){
    var=${1}
    shift
    eval "${var}=\"\$${var} ${*}\""
}

filter_extension(){
  var=$1
  shift
  for ext in $*; do
    eval "var2=$(echo \$${var})"
    tmp=$(echo ${var2} | sed -e "s/\(.*\)\(\.${ext}\$\)/\1/g")
    eval "${var}=${tmp}"
  done
  unset tmp
  unset var
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

decode(){
  #TODO handle /dev/null output and optional log
  dec_arg="-i ${1} -o ${2}"
  $DECODER ${dec_arg} 2> ${3}
  return $?
}

log_success(){
    echo -e $GREEN${name}
    echo -e Computed MD5:'\t'${out_md5} $NC
    rm -f ${log_file}.log
    echo
}

log_failure(){
    echo ${name} >> failed.txt
    echo -e $RED${name}
    echo -e Computed  MD5:'\t'${out_md5} $NC
    echo -e Reference MD5:'\t'${ref_md5} $NC
    echo
    cat ${log_file}
}

check_md5sum(){
  src_dir=$(dirname ${file})
  md5_file="${src_dir}/${name}.md5"

  out_md5=$(md5sum ${yuv_file} | grep -o '[0-9,a-f]*\ ')
  ref_md5=$(cat    ${md5_file} | grep -o '[0-9,a-f]*\ ')

  test "${out_md5}" = "${ref_md5}"
  return $?
}

increment(){
    eval "((${1}=${1}+1))"
}

handle_decoding_error(){
    increment nb_error
    log_error "Error while decoding ${file}"
    cat ${log_file}
    error="Decoder issue"
}

handle_md5sum_mismatch(){
    increment nb_error && log_failure
    error="MD5 mismatch"
}

has_error(){
   test -n "$error"
   return $?
}

cleanup_onfail(){
  unset error
  rm -f ${yuv_file}
}

cleanup_onsuccess(){
  unset error
  rm -f ${log_file}
  rm -f ${yuv_file}
}

# Construct list of files based on extension rules
for ext in ${ext_list}; do
  append file_list $(find ${STREAM} -name "*.${ext}")
done

rm -f failed.txt

tmp_dir=$STREAM
nb_error=0

for file in ${file_list}; do

  name=$(basename ${file})

  filter_extension name ${ext_list}

  yuv_file="${tmp_dir}/${name}.yuv"
  log_file="${tmp_dir}/${name}.log"

  decode ${file} ${yuv_file} ${log_file} || handle_decoding_error

  has_error && cleanup_onfail && continue

  check_md5sum ${yuv_file} ${file} && log_success || handle_md5sum_mismatch

  has_error && cleanup_onfail && continue

  cleanup_onsuccess

done

exit $nb_error
