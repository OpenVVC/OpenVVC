#!/bin/sh

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

ERROR_LOG_FILE="error.log"

log_error(){
    echo "${*}" >> ${ERROR_LOG_FILE}
}

die(){
    echo -e "$RED${*}$NC"
    exit 1
}

if [ $# -lt 2 ];  then
    print_usage
    die "Illegal number of parameters"
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

if [ $# -eq 3 ];  then
  STREAMLIST=$(curl --silent $URL | grep -o -E '"([[:alnum:]]+_)+([[:alnum:]]+.266)"' | sed 's/\"/\ /g')
fi

for file in $STREAMLIST
do
  if [ ! -f "$STREAM/$file" ]; then
    echo Downloading $file ...
    curl --silent $URL/$file --output $STREAM/$file
    echo Downloading ${file%.266}.md5 ...
    curl --silent $URL/${file%.266}.md5 --output $STREAM/${file%.266}.md5
  fi
done

log_status(){
  status=$1
  color=$2
  printf "\r%-70.70s ${color}%+10s${NC}\n" "${name}" "${status}"
}

decode(){
  #TODO handle /dev/null output and optional log
  dec_arg="-i ${1} -o ${2} -t 8"
  $DECODER ${dec_arg} 2> ${3}
  return $?
}

log_success(){
    log_status 'PASS' $GREEN
    rm -f ${log_file}.log
}

dump_md5error(){
    cat <<EOF
${file}:
A md5 mismatch occured on file ${file}.
Reference MD5:	${ref_md5}
The decoder output has been saved to ${log_file}.
EOF
}

log_failure(){
    echo ${name} >> failed.txt
  log_status 'FAIL' $RED
    echo -e "$RED${name}$NC: See $ERROR_LOG_FILE for more info."
    dump_md5error >> ${ERROR_LOG_FILE}
}

check_md5sum(){
  out_md5=$(md5sum ${yuv_file} | grep -o '[0-9,a-f]*\ ')

  src_dir=$(dirname ${file})
  md5_file="${src_dir}/${name}.md5"
  if [ -f $md5_file ] ; then
      ref_md5=$(cat    ${md5_file} | grep -o '[0-9,a-f]*\ ')
          test "${out_md5}" = "${ref_md5}" || handle_md5sum_mismatch
  else
    error="No MD5 file"
    log_status "MISSING MD5" $RED
    log_error "${name}"
    log_error "Could not find a md5 reference."
    return 1
  fi
}

increment(){
    eval "${1}=$((${1} + 1))"
}

handle_decoding_error(){
    retval=$?
    increment nb_error
    case $retval in
      "139")
        er='SEGFAULT'
        ;;
      "134")
        er='STACK'
        ;;
      *)
        er='OTHER'
        ;;
    esac

    log_status "$er" $RED

    log_error "Error while decoding ${file}"

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

find_in_list(){
   exp=$1
   shift
   for v in $*; do
       eval "case '$v' in $exp) return 0;; esac"
   done
   return 1
}

# Construct list of files based on extension rules
for ext in ${ext_list}; do
  append file_list $(find ${STREAM} -name "*.${ext}" | sort)
done

#TODO remove this one
rm -f failed.txt
rm -f ${ERROR_LOG_FILE}

tmp_dir=$(mktemp -d -t ovnreg_XXXX)

nb_error=0

nb_files="$(echo "$file_list" | wc -w)"
echo ${nb_files}
file_id=0

short_name(){
    printf "%.$2s" "$1"
}

for file in ${file_list}; do

  name=$(basename ${file})

  increment file_id

  printf "%-62.62s %s" "Processing $(short_name "$name" 30)..." "$file_id of $nb_files"

  filter_extension name ${ext_list}

  yuv_file="${tmp_dir}/${name}.yuv"
  log_file="${tmp_dir}/${name}.log"

  decode ${file} ${yuv_file} ${log_file} || handle_decoding_error

  has_error && cleanup_onfail && continue

  check_md5sum ${yuv_file} ${file}

  has_error && cleanup_onfail && continue

  log_success
  cleanup_onsuccess

done

exit $nb_error
