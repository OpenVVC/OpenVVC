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
  [ -z $no_progress ] && printf "\r%-70.70s ${color}%+10s${NC}\n" "${name}" "${status}"
}

decode(){
  dec_arg="-i ${1} -o ${2} -t 8"
  $($($DECODER ${dec_arg} 2> ${3}) 2> /dev/null)
  return $?
}

log_success(){
    echo ${name} >> success.txt
    append pass_list ${name}
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
  dump_md5error >> ${ERROR_LOG_FILE}
}

on_missing_md5_file(){
    increment nb_error
    log_status "NO MD5 REF" $RED
    log_error "${name}"
    log_error "Could not find a md5 reference."
    append nmd5_list ${name}
    return 1
}

find_md5_file() {
  src_dir=$(dirname ${1})
  md5_file="${src_dir}/${name}.md5"
  [ -f $md5_file ] || on_missing_md5_file
}

check_md5sum(){
  out_md5=$(md5sum ${yuv_file} | cut -f 1 -d ' ')
  ref_md5=$(cat    ${md5_file} | cut -f 1 -d ' ' | tr ABCDEF abcdef | sed 's/[^a-f0-9]//g')
  test "${out_md5}" = "${ref_md5}" || on_md5sum_mismatch
}

increment(){
    eval "${1}=$((${1} + 1))"
}

on_decoding_error(){
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
    append segf_list ${name}

    error="Decoder issue"
}

on_md5sum_mismatch(){
    increment nb_error
    log_failure
    error="MD5 mismatch"
    append fail_list ${name}
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
  append file_list $(find ${STREAM} -maxdepth 1 -name "*.${ext}" | sort)
done

short_name(){
    printf "%.$2s" "$1"
}

print_summary(){
    printf "$fail_list \n" | sed 's/ /\n/g' | column -c $(tput cols)
}

gen_file(){
lfile="$1"
shift
for v ; do
    printf "$v\n" >> ${lfile}
done
}

log_cleanup() {
	rm $tmp_fail;
	rm $tmp_pass;
	rm $tmp_segf;
	rm $tmp_nmd5;
	rm -rf $tmp_dir;
}

on_interrupt() {
	printf "\r"
	tput el
	log_cleanup ;
	exit 0;
}

trap on_interrupt INT

#TODO remove this one
rm -f failed.txt
rm -f ${ERROR_LOG_FILE}

nb_files="$(echo "$file_list" | wc -w)"
tmp_dir=$(mktemp -d -t ovnreg_XXXX)
tmp_fail=$(mktemp -p . -t .ovnreg_failXXX)
tmp_pass=$(mktemp -p . -t .ovnreg_passXXX)
tmp_segf=$(mktemp -p . -t .ovnreg_segfXXX)
tmp_nmd5=$(mktemp -p . -t .ovnreg_nmd5XXX)

nb_error=0
file_id=0
no_progress=''
keep_log=''

for file in ${file_list}; do

  name=$(basename ${file})

  increment file_id

  [ -z $no_progress ] && printf "%-62.62s %s" "Processing $(short_name "$name" 30)..." "$file_id of $nb_files"

  filter_extension name ${ext_list}

  yuv_file="${tmp_dir}/${name}.yuv"
  log_file="${tmp_dir}/${name}.log"

  find_md5_file ${file} || continue

  decode ${file} ${yuv_file} ${log_file} || on_decoding_error

  has_error && cleanup_onfail && continue

  check_md5sum ${yuv_file} ${file}

  has_error && cleanup_onfail && continue

  log_success
  cleanup_onsuccess

done

gen_file ${tmp_fail} ${fail_list}
gen_file ${tmp_pass} ${pass_list}
gen_file ${tmp_segf} ${segf_list}
gen_file ${tmp_nmd5} ${nmd5_list}

append fail_list ${nmd5_list} ${segf_list}

printf "${RED}Detected ${nb_error} errors in ${nb_files} files$NC: See $ERROR_LOG_FILE for more info.\n"

[ -z $summary ] && print_summary

[ -z $keep_log ] && rm -r $tmp_dir

exit $nb_error
