#!/usr/bin/bash -x
#set -o errexit
set -o pipefail

command_buffer=()
flush_command_buffer(){
  if [[ "${#command_buffer[@]}" != 0 ]] ; then
    local com
    local i
    for ((i=${#command_buffer[@]}-1;i>=0;i--)) ; do
      echo "CBUFFER -> ${command_buffer[$i]}"
      eval "${command_buffer[$i]}"
     done
  fi
  command_buffer=()
}
trap flush_command_buffer EXIT
push_command(){
  command_buffer[${#command_buffer[@]}]="${1}"
}

mv_line(){
  in_file=$1
  out_file=$2
  line=$3
  
  echo "REMOVING ${line} from ${in_file}"
  sed -i "/.*${line}.*/d" "${in_file}"
  echo "${line}" >> ${out_file}
}
