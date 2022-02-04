#!/bin/bash -x

dir=$(pwd)

{{#wilson_flows}}
store_dir="${dir}_{{wtime}}_{{wsteps}}"
test -d ${store_dir} || mkdir -p ${store_dir}
for file in ${store_dir}/*
do
  unlink ${file}
done
for file in ${dir}/output/*{{> gauge_id}}.lime
do
  filename=${file##*/}
  ln -s ${file} ${store_dir}/${filename%%_0_*}.lime
done

{{/wilson_flows}}
