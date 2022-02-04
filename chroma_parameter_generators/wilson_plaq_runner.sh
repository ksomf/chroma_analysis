#!/bin/bash -x

set -e

dest=$(pwd)
dir=$(mktemp -d)
cd $dir

tag=${dir##*.}
$@
mv "${dir}/action_shift.yml" "${dest}/trajrun${tag}_as.yml"
mv "${dir}/action_shift.dat" "${dest}/trajrun${tag}_as.dat"

cd $dest
rm -r $dir
