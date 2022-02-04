#!/bin/bash -x

set -e
files=${@:2}
beta=$1

rm *.yml | :
#echo $files
xargs -n 1 -P 2 "${SRC_DIR}/chroma_parameter_generators/wilson_plaw_runner.sh" "${TODO}/chroma_wilson_shift" "$beta" <<< $files
#xargs -n 100 -P 16 "${SRC_DIR}/chroma_parameter_generators/wilson_plaw_runner.sh" "$beta" <<< $files
python "${SRC_DIR}/chroma_parameter_generators/merge_plaq_results.py" as.yml trajrun*.yml
