#!/bin/bash -x

source ../utils/compiler_flags || return 1

COMPILER_FLAGS="${COMPILER_FLAGS} -DhlSSE -DhlSIMD"
COMPILER_FLAGS="${COMPILER_FLAGS} -DhlDEBUG "
ZERO_COMPILER_FLAGS="${COMPILER_FLAGS} -DhlZERO_MOM"
ONE_COMPILER_FLAGS="${COMPILER_FLAGS} -DhlONE_MOM"
SIMPLE_COMPILER_FLAGS="${COMPILER_FLAGS} -DhlSIMPLE_MOM"
WILSON_ACTION_CCF="${COMPILER_FLAGS}"
CLOVER_GMUNU_CCF="${COMPILER_FLAGS} -DCLOVER_GMUNU" # -DhlAVX512"

BUILD_DIR='../bin'
SRC_DIR=$(pwd)

test -d ${BUILD_DIR} || mkdir -p ${BUILD_DIR}
pushd ${BUILD_DIR} > /dev/null

return_result=0
${COMPILER} ${COMPILER_FLAGS}          -x c ${SRC_DIR}/chroma_extract.c -o chroma_extract            ${LINKER_FLAGS} || return_result=1
${COMPILER} ${ZERO_COMPILER_FLAGS}     -x c ${SRC_DIR}/chroma_extract.c -o chroma_extract_zero       ${LINKER_FLAGS} || return_result=1
${COMPILER} ${ONE_COMPILER_FLAGS}      -x c ${SRC_DIR}/chroma_extract.c -o chroma_extract_one        ${LINKER_FLAGS} || return_result=1
${COMPILER} ${SIMPLE_COMPILER_FLAGS}   -x c ${SRC_DIR}/chroma_extract.c -o chroma_extract_simple     ${LINKER_FLAGS} || return_result=1

${COMPILER} ${WILSON_ACTION_CCF}  -x c ${SRC_DIR}/plaq_shift.c   -o chroma_wilson_shift ${LINKER_FLAGS} || return_result=1
${COMPILER} ${CLOVER_GMUNU_CCF}   -x c ${SRC_DIR}/plaq_shift.c   -o chroma_clover_gmunu ${LINKER_FLAGS} || return_result=1

popd > /dev/null
