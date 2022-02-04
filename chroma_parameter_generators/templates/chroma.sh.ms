#!/bin/bash -x

source hl.sh || exit 1

run_dir=$(pwd)
filelist="${run_dir}/filelist"
inproglist="${run_dir}/inproglist"
donelist="${run_dir}/donelist"
errorlist="${run_dir}/errorlist"
readylist="${run_dir}/readylist"
output_dir="${run_dir}/output"
fail_dir="${run_dir}/failed"
lock_dir="${run_dir}/lock.tmp"
{{#jobfs}}tmp_dir="${PBS_JOBFS}/tmp${PBS_JOBID}_$RANDOM"{{/jobfs}}
{{#jureca}}tmp_dir="$WORK/tmp${SLURM_JOB_ID}_$RANDOM"{{/jureca}}
{{#isaac}}tmp_dir="$WORK/tmp/${PBS_JOBID}_$RANDOM"{{/isaac}}
{{#phoenix}}tmp_dir="$TMPFS/$HOSTNAME/${SLURM_JOB_ID}${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}_$RANDOM"{{/phoenix}}
{{#hlrn}}tmp_dir="${run_dir}/tmp/${PBS_JOBID}_$RANDOM"{{/hlrn}}
{{#cpu_contraction}}
tmp_dir="$1"
tag="$2"
{{/cpu_contraction}}


test -f ${filelist}   || (echo "SETUP FILELIST NOT FOUND"; exit 1)
test -f ${inproglist} || touch ${inproglist}
test -f ${donelist}   || touch ${donelist}
test -f ${errorlist}  || touch ${errorlist}
test -f ${readylist}  || touch ${readylist}
test -f ${lock_dir}   || touch ${lock_dir}
test -d ${output_dir} || mkdir --parents "${output_dir}"
test -d ${fail_dir}   || mkdir --parents "${fail_dir}"
test -d ${tmp_dir}   || mkdir --parents "${tmp_dir}"

set -eu
{{#cpu_contraction}}
if [ -z "$3" ]; then exit 1; fi
mv "$3"/* "${tmp_dir}"
push_command 'rmdir "$3"'
{{/cpu_contraction}}
push_command 'rm -r "${tmp_dir}"'

{{#raijin}}
module purge
module load gcc/4.9.0 pbs openmpi/2.1.2 python/2.7.6 cuda/9.0
export CUDA_HOME="${CUDA_BASE}"
ENVIRON="-envvar OMPI_COMM_WORLD_LOCAL_RANK"
CALLER="mpiexec -np {{run_cpus}}"
CPU_EXE="${CHROMA_INSTALL}/ac_chroma_cpu_node_{{cpu_precision}}/bin/chroma"
GPU_EXE="${CHROMA_INSTALL}/jz_chroma_gpu_node_{{gpu_precision}}/bin/chroma"
GPU_CALLER="mpiexec -np {{gpus}}"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${CHROMA_INSTALL}/ac_chroma_cpu_node_{{cpu_precision}}:$CUDA_HOME/lib64:$CUDA_HOME/nvvm/lib64"
#{{#knl_cpus}}
#module purge
#module load intel-cc/16.0.3.210 pbs intel-mpi/5.1.3.210 python/2.7.6 intel-mkl/16.0.3.210 intel-fc/16.0.3.210
#export OMP_NUM_THREADS={{knl_logical_cpus}}
#CALLER="mpirun -n {{knl_logical_cpus}} -ppn {{knl_cpn}} "
#CPU_EXE="${CHROMA_INSTALL}/acknl_chroma_knl_node_double/bin/chroma"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${CHROMA_INSTALL}/acknl_chroma_knl_node_double/lib:$CUDA_HOME/lib64:$CUDA_HOME/nvvm/lib64"
#{{/knl_cpus}}
{{/raijin}}
{{#hlrn}}
module switch PrgEnv-cray PrgEnv-gnu
module load python
CPU_EXE="${CHROMA_INSTALL}/install-chroma-double-AlexGit-gnu5/chroma_qcdsf/bin/chroma"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CHROMA_INSTALL}/install-chroma-double-AlexGit-gnu5/chroma_qcdsf/lib
CALLER="aprun -n {{run_cpus}} -N {{cpn}} ${CPU_EXE} -geom {{geometry}}"
{{/hlrn}}
{{#phoenix}}
module purge
module load GCC/4.9.3-binutils-2.25 OpenMPI/1.6.5-GNU-4.9.3-2.25 CUDA/6.5.14 Autoconf/2.69-GNU-4.9.3-2.25
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${CHROMA_INSTALL}/ac_chroma_cpu_node_double:$CUDA_HOME/lib64:$CUDA_HOME/nvvm/lib64"
ENVIRON="-envvar OMPI_COMM_WORLD_LOCAL_RANK"
CPU_EXE="${CHROMA_INSTALL}/ac_chroma_cpu_node_double/bin/chroma"
GPU_EXE="${CHROMA_INSTALL}/jz_chroma_gpu_node_single/bin/chroma"
GPU_CALLER="mpiexec -np 4"
CALLER="mpiexec -np 16"
{{/phoenix}}



{{^cpu_contraction}}
{
  flock -w 10 -x 200 
  tag=$({{#srcdir}}python {{.}}/chroma_parameter_generators/{{/srcdir}}generate_template_for_next_trajectory.py "${tmp_dir}" "{{sources}}" "{{name}}" "{{geometry}}{{^geometry}}no geometry set{{/geometry}}") || exit 1
} 200>${lock_dir}
{{/cpu_contraction}}

echo "RUNNING WITH ${tag}"
{{#gpu_inversion}}
intermediate_output_dir="${run_dir}/tmp/${tag}"
test -d "${intermediate_output_dir}" || mkdir --parents "${intermediate_output_dir}"
{{/gpu_inversion}}

  
result_files=()
input_files=()
output_files=()
{{^unified}}
for template in input_{{name}}*.ms
do
    # Get head of input template and create filenames of input and output files
    input_name="${template%%.*}"
    input_filename="${tmp_dir}/${input_name}_${tag}.xml"
    output_filename="${tmp_dir}/output${input_name##input}_${tag}.xml"
    input_files+=("${input_filename}")
    output_files+=("${output_filename}")
    {{^gpu_inversion}}
    result_files+=("${input_filename}")
    result_files+=("${output_filename}")
    {{/gpu_inversion}}
done
{{/unified}}
{{#unified}}
template="input_{{name}}_unified.xml.ms"
input_name="${template%%.*}"
input_filename="${tmp_dir}/${input_name}_${tag}.xml"
output_filename="${tmp_dir}/output${input_name##output}_${tag}.xml"
input_files+=("${input_filename}")
output_files+=("${output_filename}")
result_files+=("${input_filename}")
result_files+=("${output_filename}")
{{/unified}}

{{^gpu_inversion}}
{{#slices}}
result_files+=("${tmp_dir}/baryon_${tag}_{{#slice}}{{.}}{{/slice}}{{^slice}}0.000{{/slice}}.xml")
result_files+=("${tmp_dir}/baryon_${tag}_{{#slice}}{{.}}{{/slice}}{{^slice}}0.000{{/slice}}.lime")
result_files+=("${tmp_dir}/meson_${tag}_{{#slice}}{{.}}{{/slice}}{{^slice}}0.000{{/slice}}.xml")
result_files+=("${tmp_dir}/meson_${tag}_{{#slice}}{{.}}{{/slice}}{{^slice}}0.000{{/slice}}.lime")
{{/slices}}
{{^slices}}
result_files+=("${tmp_dir}/baryon_${tag}_0.000.xml")
result_files+=("${tmp_dir}/baryon_${tag}_0.000.lime")
{{/slices}}
{{#wilson_flows}}
result_files+=("${tmp_dir}/${tag}_{{> gauge_id}}.lime")
{{/wilson_flows}}

result_destination="${fail_dir}/"
{{/gpu_inversion}}
after_run_list="${errorlist}"

{{^cpu_contraction}}
push_command 'mv_line "${inproglist}" "${after_run_list}" "${tag}"'
{{/cpu_contraction}}
{{#cpu_contraction}}
push_command 'mv_line "${readylist}" "${after_run_list}" "${tag}"'
{{/cpu_contraction}}
{{^gpu_inversion}}
push_command 'mv "${result_files[@]}" "${result_destination}" || :'
{{/gpu_inversion}}

{{^unified}}
chrome_exit_status=0
echo "STARTING UNIFIED RUN WITH {{nodes}}x({{cpn}}{{#gpus}},{{gpus}}{{/gpus}})"
{{^cpu_contraction}}
{{#gpus}}
${GPU_CALLER} "${GPU_EXE}" ${ENVIRON} -geom 1 1 1 {{.}} -i "${input_files[0]}" -o "${output_files[0]}"
{{/gpus}}
{{^gpus}}
${CALLER} "${CPU_EXE}" -geom {{geometry}} -i "${input_files[0]}" -o "${output_files[0]}"
{{/gpus}}
chroma_exit_status="$?"
{{/cpu_contraction}}
{{^gpu_inversion}}{{#combine_runs}}
if [[ "${chroma_exit_status}" == 0 ]]
then
  echo "STARTING CHROMA COMBINE RUN {{input_file_index}}"
  ${CALLER} "${CPU_EXE}" -i "${input_files[{{input_file_index}}]}" -o "${output_files[{{input_file_index}}]}"
  chroma_exit_status="$?"
fi
{{/combine_runs}}{{/gpu_inversion}}
{{/unified}}
{{#unified}}
echo "STARTING UNIFIED RUN WITH {{nodes}}x({{cpn}}{{#gpus}},{{gpus}}{{/gpus}})"
${CALLER} "${CPU_EXE}" -geom {{geometry}} -i "${input_files[0]}" -o "${output_files[0]}"
chroma_exit_status="$?"
{{/unified}}

if [[ "${chroma_exit_status}" == 0 ]]
then
  echo "Chroma ran successfully"
  {{^gpu_inversion}}
  result_destination="${output_dir}"
  after_run_list="${donelist}"
  {{/gpu_inversion}}
  {{#gpu_inversion}}
  after_run_list="${readylist}"
  ls "${tmp_dir}"
  mv "${tmp_dir}"/* "${intermediate_output_dir}"
  sbatch jobber_cpu_contraction "${tmp_dir}" "${tag}" "${intermediate_output_dir}"
  {{/gpu_inversion}}
else
  echo "Chroma ran badly"
  exit 1
fi
