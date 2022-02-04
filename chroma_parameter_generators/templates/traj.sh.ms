#!/bin/bash -x
{{#pbs}}
#PBS -N {{job_name}}
#PBS -l walltime={{time}}
#PBS -l mem={{cpus}}GB
#PBS -l jobfs=5GB
#PBS -l ncpus={{cpus}}
#PBS -q {{partition}}
{{#project}}#PBS -P {{.}}{{/project}}
#PBS -l wd
{{/pbs}}
{{#slurm}}
#SBATCH --job-name=gather_{{job_name}}
#SBATCH --nodes=1
#SBATCH --ntasks={{cpus}}
#SBATCH --ntasks-per-node={{cpus}}
#SBATCH --time={{time}}
#SBATCH --partition={{partition}}
{{/slurm}}


{{#raijin}}
module purge
module load gcc/5.2.0
module load pbs
module load openmpi/1.6.3
module load python/2.7.6
module load cuda/6.5
{{/raijin}}

set -e 
WILSON_EXE="chroma_wilson_shift"
dir=$(pwd)
shift_dir="${dir}/shift"

test -d "${shift_dir}" || mkdir -p "${shift_dir}"
cd ${shift_dir}

{{#wilson_flows}}
rm trajrun* || :
xargs -n 85 -P {{cpus}} "${SRC_DIR}/chroma/wilson_plaq_runner.sh" "chroma_clover_gmunu" "{{beta}}" <<< $(ls ${dir}_{{wtime}}_{{wsteps}}/*.lime)
python "${SRC_DIR}/chroma_parameter_generators/merge_plaq_results.py" as_{{wtime}}_{{wsteps}}.yml trajrun*.yml
#mean_plaq=$(hl_json_extract.py as_{{wtime}}_{{wsteps}}.yml "mean_plaq")
#echo $mean_plaq
{{/wilson_flows}}

rm trajrun* || :
