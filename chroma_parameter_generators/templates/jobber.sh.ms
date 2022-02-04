#!/bin/bash -x
{{#pbs}}
#PBS -l walltime={{walltime}}
#PBS -q {{queue}}
#PBS -l {{^usencpus}}nodes={{nodes}}:ppn={{cpn}}{{/usencpus}}{{#usencpus}}ncpus={{cpus}}{{/usencpus}}
{{#gpus}}#PBS -l ngpus={{gpus}}{{/gpus}}
#PBS -N {{job_name}}
{{#jobfs}}#PBS -l jobfs={{jobfs}}GB{{/jobfs}}
{{#memory}}#PBS -l mem={{memory}}GB{{/memory}}
{{#project}}#PBS -P {{project}}{{/project}}
{{#account}}#PBS -A {{account}}{{/account}}
{{#feature}}#PBS -l feature={{feature}}{{/feature}}
#PBS -j oe
{{/pbs}}
{{#slurm}}
#SBATCH --job-name={{job_name}}
#SBATCH --nodes={{nodes}}
#SBATCH --time={{walltime}}
#SBATCH --output=slurm_{{name}}_{{job_name}}_%j.out
#SBATCH --ntasks={{cpus}}
#SBATCH --partition={{partition}}
#SBATCH --gres={{gres}}
#SBATCH --qos={{qos}}
#SBATCH --mem-per-cpu={{mpc}}GB
{{/slurm}}

repeats={{repeats}}

{{#pbs}}
cd $PBS_O_WORKDIR
{{/pbs}}

echo "STARTING JOB ${SLURM_JOB_NAME}${PBS_JOBNAME} : {{job_name}}"
echo "CONFIGURATION {{nodes}}x({{cpn}}{{#gpus}},{{gpus}}{{/gpus}})"

for ((i=1; i<=repeats; ++i))
do
  ./chroma_{{name}}
  
   if [[ ! "$?" == 0 ]]
   then
     echo "FINISHED RUNS"
     exit 1
   fi
done

completed_seconds=$SECONDS
time_minutes=$(($completed_seconds / 60))
echo "Job completed in ${time_minutes} minutes."

{{#resubmit}}
{{#pbs}}qsub{{/pbs}}{{#slurm}}sbatch{{/slurm}} jobber_{{name}}
{{/resubmit}}
