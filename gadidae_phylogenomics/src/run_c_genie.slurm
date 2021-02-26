#!/bin/bash

# m_matschiner Sat Sep 1 18:27:33 CEST 2018
#
# Job name:
#SBATCH --job-name=cgenie
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G

## Set up the job environment.
source /cluster/bin/jobsetup
module load python3/3.5.0
module load gsl/2.2

## Feedback.
echo
echo "run_c_genie.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
tree=${1}
shift
res=${1}
shift
gen_time=${1}
shift
rec_map=${1}
shift

# Run c-genie.
if [ -z ${1} ]
then
    python3 c-genie ${tree} ${res} -g ${gen_time} -m ${rec_map}
else
#    for file in "$@"; do
#	fastqs_with_relative_path+=(${file})
#    done
    python3 c-genie ${tree} ${res} -g ${gen_time} -m ${rec_map} "$@"
fi
