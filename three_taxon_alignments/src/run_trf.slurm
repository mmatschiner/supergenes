#!/bin/bash

# m_matschiner Wed Apr 4 13:18:02 CEST 2018

# Job name:
#SBATCH --job-name=wgatrf
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G

## Set up the job environment.
source /cluster/bin/jobsetup
module load trf/4.07b

## Feedback.
echo
echo "run_trf.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
in_file_with_relative_path=${1}
in_file_with_absolute_path=`readlink -f ${in_file_with_relative_path}`
in_file=`basename ${in_file_with_relative_path}`
out_file_with_relative_path=${2}
out_file_with_absolute_path=`readlink -f ${out_file_with_relative_path}`

# Go to the temporary directory.
cd ${SCRATCH}

# Run trf.
trf ${in_file_with_absolute_path} 2 7 7 80 10 50 2000 -m -h

# Copy the output file back to the data directory.
cp ${in_file}.2.7.7.80.10.50.2000.mask ${out_file_with_absolute_path}

