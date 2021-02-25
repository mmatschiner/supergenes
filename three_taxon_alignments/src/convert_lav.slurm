#!/bin/bash

# m_matschiner Sat Jan 6 16:28:54 CET 2018

# Job name:
#SBATCH --job-name=wgacnv
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G

## Set up the job environment.
source /cluster/bin/jobsetup
module load lastz/1.04

## Feedback.
echo
echo "convert_lav.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
lav_file=${1}
psl_file=${2}

# Convert the lav file to psl format.
../bin/kent/lavToPsl ${lav_file} ${psl_file}


