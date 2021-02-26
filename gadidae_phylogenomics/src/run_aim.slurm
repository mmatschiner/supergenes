#!/bin/bash

# m_matschiner Tue Oct 16 01:04:58 CEST 2018

# Job name:
#SBATCH --job-name=aim
#
# Wall clock limit:
#SBATCH --time=168:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#
# Outfile:
#SBATCH --output=aim.out

## Set up the job environment.
source /cluster/bin/jobsetup

# Run aim.
n_logs=`ls *.log | wc -l`
if (( ${n_logs} > 0 ))
then
    beast -resume -seed ${RANDOM} aim.xml
else
    beast -seed ${RANDOM} aim.xml
fi
