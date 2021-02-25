#!/bin/bash

# m_matschiner Fri Aug 17 11:40:22 CEST 2018

# Job name:
#SBATCH --job-name=beast2
#
# Wall clock limit:
#SBATCH --time=168:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#
# Outfile:
#SBATCH --output=91genes.out

## Set up the job environment.
source /cluster/bin/jobsetup
module load beast2/2.5.0

# Run beast.
beast -seed ${RANDOM} -beagle 91genes.xml