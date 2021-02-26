#!/bin/bash

# m_matschiner Fri Jul 5 22:54:11 CEST 2019

# Job name:
#SBATCH --job-name=snapp
#
# Project:
# Wall clock limit:
#
#SBATCH --time=100:00:00
#
# Processor and memory usage:
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#
# Outfile:
#SBATCH --output=QQQQQQ.out

## Set up job environment:
## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load Beast/2.5.2-GCC-8.2.0-2.31.1
module list

## Run or resume beast analysis.
if [ ! -f QQQQQQ.log ]
then
    beast -threads 8 -seed ${RANDOM} QQQQQQ.xml
else
    beast -threads 8 -seed ${RANDOM} -resume QQQQQQ.xml
fi
