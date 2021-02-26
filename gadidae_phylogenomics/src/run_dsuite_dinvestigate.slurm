#!/bin/bash

# m_matschiner Thu Oct 31 20:36:07 CET 2019

# Job name:
#SBATCH --job-name=dsuite
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#
## Set up the job environment.
source /cluster/bin/jobsetup

# Load modules.
module load dsuite/0.1.r5

# Get the command-line arguments.
vcfgz=${1}
samples=${2}
trios=${3}
res_dir=${4}
window_string="50,25"
res_file_ending=`echo __${window_string}.txt | tr "," "_"`

# Run dsuite's dinvestigate function.
Dsuite Dinvestigate -w ${window_string} ${vcfgz} ${samples} ${trios}

# Move the results to the result directory.
mv *${res_file_ending} ${res_dir}

# Clean up.
if [ ${trios:0:3} == "tmp" ]
then
    rm -f ${trios}
fi

