#!/bin/bash

# m_matschiner Thu Jun 18 12:04:02 CEST 2020

# Job name:
#SBATCH --job-name=subsamp
#
# Wall clock limit:
#SBATCH --time=2:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=1

## Set up job environment:
## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load SAMtools/1.10-GCC-8.3.0
module list

## Feedback.
echo
echo "subsample.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
in_bam=${1}
out_bam=${2}
sampling_proportion=${3}

# Subsample the bam file with samtools.
samtools view -s ${sampling_proportion} -b -o ${out_bam} ${in_bam}
samtools index ${out_bam}
