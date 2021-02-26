#!/bin/bash

# m_matschiner Wed Sep 26 17:34:43 CEST 2018
#
# Job name:
#SBATCH --job-name=samclp
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G

## Set up the job environment.
source /cluster/bin/jobsetup
module load samtools/1.4
module load samclip/v0.2

## Feedback.
echo
echo "remove_softclipped_reads.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
bam=${1}
ref=${2}
bam_red=${3}

# Remove soft-clipped reads.
samtools view -h ${bam} | samclip --ref ${ref} | samtools sort > ${bam_red}
samtools index ${bam_red}