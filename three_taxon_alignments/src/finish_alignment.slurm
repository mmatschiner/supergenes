#!/bin/bash

# m_matschiner Fri Jun 29 23:38:42 CEST 2018

# Job name:
#SBATCH --job-name=fnshal
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G

## Set up the job environment.
source /cluster/bin/jobsetup
module load ruby/2.1.5

## Feedback.
echo
echo "finish_alignment.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
refined_alignment=${1}
mapping_alignment=${2}
finished_alignment=${3}

# Finish the alignment by removing sites that open gaps in the reference, and by including only one sequence for gadMor_Stat and melAeg, only in regions where at least 100 bp are exactly identical in all three sequences representing this species. Also exclude sites that differ in the alignment based on mapping.
ruby finish_alignment.rb ${refined_alignment} ${mapping_alignment} ${finished_alignment}
