#!/bin/bash

# m_matschiner Wed Jun 27 13:46:17 CEST 2018

# Job name:
#SBATCH --job-name=thread
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=60G

## Set up the job environment.
source /cluster/bin/jobsetup
module load ruby/2.1.5

## Feedback.
echo
echo "thread_alignment.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
merged_alignment=${1}
threaded_alignment=${2}

# Thread the alignment according to the reference sequences.
ruby thread_alignment.rb ${merged_alignment} ${threaded_alignment}


