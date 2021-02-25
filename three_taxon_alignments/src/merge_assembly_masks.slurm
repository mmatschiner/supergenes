#!/bin/bash

# m_matschiner Sat Jan 6 16:28:54 CET 2018

# Job name:
#SBATCH --job-name=wgamrg
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=2G

## Set up the job environment.
source /cluster/bin/jobsetup
module load ruby/2.1.5

## Feedback.
echo
echo "merge_assembly_masks.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
rm_file=${1}

# Merge the versions of the assemblies that are masked with either repeatmasker or trf.
ruby merge_assembly_masks.rb ${rm_file}

