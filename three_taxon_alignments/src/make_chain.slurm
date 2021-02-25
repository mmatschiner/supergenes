#!/bin/bash

# m_matschiner Sat Jan 6 16:28:54 CET 2018

# Job name:
#SBATCH --job-name=wgamkc
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=1G

## Set up the job environment.
source /cluster/bin/jobsetup

## Feedback.
echo
echo "make_chain.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
psl_file=${1}
chain_file=${2}
target_genome_in_2bit=${3}
query_genome_in_2bit=${4}

# Chain the alignments of the psl file.
../bin/kent/axtChain -linearGap=medium -psl ${psl_file} ${target_genome_in_2bit} ${query_genome_in_2bit} ${chain_file}


