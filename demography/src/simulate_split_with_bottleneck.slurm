#!/bin/bash

# michaelm Thu Jul 9 15:39:50 CEST 2020

# Job name:
#SBATCH --job-name=msprime
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=1

## Set up job environment:
## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load Python/3.8.2-GCCcore-9.3.0
module load GSL/2.6-GCC-9.3.0
module list

## Feedback.
echo
echo "simulate_split_with_bottleneck.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
recombination_rate=${1}
mutation_rate_per_yr=${2}
generation_time=${3}
chr_length=${4}
res=${5}

# Run simulations with msprime
echo -e "bottleneck_time\tpop_size\tn_var1\tn_var2\tpi1\tpi2" > ${res}
for bottleneck_time in 1000 3000 10000 30000 100000 300000 1000000 3000000
do
    for population_size in 1000 3000 10000 30000 100000
    do
        for rep in {1..20}
        do
            python simulate_split_with_bottleneck.py ${recombination_rate} ${mutation_rate_per_yr} ${generation_time} ${bottleneck_time} ${population_size} ${chr_length} >> ${res}
        done
    done
done
