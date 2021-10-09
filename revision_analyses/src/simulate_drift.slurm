#!/bin/bash

# michaelm Sat Oct 9 11:44:29 CEST 2021
#
# This script runs simulations of drift to determine how
# quickly either loss or fixation is reached.
#
# Start with 'sbatch simulate_drift.slurm 100000 10000 out.txt'

# Job name:
#SBATCH --job-name=simdrift
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G

# Set up job environment.
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load Ruby/2.7.2-GCCcore-9.3.0
module list

# Get command line arguments and set variables.
pop_size=${1}
n_simulations=${2}
out_file=${3}

# Simulate drift.
ruby simulate_drift.rb ${pop_size} ${n_simulations} > ${out_file}
