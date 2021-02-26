#!/bin/bash

# michaelm Mon Aug 17 16:17:25 CEST 2020

# Job name:
#SBATCH --job-name=adjanc
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=1

## Set up job environment:
## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BLAST+/2.9.0-iimpi-2019b
module load Ruby/2.7.1-GCCcore-8.3.0
module list

## Feedback.
echo
echo "adjust_ancestral.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
gzvcf=`readlink -f ${1}`
adjusted_gzvcf=`readlink -f ${2}`
reference_fasta=`readlink -f ${3}`
ancestral_fasta=`readlink -f ${4}`

# Copy the ruby script to the temporary directory.
cp adjust_ancestral.rb ${SCRATCH}

# Change to the temporary directory.
cd ${SCRATCH}

# Uncompress the gzvcf.
gunzip -c ${gzvcf} > tmp.vcf

# Run ruby script to adjust reference alleles according to ancestral sequence.
ruby adjust_ancestral.rb tmp.vcf tmp.adjusted.vcf ${reference_fasta} ${ancestral_fasta}

# Compress the adjusted vcf.
gzip -c tmp.adjusted.vcf > ${adjusted_gzvcf}

