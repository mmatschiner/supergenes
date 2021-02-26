#!/bin/bash

# m_matschiner Fri Jun 5 21:36:44 CEST 2020

# Job name:
#SBATCH --job-name=beagle
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2

## Set up job environment:
## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load Java/1.8.0_212
module list

## Feedback.
echo
echo "run_beagle.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
gzvcf=`readlink -f ${1}`
phased_gzvcf=`readlink -f ${2}`

# Copy beagle to the temporary directory.
cp ../bin/beagle.18May20.d20.jar ${SCRATCH}

# Replace "." with "./." in the vcf.
zcat ${gzvcf} | perl -pe "s/\s\.:/\t.\/.:/g" > ${SCRATCH}/tmp.vcf

# Move to the temporary directory.
cd ${SCRATCH}

# Run beagle.
java -jar -Xmx10G ./beagle.18May20.d20.jar nthreads=2 ne=10000 gt=tmp.vcf out=${phased_gzvcf%.vcf.gz}

