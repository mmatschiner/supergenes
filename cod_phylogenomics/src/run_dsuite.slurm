#!/bin/bash

# michaelm Mon Feb 8 00:59:50 CET 2021
#
# This script runs the dtrios function of dsuite.
#
# Start with
# 'sbatch run_dsuite.slurm data.vcf.gz tree.tre'

# Job name:
#SBATCH --job-name=dsuite
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=5G

## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.10.2-GCC-9.3.0

## Feedback.
echo
echo "run_dsuite.slurm was called with arguments ${@}."
echo

## Get command line arguments.
gzvcf=${1}
tre=${2}
res_dir=${3}

# Get the vcf id.
gzvcf_id=`basename ${gzvcf%.full.vcf.gz}`
            
# Get the list of samples.
if [ ${gzvcf_id} == "inversion_lg12" ]
then
    samples=`bcftools view -s ^Gadmor_lfc1,Gadmor_lfc2 -h ${gzvcf} | tail -n 1 | cut -f 10-`
else
    samples=`bcftools view -h ${gzvcf} | tail -n 1 | cut -f 10-`
fi

# Make the species sets file.
sets=${res_dir}/${gzvcf_id}.txt
rm -f ${sets}
touch ${sets}
for sample in ${samples}
do
	echo -ne "${sample}\t" >> ${sets}
	echo "${sample:0:10}" | sed 's/Borsai/Outgroup/g' >> ${sets}
done

# Run the dtrios function of dsuite.
../bin/Dsuite Dtrios -t ${tre} ${gzvcf} ${sets}
