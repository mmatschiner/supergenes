#!/bin/bash

# michaelm Thu Sep 9 16:45:59 CEST 2021
#
# This script calculates basic statistics for chromosome windows.

# Start with
# 'sbatch get_window_stats.slurm in.vcf.gz res_dir window_size'

# Job name:
#SBATCH --job-name=getwst
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1

## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.10.2-GCC-9.3.0
module load Python/3.8.2-GCCcore-9.3.0
module list

## Feedback.
echo
echo "get_window_stats.slurm was called with arguments ${@}."
echo

## Get command-line arguments.
gzvcf=`readlink -f ${1}`
lg_id=${2}
pop_group1=${3}
pop_group2=${4}
res=`readlink -f ${5}`
window_size=${6}

# Make a joint string for the two pop gruops.
joint_group="${pop_group1},${pop_group2}"

# Make a temporary reduced version of the vcf file excluding all samples not in the joint group.
bcftools view --types snps -s ${joint_group} ${gzvcf} | bcftools filter -e 'AC==0 || AC==AN' | grep -v GmG20150304 > tmp.${lg_id}.vcf

# Split the header and the body of the vcf file.
cat tmp.${lg_id}.vcf | grep "#" > tmp.${lg_id}.header.vcf
cat tmp.${lg_id}.vcf | grep -v "#" > tmp.${lg_id}.body.vcf

# Add a pop_group id to each sample name of pop_group1.
IFS=',' read -r -a array1 <<< "${pop_group1}"
for sample in "${array1[@]}"
do
    sed -i "s/${sample}/pop1_${sample}/g" tmp.${lg_id}.header.vcf
done
IFS=',' read -r -a array2 <<< "${pop_group2}"
for sample in "${array2[@]}"
do
    sed -i "s/${sample}/pop2_${sample}/g" tmp.${lg_id}.header.vcf
done

# Combine again the header and the body of the vcf.
cat tmp.${lg_id}.header.vcf tmp.${lg_id}.body.vcf > tmp.${lg_id}.vcf

# Clean up.
rm -f tmp.${lg_id}.*.vcf

# Calculate population statistics from the vcf.
python run_sliding_windows_over_vcf.py tmp.${lg_id}.vcf ${res} -w ${window_size} -c ${lg_id} -p pop1 pop2

# Clean up.
rm -f tmp.${lg_id}.vcf
