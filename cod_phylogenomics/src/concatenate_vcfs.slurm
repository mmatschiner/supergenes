#!/bin/bash

# m_matschiner Sat Jan 6 16:28:54 CET 2018

# Job name:
#SBATCH --job-name=ccall
#
# Wall clock limit:
#SBATCH --time=20:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G

## Set up the job environment.
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-foss-2018b
module load GATK/4.1.0.0-foss-2018b-Python-3.6.6
module list


## Feedback.
echo
echo "concatenate_vcfs.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
ref=${1}
gzvcf_dir=${2}
concatenated_vcfgz=${3}

# Concatenate vcf files of all linkage groups.
gatk GatherVcfs \
    -R ${ref} \
    -I ${gzvcf_dir}/LG01.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG02.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG03.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG04.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG05.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG06.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG07.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG08.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG09.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG10.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG11.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG12.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG13.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG14.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG15.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG16.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG17.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG18.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG19.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG20.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG21.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG22.filtered.vcf.gz \
    -I ${gzvcf_dir}/LG23.filtered.vcf.gz \
    -O ${concatenated_vcfgz}

# Index the concatenate vcf.
bcftools index ${concatenated_vcfgz}
