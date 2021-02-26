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
source /cluster/bin/jobsetup

## Feedback.
echo
echo "concatenate_vcfs.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
ref=${1}
concatenated_vcfgz=${2}

# Concatenate vcf files of all linkage groups.
java -Xmx10G -cp /cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
    -R ${ref} \
    -V ../res/gatk/LG01.masked.variable.vcf.gz \
    -V ../res/gatk/LG02.masked.variable.vcf.gz \
    -V ../res/gatk/LG03.masked.variable.vcf.gz \
    -V ../res/gatk/LG04.masked.variable.vcf.gz \
    -V ../res/gatk/LG05.masked.variable.vcf.gz \
    -V ../res/gatk/LG06.masked.variable.vcf.gz \
    -V ../res/gatk/LG07.masked.variable.vcf.gz \
    -V ../res/gatk/LG08.masked.variable.vcf.gz \
    -V ../res/gatk/LG09.masked.variable.vcf.gz \
    -V ../res/gatk/LG10.masked.variable.vcf.gz \
    -V ../res/gatk/LG11.masked.variable.vcf.gz \
    -V ../res/gatk/LG12.masked.variable.vcf.gz \
    -V ../res/gatk/LG13.masked.variable.vcf.gz \
    -V ../res/gatk/LG14.masked.variable.vcf.gz \
    -V ../res/gatk/LG15.masked.variable.vcf.gz \
    -V ../res/gatk/LG16.masked.variable.vcf.gz \
    -V ../res/gatk/LG17.masked.variable.vcf.gz \
    -V ../res/gatk/LG18.masked.variable.vcf.gz \
    -V ../res/gatk/LG19.masked.variable.vcf.gz \
    -V ../res/gatk/LG20.masked.variable.vcf.gz \
    -V ../res/gatk/LG21.masked.variable.vcf.gz \
    -V ../res/gatk/LG22.masked.variable.vcf.gz \
    -V ../res/gatk/LG23.masked.variable.vcf.gz \
    -out ${concatenated_vcfgz} \
    --assumeSorted

