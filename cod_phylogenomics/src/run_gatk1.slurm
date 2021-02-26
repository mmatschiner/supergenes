#!/bin/bash

# michaelm Wed Feb 10 15:46:50 CET 2021
#
# This script performs per individual variant calling with GATK's
# HaplotypeCaller.
#
# Start in parallel for each individual with
# 'sbatch run_gatk1.slurm reference.fasta XXX.merged.sorted.dedup.realn.bam XXX.gvcf'

# Job name:
#SBATCH --job-name=gatk1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G

# Define a function to print the time.
function print_time {
    END=$(date +%s)
    DIFF=$(( $END - $1 ))
    dd=$(echo "$DIFF/86400" | bc)
    dt2=$(echo "$DIFF-86400*$dd" | bc)
    dh=$(echo "$dt2/3600" | bc)
    dt3=$(echo "$dt2-3600*$dh" | bc)
    dm=$(echo "$dt3/60" | bc)
    ds=$(echo "$dt3-60*$dm" | bc)
    if [ $dd -gt 0 ]; then
        echo "done in ${dd} days and ${dh} hours."
    elif [ $dh -gt 0 ]; then
        echo "done in ${dh} hours and ${dm} minutes."
    elif [ $dm -gt 0 ]; then
        echo "done in ${dm} minutes and ${ds} seconds."
    else
        echo "done in ${ds} seconds."
    fi
}

## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-intel-2018b
module swap GCCcore/7.3.0 GCCcore/8.2.0
module swap zlib/1.2.11-GCCcore-7.3.0 zlib/1.2.11-GCCcore-8.2.0
module swap binutils/2.30-GCCcore-7.3.0 binutils/2.31.1-GCCcore-8.2.0
module swap bzip2/1.0.6-GCCcore-7.3.0 bzip2/1.0.6-GCCcore-8.2.0
module swap XZ/5.2.4-GCCcore-7.3.0 XZ/5.2.4-GCCcore-8.2.0
module swap cURL/7.60.0-GCCcore-7.3.0 cURL/7.63.0-GCCcore-8.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
module list

## Feedback.
echo
echo "run_gatk1.slurm was called with arguments ${@}."
echo

## Get command line arguments.
ref_with_relative_path=$1
bam_with_relative_path=$2
gvcf_with_relative_path=$3
log_with_relative_path=$4
chromosome_id=$5

# Set variables based on arguments.
ref_with_absolute_path=`readlink -e ${ref_with_relative_path}`
bam_with_absolute_path=`readlink -e ${bam_with_relative_path}`
gvcf_with_absolute_path=`readlink -f ${gvcf_with_relative_path}`
gvcf_base=`basename ${gvcf_with_absolute_path}`
log_with_absolute_path=`readlink -f ${log_with_relative_path}`

# Change to working directory.
cd $SCRATCH

# Index the bam file.
samtools index ${bam_with_absolute_path}

# Run gatk's haplotypecaller.
START=$(date +%s)
echo -n "Calling variants in bam file ${bam_with_absolute_path} with GATK's HaplotypeCaller... "
gatk HaplotypeCaller \
    -I ${bam_with_absolute_path} \
    -O ${gvcf_base}.g.vcf \
    -R ${ref_with_absolute_path} \
    -L ${chromosome_id} \
    --emit-ref-confidence BP_RESOLUTION \
    --min-base-quality-score 15 \
    --pcr-indel-model NONE &> run_gatk1.${gvcf_base}.log
print_time $START

# Compress the variant call vcf file.
START=$(date +%s)
echo -n "Compressing file ${gvcf_base}.g.vcf... "
bgzip ${gvcf_base}.g.vcf
print_time $START

# Index the compressed variant call vcf file.
START=$(date +%s)
echo -n "Indexing file ${gvcf_base}.g.vcf.gz..."
tabix -p vcf ${gvcf_base}.g.vcf.gz
print_time $START

# Copy the compressed variant call vcf file back to the analysis directory.
START=$(date +%s)
echo -n "Copying file ${gvcf_base}.g.vcf.gz to the output directory... "
cp ${gvcf_base}.g.vcf.gz ${gvcf_with_absolute_path%.gz}.gz
print_time $START

# Copy the index file back to the analysis directory.
START=$(date +%s)
echo -n "Copying file ${gvcf_base}.g.vcf.gz.tbi to the output directory... "
cp ${gvcf_base}.g.vcf.gz.tbi ${gvcf_with_absolute_path%.gz}.gz.tbi
print_time $START

# Copy the log file to the log directory.
START=$(date +%s)
echo -n "Copying file run_gatk1.${gvcf_base}.log back to the submission directory... "
cp run_gatk1.${gvcf_base}.log ${log_with_absolute_path}
print_time $START
echo
