#!/bin/bash

# m_matschiner Sun Jul 15 23:51:14 CEST 2018

# This script produces a fasta consensus sequence from a bam
# file.
# Start e.g. with
# 'bash convert_bam_to_fasta.slurm in.bam reference.fasta out.fasta'

# Job name:
#SBATCH --job-name=bam2fas
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G

## Set up job environment.
source /cluster/bin/jobsetup
module load bcftools
module load samtools
module load vcftools
module load picard-tools
module load seqtk

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

# Get the command line arguments.
bam_with_relative_path=$1
ref_with_relative_path=$2
fasta_with_relative_path=$3

## Set variables based on arguments.
bam_with_absolute_path=`readlink -e ${bam_with_relative_path}`
bam=`basename ${bam_with_absolute_path}`
bam_absolute_path=${bam_with_absolute_path%${bam}}
bam_extension=${bam##*.}
bam_base=${bam%.$bam_extension}
ref_with_absolute_path=`readlink -e ${ref_with_relative_path}`
ref=`basename ${ref_with_absolute_path}`
ref_absolute_path=${ref_with_absolute_path%${ref}}
ref_extension=${ref##*.}
ref_base=${ref%.$ref_extension}

## Make pilup file, recalculate BAQ and exclude reads with mapping quality <30 and base phred score <30. Dont exclude orphans (no read pair), and output uncompressed.
echo "Running mpileup... "
START=$(date +%s)
samtools mpileup -E -q 30 -Q 30 -A -u -f ${ref_absolute_path}${ref} ${bam_absolute_path}${bam} > ${bam_base}.pileup
print_time $START

## Create VCF consensus file with bcftools.
echo "Calling consensus sequence to VCF... "
START=$(date +%s)
bcftools call -c -o ${bam_base}.vcf ${bam_base}.pileup
print_time $START

## Remove all sites where these is deletions in the mapped species.
echo "Removing deleted sites..."
START=$(date +%s)
sed -i '/DP=0/d' ${bam_base}.vcf
print_time $START

## Convert the VCF file to fastq format.
echo "Converting VCF to fastq format... "
vcfutils.pl vcf2fq ${bam_base}.vcf > ${bam_base}.fq
print_time $START

## Convert fastq to fasta format.
echo "Converting fastq to fasta format... "
seqtk seq -A ${bam_base}.fq > ${fasta}
print_time $START

## Clean up.
rm ${bam_base}.pileup
rm ${bam_base}.vcf
rm ${bam_base}.fq
