#!/bin/bash

# m_matschiner Mon Sep 24 13:29:03 CEST 2018
#
# This script performs joint variant calling with GATK's
# GenotypeGVCFs, using compressed vcf files as input.
#
# Start with
# 'sbatch run_gatk2.slurm reference.fasta *.g.vcf.gz'

# Job name:
#SBATCH --job-name=gatk2
#
# Wall clock limit:
#SBATCH --time=24:00:00 ## 48 or 24.
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=60G ## 60 or 4.
##           #SBATCH --cpus-per-task=16 ## in or out.

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
module load BCFtools/1.9-foss-2018b
module load SAMtools/1.9-foss-2018b
module load GATK/4.1.0.0-foss-2018b-Python-3.6.6
module load picard/2.18.27-Java-1.8
module list

## Feedback.
echo
echo "run_gatk2.slurm was called with arguments ${@}."
echo

## Get command line arguments.
chromosome_id=$1
shift
ref_with_relative_path=$1
shift
vcf_with_relative_path=$1
shift
log_with_relative_path=$1
shift
compressed_gvcfs_with_absolute_paths=()
file_count=0
for file_name in "$@"
do
    file_count=$((file_count+1))
    compressed_gvcf_with_absolute_path=`readlink -f ${file_name}`
    compressed_gvcfs_with_absolute_paths+=(${compressed_gvcf_with_absolute_path})
done
variant_string=`for gvcf in "${compressed_gvcfs_with_absolute_paths[@]}"; do echo "-V ${gvcf} "; done`

## Set variables based on arguments.
vcf_with_absolute_path=`readlink -f ${vcf_with_relative_path}`
log_with_absolute_path=`readlink -f ${log_with_relative_path}`
ref=`basename ${ref_with_relative_path}`

## Copy the reference to the work directory.
START=$(date +%s)
echo -n "Copying file ${ref}... "
cp ${ref_with_relative_path} $SCRATCH
print_time $START

## If necessary, prepare the index file for the reference.
if [ ! -f ${ref_with_relative_path}.fai ]
then
    START=$(date +%s)
    echo -n "Preparing index file for ${ref}... "
    samtools faidx ${ref_with_relative_path}
    print_time $START
fi

## Copy the index of the reference to the work directory.
START=$(date +%s)
echo -n "Copying file ${ref}.fai... "
cp ${ref_with_relative_path}.fai $SCRATCH
print_time $START

## If necessary, prepare the dictionary file for the reference.
if [ ! -f ${ref_with_relative_path%.fasta}.dict ]
then
    START=$(date +%s)
    echo -n "Preparing dictionary file for ${ref}... "
    java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
	R=${ref_with_relative_path} \
	O=${ref_with_relative_path%.fasta}.dict \
	URI=./${ref_with_relative_path}
    print_time $START
fi

## Copy the dictionary of the reference to the work directory. 
START=$(date +%s)
echo -n "Copying file ${ref%.fasta}.dict... "
cp ${ref_with_relative_path%.fasta}.dict $SCRATCH
print_time $START

## Change to working directory.
cd $SCRATCH

## Run GATK's GenomicsDBImport.
START=$(date +%s)
echo -n "Combining ${file_count} gvcf files with GATK's GenomicsDBImport... "
gatk GenomicsDBImport \
  ${variant_string} \
  -L ${chromosome_id} \
  --genomicsdb-workspace-path my_database &> gatk2.log
print_time $START

## Run GATK's GenotypeGVCFs.
START=$(date +%s)
echo -n "Joint variant calling for ${file_count} gvcf files with GATK's GenotypeGVCFs... "
variant_string=`for gvcf in "${compressed_gvcfs_with_absolute_paths[@]}"; do echo "--variant ${gvcf} "; done`
gatk GenotypeGVCFs \
  -R ${ref} \
  -V gendb://my_database \
  -new-qual \
  -O gatk.vcf.gz # &> gatk2.log
print_time $START

## Index the vcf.gz file.
START=$(date +%s)
echo -n "Indexing the vcf file... "
bcftools index gatk.vcf.gz
print_time $START

## Copy the variant call vcf.gz file and its index file back to the analysis directory.
START=$(date +%s)
echo -n "Copying the compressed file to the results directory... "
cp gatk.vcf.gz ${vcf_with_absolute_path}
print_time $START
START=$(date +%s)
echo -n "Copying the index file to the results directory... "
cp gatk.vcf.gz.csi ${vcf_with_absolute_path}.csi
print_time $START

## Copy the log file back to the submission directory.
START=$(date +%s)
echo -n "Copying the log file to the log directory... "
cp gatk2.log ${log_with_absolute_path}
print_time $START
echo
