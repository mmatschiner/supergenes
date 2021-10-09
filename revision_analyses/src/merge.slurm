#!/bin/bash

# m_matschiner Fri Nov 22 11:42:08 CET 2019
#
# This script merges multiple bam files for the same sample/
# individual, and performs deduplication and indel realignment
# on the merged file.
#
# Start in parallel for each individual with
# 'sbatch merge.sh merged.bam reference.fasta out.log source1.bam source2.bam'

# Job name:
#SBATCH --job-name=merge
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=7680M

# Project.
#SBATCH --account=nn9244k

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

## Set up job environment.
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BWA/0.7.17-GCC-8.2.0-2.31.1 # bwa/0.7.12
module load SAMtools/1.9-GCC-8.2.0-2.31.1 # samtools/1.3.1
module load picard/2.18.27-Java-1.8 # picard/2.6.0

## Feedback.
echo
echo "merge.slurm was called with arguments ${@}."
echo

## Get command line arguments.
merged_bam=`readlink -f ${1}`
shift
ref=`readlink -f ${1}`
shift
log=`readlink -f ${1}`
shift
bams=()
for bam_w_rel_path in "$@"; do
    bam=`readlink -f ${bam_w_rel_path}`
    bams+=(${bam})
done

## Copy input files to the work directory.
file_count=0
for bam in ${bams[@]}; do
    file_count=$((file_count+1))
    START=$(date +%s)
    echo -n "Copying file ${bam}... "
    cp ${bam} ${SCRATCH}
    print_time ${START}
done

## Change to working directory.
cd ${SCRATCH}

## Merging bam files with picard-tools.
START=$(date +%s)
echo -n "Merging ${file_count} bam files with picard-tools... "
input_string=`for bam in *.bam ; do echo -n "INPUT=${bam} " ; done`
java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
    ${input_string} \
    OUTPUT=t1.bam &> ${log}
if [ ! -f t1.bam ]; then
    echo
    echo "ERROR: File t1.bam does not exist."
    echo "ERROR: Something must have gone wrong with file merging."
    exit 1
fi
print_time ${START}

## Sorting the bam file with samtools.
START=$(date +%s)
echo -n "Sorting the bam file t1.bam with samtools... "
samtools sort -@ 8 -m 7000M -o t2.bam t1.bam 2> ${log}
rm t1.bam
if [ ! -f t2.bam ]; then
    echo
    echo "ERROR: File t2.bam does not exist."
    echo "ERROR: Something must have gone wrong with file sorting."
    exit 1
fi
print_time ${START}

## Marking duplicates in the bam file with picard-tools.
START=$(date +%s)
echo -n "Marking duplicates in t2.bam... "
java -Xmx7g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    I=t2.bam \
    O=t3.bam \
    M=t3.dupmetrics.txt &> ${log}
rm t2.bam
if [ ! -f t3.bam ]; then
    echo
    echo "ERROR: File t3.bam does not exist."
    echo "ERROR: Something must have gone wrong with duplicate marking."
    exit 1
elif [ ! -f t3.dupmetrics.txt ]; then
    echo
    echo "ERROR: File t3.dupmetrics.txt does not exist."
    echo "ERROR: Something must have gone wrong with duplicate marking."
    exit 1
fi
print_time ${START}

## Index the bam file file with samtools.
START=$(date +%s)
echo -n "Indexing bam file t3.bam with samtools... "
samtools index t3.bam > ${log}
print_time ${START}

## Validate the bam file with picard-tools.
echo -n "Validating bam file t3.bam... "
java -Xmx7g -jar $EBROOTPICARD/picard.jar ValidateSamFile \
    MAX_OPEN_TEMP_FILES=1000 \
    I=t3.bam \
    O=t3.bam.valid &> ${log}
if [ ! -f t3.bam.valid ]; then
    echo
    echo "ERROR: File t3.bam.valid does not exist."
    echo "ERROR: Something must have gone wrong with bam file validation."
    exit 1
elif [ "$(cat t3.bam.valid)" != "No errors found" ]; then
    echo
    echo "ERROR: File t3.bam is not valid:"
    echo "$(cat t3.bam.valid)"
    exit 1
fi
print_time ${START}

## Copying the bam file back to the analysis directory.
START=$(date +%s)
echo -n "Copying file t3.bam to ${merged_bam}... "
cp t3.bam ${merged_bam}
cp t3.bam.bai ${merged_bam}.bai
print_time ${START}
