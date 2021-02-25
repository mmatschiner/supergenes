#!/bin/bash

# m_matschiner Mon Nov 4 09:55:39 CET 2019
#
# This script maps resequencing data in fastq.gz format against
# a reference genome. Before running this script, the reference
# file should be indexed with script prepare.sh.
#
# Start in parallel for each individual with
# 'sbatch map.sh reference.fasta XXX_R1.fastq.gz XXX_R2.fastq.gz XXX.bam XXX.log'

# Job name:
#SBATCH --job-name=map
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3936M

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
source /cluster/bin/jobsetup
module load bwa/0.7.12
module load samtools/1.3.1
module load picard/2.6.0 # Has no effect, as picard is a jar file. Full path specified below.
module load gatk/3.7 # Has no effect, as gatk is a jar file. Full path specified below.

## Feedback.
echo
echo "map.sh was called with arguments ${@}."
echo

## Get command line arguments.
fastq1=`readlink -f ${1}`
fastq2=`readlink -f ${2}`
ref=`readlink -f ${3}`
bam=`readlink -f ${4}`
log=`readlink -f ${5}`
sample_id=`basename ${fastq1%_R1.fastq.gz}`
fastq1_id=`basename ${fastq1%.fastq.gz}`
fastq2_id=`basename ${fastq2%.fastq.gz}`

## Change to working directory.
cd ${SCRATCH}

## Map against the reference with bwa mem.
START=$(date +%s)
echo -n "Mapping against reference ${ref} with bwa mem... "
bwa mem -t 16 -M -R "@RG\tID:${sample_id}\tSM:${sample_id}" ${ref} ${fastq1} ${fastq2} > t1.sam 2> ${log}
print_time ${START}

## Converting sam to bam format with samtools.
START=$(date +%s)
echo -n "Converting sam to bam format with samtools... "
samtools view -hbS t1.sam > t1.bam 2>> ${log}
rm t1.sam
if [ ! -f t1.bam ]; then
    echo
    echo "ERROR: File t1.bam does not exist."
    echo "ERROR: Something must have gone wrong with mapping."
    exit 1
fi
print_time $START

## Sort the bam file with samtools.
START=$(date +%s)
echo -n "Sorting the bam file t1.bam with samtools... "
samtools sort -@ 16 -m 3500M -o t2.bam t1.bam 2> ${log}
rm t1.bam
if [ ! -f t2.bam ]; then
    echo
    echo "ERROR: File t2.bam does not exist."
    echo "ERROR: Something must have gone wrong with file sorting."
    exit 1
fi
print_time $START

## Mark duplicates in the bam file with picard-tools.
START=$(date +%s)
echo -n "Marking duplicates in bam file t2.bam... "
java -jar /projects/cees/bin/picard/2.6.0/picard.jar MarkDuplicates \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    I=t2.bam \
    O=t3.bam \
    M=${sample_id}.dupmetrics.txt &>> ${log}
rm t2.bam
if [ ! -f t3.bam ]; then
    echo
    echo "ERROR: File t3.bam does not exist."
    echo "ERROR: Something must have gone wrong with duplicate marking."
    exit 1
fi
print_time ${START}

## Index the bam file file samtools.
START=$(date +%s)
echo -n "Indexing bam file t3.bam with samtools... "
samtools index t3.bam
print_time ${START}

## Run GATK's RealignerTargetCreator as a first step to realign indels.
START=$(date +%s)
echo -n "Creating realignment targets for bam file t3.bam with GATK... "
java -jar /cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${ref} \
    -I t3.bam \
    -o t3.bam.intervals \
    --validation_strictness STRICT &>> ${log}
if [ ! -f t3.bam.intervals ]; then
    echo
    echo "ERROR: File t3.bam.intervals does not exist."
    echo "ERROR: Something must have gone wrong with realignment target creation."
    exit 1
fi
print_time ${START}

## Run GATK's IndelRealigner as a second step to realign indels.
START=$(date +%s)
echo -n "Realigning indels in bam file t3.bam with GATK... "
java -jar /cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${ref} \
    -I t3.bam \
    -targetIntervals t3.bam.intervals \
    -o t4.bam &>> ${log}
rm t3.bam
if [ ! -f t4.bam ]; then
    echo
    echo "ERROR: File t4.bam does not exist."
    echo "ERROR: Something must have gone wrong with indel realignment."
    exit 1
fi
print_time ${START}

## Add read groups.
START=$(date +%s)
echo -n "Adding read groups with picard-tools... "
java -jar /projects/cees/bin/picard/2.6.0/picard.jar AddOrReplaceReadGroups \
    I=t4.bam \
    O=t5.bam \
    RGID=${sample_id} \
    RGLB=${sample_id} \
    RGPL=illumina \
    RGPU=${sample_id} \
    RGSM=${sample_id} \
    VALIDATION_STRINGENCY=LENIENT &>> ${log}
if [ ! -f t5.bam ]; then
    echo
    echo "ERROR: File t5.bam does not exist."
    echo "ERROR: Something must have gone wrong with read group adding."
    exit 1
fi
print_time ${START}

## Index the bam file file samtools.
START=$(date +%s)
echo -n "Indexing bam file t5.bam with samtools... "
samtools index t5.bam &>> ${log}
print_time ${START}

## Copy the bam file and the dupmetrics file back to the analysis directory.
START=$(date +%s)
echo -n "Copying file ${bam} back to the analysis directory... "
cp t5.bam ${bam}
print_time ${START}
START=$(date +%s)
echo -n "Copying file ${bam}.bai back to the analysis directory... "
cp t5.bam.bai ${bam}.bai
print_time ${START}
