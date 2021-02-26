#!/bin/bash

# m_matschiner Fri Sep 1 17:09:12 CEST 2017

# Job name:
#SBATCH --job-name=mask
#
# Wall clock limit:
#SBATCH --time=100:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G

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

## Set up the job environment.
source /cluster/bin/jobsetup
module load bcftools/1.6

## Feedback.
echo
echo "remove_monomorphic_sites.slurm was called with arguments ${@}."
echo

## Get command line arguments.
vcfgz=`readlink -f ${1}`
vcfgzcsi=`readlink -f ${vcfgz}.csi`
ref=`readlink -f ${2}`
variable_vcfgz=`readlink -f ${3}`
log=`readlink -f ${4}`

## Make sure the input file and its index file exist.
if [ ! -f ${vcfgz} ]
then
    echo
    echo "ERROR: The file ${vcfgz} does not exist!"
    exit 1
fi
if [ ! -f ${vcfgzcsi} ]
then
    START=$(date +%s)
    echo -n "Indexing file ${vcfgz}... "
    bcftools index ${vcfgz}
    if [ ! -f ${vcfgzcsi} ]
    then
	echo
	echo "ERROR: The file ${vcfgzcsi} does not exist!"
	echo "ERROR: Something must have gone wrong with indexing."
	exit 1
    fi
    print_time $START
fi

## Exit if the output file exists already.
if [ -f ${variable_vcfgz} ]
then
    echo
    echo "ERROR: The output file ${variable_vcfgz} exists already."
    exit 1
fi

## Remove the log file if it exists already.
rm -f ${log}

## Touch the index file to turn off the bcftools warning that it is older.
touch ${vcfgzcsi}

## Change to the temporary directory.
cd ${SCRATCH}

## Removing sites with low quality with bcftools.
START=$(date +%s)
echo -n "Removing monomorphic sites... "
bcftools view --min-ac=1 -e 'AC==0 || AC==AN' -o t1.vcf ${vcfgz} &> ${log}
if [ ! -f t1.vcf ]; then
  echo
  echo "ERROR: File t1.vcf does not exist."
  echo "ERROR: Something must have gone wrong with filtering for paralogy."
  exit 1
fi
length_t0=`bcftools view -H ${vcfgz} | wc -l`
length_t1=`bcftools view -H t1.vcf | wc -l`
print_time $START
echo "-> Removed $(( ${length_t0} - ${length_t1} )) monomorphic sites (before: ${length_t0}, after: ${length_t1})."

## Compress with bcftools.
START=$(date +%s)
echo -n "Compressing file ${variable_vcfgz}... "
bcftools view -O z -o ${variable_vcfgz} t1.vcf 2> ${log}
if [ ! -f ${variable_vcfgz} ]; then
  echo
  echo "ERROR: File ${variable_vcfgz} does not exist."
  echo "ERROR: Something must have gone wrong with compressing."
  exit 1
fi
rm t1.vcf
print_time $START

## Index with bcftools.
START=$(date +%s)
echo -n "Indexing file ${variable_vcfgz}... "
bcftools index ${variable_vcfgz} 2> ${log}
if [ ! -f ${variable_vcfgz}.csi ]; then
  echo
  echo "ERROR: File ${variable_vcfgz}.csi does not exist."
  echo "ERROR: Something must have gone wrong with indexing."
  exit 1
fi
print_time $START
