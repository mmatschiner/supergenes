#!/bin/bash

# m_matschiner Tue Apr 24 12:21:44 CEST 2018

# Job name:
#SBATCH --job-name=abba
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=50G

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
module load python2/2.7.10
module load bcftools/1.6
module load R/3.4.4

## Feedback.
echo
echo "calculate_abba_baba.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
gzvcf_prefix_with_relative_path=${1}
gzvcf_prefix_with_absolute_path=`readlink -f ${gzvcf_prefix_with_relative_path}`
species=${2}
spc_o=`echo ${species} | cut -d "," -f 4`
spc_p3=`echo ${species} | cut -d "," -f 3`
spc_p2=`echo ${species} | cut -d "," -f 2`
spc_p1=`echo ${species} | cut -d "," -f 1`
summary_file_with_relative_path=${3}
summary_file_with_absolute_path=`readlink -f ${summary_file_with_relative_path}`

# Feedback.
echo "Input vcf file prefix: ${gzvcf_prefix_with_absolute_path}."
echo "Species string: ${species}."
echo

# Copy r scripts to the temporary directory.
cp calculate_abba_baba.r ${SCRATCH}

# Move to the temporary directory.
cd ${SCRATCH}

# Feedback.
echo "Outgroup: ${spc_o}"
echo "P3: ${spc_p3}"
echo "P2: ${spc_p2}"
echo "P1: ${spc_p1}"
echo

# Make a list of approximate chromosome lengths, using the locations of the last snps.
echo -n "Making a list of chromosome lengths..."
START=$(date +%s)
rm -f lg_lengths.txt
for n in `seq -w 1 23`
do
    echo -ne "LG${n}\t" >> lg_lengths.txt
    zcat ${gzvcf_prefix_with_absolute_path}${n}.vcf.gz | tail -n 1 | cut -f 2 >> lg_lengths.txt
done
print_time $START

# Make a combined vcf for bi-allelic snps without missing data in the four species.
echo -n "Making a combined vcf file for all chromosomes..."
START=$(date +%s)
bcftools view -h -s "${spc_p1},${spc_p2},${spc_p3},${spc_o}" ${gzvcf_prefix_with_absolute_path}01.vcf.gz | grep -v GmG20150304_scaffold > tmp.vcf
for n in `seq -w 1 23`
do
    bcftools view -s "${spc_p1},${spc_p2},${spc_p3},${spc_o}" -a --min-ac=1 -v snps -m 2 -M 2 ${gzvcf_prefix_with_absolute_path}${n}.vcf.gz | bcftools view -i 'AN=8' | bcftools view -H -e 'AC==AN' >> tmp.vcf
done
bcftools view -O z -o tmp.vcf.gz tmp.vcf
rm tmp.vcf
print_time $START

# Get simon martin's genomics_general code repository.
echo -n "Downloading simon martin's code..."
START=$(date +%s)
git clone https://github.com/simonhmartin/genomics_general.git
PYTHONPATH=${PYTHONPATH}:./genomics_general
print_time $START

# Convert the vcf to genotype format.
echo -n "Converting the vcf to geno format..."
START=$(date +%s)
python genomics_general/VCF_processing/parseVCF.py -i tmp.vcf.gz | gzip > tmp.geno.gz
print_time $START

# Make a samples file.
rm -f samples.txt
echo -e "${spc_p1}\t${spc_p1}" >> samples.txt
echo -e "${spc_p2}\t${spc_p2}" >> samples.txt
echo -e "${spc_p3}\t${spc_p3}" >> samples.txt
echo -e "${spc_o}\t${spc_o}" >> samples.txt

# Calculate allele frequencies.
echo -n "Calculating allele frequencies..."
START=$(date +%s)
python genomics_general/freq.py -g tmp.geno.gz -p ${spc_p1} -p ${spc_p2} -p ${spc_p3} -p ${spc_o} --popsFile samples.txt --target derived 2> /dev/null | grep -v nan | gzip > allele_freqs.tsv.gz
print_time $START

# Calculate abba-baba statistics across the whole genome.
echo -n "Calculating abba-baba statistics across the whole genome..."
START=$(date +%s)
Rscript calculate_abba_baba.r allele_freqs.tsv.gz abba_baba.summary.txt ${spc_p1} ${spc_p2} ${spc_p3} ${spc_o} lg_lengths.txt
print_time $START

# Copy the abba-baba summary file back to the results directory.
echo -n "Copying the result file to ${summary_file_with_absolute_path}..."
START=$(date +%s)
cp abba_baba.summary.txt ${summary_file_with_absolute_path}
print_time $START
