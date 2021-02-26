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
module load vcftools/0.1.14.zlib.1.2.8
module load ruby/2.1.5

## Feedback.
echo
echo "apply_mask.slurm was called with arguments ${@}."
echo

## Get command line arguments.
vcfgz_with_relative_path=${1}
vcfgzcsi_with_relative_path="${vcfgz_with_relative_path}.csi"
ref_with_relative_path=${2}
mappability_bed_with_relative_path=${3}
masked_vcfgz_with_relative_path=${4}
indel_bed_with_relative_path="${masked_vcfgz_with_relative_path%.vcf.gz}.indels.bed"
log_with_relative_path=${5}

## Make sure the input file and its index file exist.
if [ ! -f ${vcfgz_with_relative_path} ]
then
  echo
  echo "ERROR: The file ${vcfgz_with_relative_path} does not exist!"
  exit 1
fi
if [ ! -f ${vcfgzcsi_with_relative_path} ]
then
  START=$(date +%s)
  echo -n "Indexing file ${vcfgz_with_relative_path}... "
  bcftools index ${vcfgz_with_relative_path}
  if [ ! -f ${vcfgzcsi_with_relative_path} ]
  then
    echo
    echo "ERROR: The file ${vcfgzcsi_with_relative_path} does not exist!"
    echo "ERROR: Something must have gone wrong with indexing."
    exit 1
  fi
  print_time $START
fi

## Get the absolute file paths.
vcfgz=`basename ${vcfgz_with_relative_path}`
lg=${vcfgz%.vcf.gz}
vcfgz_with_absolute_path=`readlink -e ${vcfgz_with_relative_path}`
vcfgzcsi=`basename ${vcfgzcsi_with_relative_path}`
vcfgzcsi_with_absolute_path=`readlink -e ${vcfgzcsi_with_relative_path}`
ref_with_absolute_path=`readlink -e ${ref_with_relative_path}`
mappability_bed=`basename ${mappability_bed_with_relative_path}`
masked_vcfgz=`basename ${masked_vcfgz_with_relative_path}`
masked_vcfgz_with_absolute_path=`readlink -f ${masked_vcfgz_with_relative_path}`
indel_bed=`basename ${indel_bed_with_relative_path}`
indel_bed_with_absolute_path=`readlink -f ${indel_bed_with_relative_path}`
log=`basename ${log_with_relative_path}`
log_with_absolute_path=`readlink -f ${log_with_relative_path}`

## Exit if the output file exists already.
if [ -f ${masked_vcfgz_with_absolute_path} ]
then
  echo
  echo "ERROR: The output file ${masked_vcfgz_with_absolute_path} exists already."
  exit 1
fi

## Remove the log file if it exists already.
rm -f ${log_with_relative_path}

## Copy ruby scripts to the temporary directory.
cp make_indel_mask.rb ${SCRATCH}
cp replace_spanning_deletions.rb ${SCRATCH}

## Copy the mappability mask to the temporary directory.
echo -e "chrom\tchromStart\tchromEnd" > ${SCRATCH}/${mappability_bed}
cat ${mappability_bed_with_relative_path} | sed 's/lg/LG/g' >> ${SCRATCH}/${mappability_bed}

## Touch the index file to turn off the bcftools warning that it is older.
touch ${vcfgzcsi_with_absolute_path}

## Change to the temporary directory.
cd ${SCRATCH}

## Removing sites with low quality with bcftools.
START=$(date +%s)
echo -n "Removing sites with potential paralogy... "
bcftools filter -i 'DP>40 && DP<400' -o t1.vcf ${vcfgz_with_absolute_path} &> ${log_with_absolute_path}
# bcftools filter -i 'DP>40 && DP<400 && InbreedingCoeff>-0.5' -o t1.vcf ${vcfgz_with_absolute_path} &> ${log_with_absolute_path}
if [ ! -f t1.vcf ]; then
  echo
  echo "ERROR: File t1.vcf does not exist."
  echo "ERROR: Something must have gone wrong with filtering for paralogy."
  exit 1
fi
length_t0=`bcftools view -H ${vcfgz_with_absolute_path} | wc -l`
length_t1=`bcftools view -H t1.vcf | wc -l`
print_time $START
echo "-> Removed $(( ${length_t0} - ${length_t1} )) sites (before: ${length_t0}, after: ${length_t1})."

## Removing sites with low quality with bcftools.
START=$(date +%s)
echo -n "Removing sites with low quality... "
bcftools filter -e 'ReadPosRankSum<-0.5 || MQRankSum<-0.5' -o t2.vcf t1.vcf &> ${log_with_absolute_path}
if [ ! -f t2.vcf ]; then
  echo
  echo "ERROR: File t2.vcf does not exist."
  echo "ERROR: Something must have gone wrong with filtering for quality."
  exit 1
fi
length_t2=`bcftools view -H t2.vcf | wc -l`
rm t1.vcf
print_time $START
echo "-> Removed $(( ${length_t1} - ${length_t2} )) sites (before: ${length_t1}, after: ${length_t2})."

## Normalize indels.
START=$(date +%s)
echo -n "Normalizing indels... "
bcftools norm -f ${ref_with_absolute_path} -o t3.vcf t2.vcf &> ${log_with_absolute_path}
if [ ! -f t3.vcf ]; then
  echo
  echo "ERROR: File t3.vcf does not exist."
  echo "ERROR: Something must have gone wrong with indel normalization."
  exit 1
fi
rm t2.vcf
print_time $START

## Separate SNP-only and indel sites.
START=$(date +%s)
echo -n "Separating indel and no-indel sites... "
bcftools view -e 'TYPE~"indel"' -o t3.no_indels.vcf t3.vcf
bcftools view -i 'TYPE~"indel"' -o t3.indels.vcf t3.vcf
if [ ! -f t3.indels.vcf ]; then
  echo
  echo "ERROR: File t3.indels.vcf does not exist."
  echo "ERROR: Something must have gone wrong with extracting indel sites."
  exit 1
fi
if [ ! -f t3.no_indels.vcf ]; then
  echo
  echo "ERROR: File t3.no_indels.vcf does not exist."
  echo "ERROR: Something must have gone wrong with removing indel sites."
  exit 1
fi
print_time $START

## Remove indel sites with minor allele count of 2 or less. 
START=$(date +%s)
echo -n "Removing indels with low minor allele count... "
bcftools filter -S . -i 'FMT/DP>=6 && FMT/GQ>=10' t3.indels.vcf | bcftools view --min-ac=3 -o t4.indels.vcf &> ${log_with_absolute_path}
if [ ! -f t4.indels.vcf ]; then
  echo
  echo "ERROR: File t4.indels.vcf does not exist."
  echo "ERROR: Something must have gone wrong with removing indel sites with low minor allele count."
  exit 1
fi
length_t3i=`bcftools view -H t3.indels.vcf | wc -l`
length_t4i=`bcftools view -H t4.indels.vcf | wc -l`
rm t3.indels.vcf
print_time $START
echo "-> Removed $(( ${length_t3i} - ${length_t4i} )) indels (before: ${length_t3i}, after: ${length_t4i})."

## Convert the indel vcf file into a list of indels.
START=$(date +%s)
echo -n "Converting the indel vcf file into a list... "
bcftools view -H t4.indels.vcf | cut -f 1,2,4,5 > t4.indels.txt
if [ ! -f t4.indels.txt ]; then
  echo
  echo "ERROR: File t4.indels.txt does not exist."
  echo "ERROR: Something must have gone wrong with conversion of the indel file.."
  exit 1
fi
rm t4.indels.vcf
print_time $START

## Make a mask for sites in proximity to indels.
START=$(date +%s)
echo -n "Making a mask for sites in proximity to indels... "
feedback=`ruby make_indel_mask.rb t4.indels.txt ${indel_bed} ${minimum} ${maximum}`
if [ ! -f ${indel_bed} ]; then
  echo
  echo "ERROR: File ${indel_bed} does not exist."
  echo "ERROR: Something must have gone wrong with making a mask for sites in proximity to indels."
  exit 1
fi
print_time $START
echo "-> ${feedback}"

## Apply the mask for sites in proximity to indels.
START=$(date +%s)
echo -n "Applying the mask for sites in proximity to indels... "
vcftools --vcf t3.no_indels.vcf --exclude-bed ${indel_bed} --recode --recode-INFO-all --out t4 2> ${log_with_absolute_path}
if [ ! -f t4.recode.vcf ]; then
  echo
  echo "ERROR: File t4.recode.vcf does not exist."
  echo "ERROR: Something must have gone wrong with applying the mask for sites in proximity to indels."
  exit 1
fi
length_t3=`bcftools view -H t3.no_indels.vcf | wc -l`
length_t4=`bcftools view -H t4.recode.vcf | wc -l`
rm t3.no_indels.vcf
print_time $START
mv t4.recode.vcf t4.vcf
echo "-> Removed $(( ${length_t3} - ${length_t4} )) sites (before: ${length_t3s}, after: ${length_t4})."

## Apply the mask for mappable sites.
START=$(date +%s)
echo -n "Applying the mask for mappable sites... "
vcftools --vcf t4.vcf --exclude-bed ${mappability_bed} --recode --recode-INFO-all --out t5 2> ${log_with_absolute_path}
if [ ! -f t5.recode.vcf ]; then
  echo
  echo "ERROR: File t5.recode.vcf does not exist."
  echo "ERROR: Something must have gone wrong with applying the mask for mappable sites."
  exit 1
fi
length_t5=`bcftools view -H t5.recode.vcf | wc -l`
rm t4.vcf
print_time $START
echo "-> Removed $(( ${length_t4} - ${length_t5} )) sites (before: ${length_t4}, after: ${length_t5})."
mv t5.recode.vcf t5.vcf

## Replace alleles coding for spanning deletions with missing data.
START=$(date +%s)
echo -n "Replace alleles coding for spanning deletions... "
ruby replace_spanning_deletions.rb t5.vcf t6.vcf
if [ ! -f t6.vcf ]; then
  echo
  echo "ERROR: File t6.vcf does not exist."
  echo "ERROR: Something must have gone wrong with replacing spanning deletions."
  exit 1
fi
rm t5.vcf
print_time $START

## Compress with bcftools.
START=$(date +%s)
echo -n "Compressing file ${masked_vcfgz}... "
bcftools view -O z -o ${masked_vcfgz} t6.vcf 2> ${log_with_absolute_path}
if [ ! -f ${masked_vcfgz} ]; then
  echo
  echo "ERROR: File ${masked_vcfgz} does not exist."
  echo "ERROR: Something must have gone wrong with compressing."
  exit 1
fi
rm t6.vcf
print_time $START

## Index with bcftools.
START=$(date +%s)
echo -n "Indexing file ${masked_vcfgz}... "
bcftools index ${masked_vcfgz} 2> ${log_with_absolute_path}
if [ ! -f ${masked_vcfgz}.csi ]; then
  echo
  echo "ERROR: File ${masked_vcfgz}.csi does not exist."
  echo "ERROR: Something must have gone wrong with indexing."
  exit 1
fi
print_time $START

## Copy the vcf.gz file, its index file, and its bed file back to the analysis directory.
START=$(date +%s)
echo -n "Copying file ${masked_vcfgz} back to the analysis directory... "
cp ${masked_vcfgz} ${masked_vcfgz_with_absolute_path}
print_time $START
START=$(date +%s)
echo -n "Copying file ${masked_vcfgz}.csi back to the analysis directory... "
cp ${masked_vcfgz}.csi ${masked_vcfgz_with_absolute_path}.csi
print_time $START
START=$(date +%s)
echo -n "Copying file ${indel_bed} back to the analysis directory... "
cp ${indel_bed} ${indel_bed_with_absolute_path}
print_time $START
