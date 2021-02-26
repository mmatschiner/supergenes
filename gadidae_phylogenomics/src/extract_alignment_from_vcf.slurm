#!/bin/bash

# m_matschiner Sun Jul 22 20:47:55 CEST 2018

# Job name:
#SBATCH --job-name=extract
#
# Wall clock limit:
#SBATCH --time=72:00:00
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
module load ruby/2.1.5

## Feedback.
echo
echo "extract_alignment_from_vcf.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
gzvcf_with_relative_path=${1}
region=${2}
from=`echo ${region} | cut -d ":" -f 2 | cut -d "-" -f 1`
to=`echo ${region} | cut -d ":" -f 2 | cut -d "-" -f 2`
align_with_relative_path=${3}
align_with_absolute_path=`readlink -f ${align_with_relative_path}`
info_with_relative_path=${4}
info_with_absolute_path=`readlink -f ${info_with_relative_path}`
script_with_absolute_path=`readlink -f extract_alignment_from_bp_res_vcf.rb`
script=`basename ${script_with_absolute_path}`
chromosome_id=`echo ${region} | cut -d ":" -f 1`
p_missing_allowed=${5}
min_likelihood=${6}

# Touch the info file to make sure that it exists.
touch ${info_with_absolute_path}

# Extracting the region from the vcf file.
bcftools view -r "${region}" ${gzvcf_with_relative_path} > ${SCRATCH}/region.vcf

# Move to the temporary directory.
cd ${SCRATCH}

# Get the number of sites in the window vcf.
n_sites_calculated=`cat ${info_with_absolute_path} | grep "n_sites:" | wc -l`
if [[ ${n_sites_calculated} == 0 ]]
then
    n_sites=`bcftools view -H region.vcf | wc -l`
    echo "n_sites:${n_sites}" >> ${info_with_absolute_path}
else
    n_sites=`cat ${info_with_absolute_path} | grep "n_sites: | cut -d ":" -f 2"`
fi

# Use a ruby script to extract the alignment if there are sufficient sites.
window_size=$(( ${to} - ${from} + 1 ))
n_sites_required=`echo "${window_size}*${p_missing_allowed}" | bc`
if (( $(echo "${n_sites} >= ${n_sites_required}" | bc -l ) ))
then
    alignment_extracted=`cat ${info_with_absolute_path} | grep "alignment_extracted:" | wc -l`
    if [[ ${alignment_extracted} == 0 ]]
    then
	ruby ${script_with_absolute_path} region.vcf ${region} ${min_likelihood} alignment.phy
	cp alignment.phy ${align_with_absolute_path}
	echo "alignment_extracted:yes" >> ${info_with_absolute_path}
    fi
fi
