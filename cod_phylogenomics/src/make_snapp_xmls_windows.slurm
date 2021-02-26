#!/bin/bash

# michaelm Thu Mar 5 14:54:16 CET 2020
#
# This script produces xmls files for snapp in sliding windows.

# Start with
# 'sbatch make_snapp_xmls.slurm in.vcf.gz res_dir window_size min_n_sites max_n_sites min_site_dist'

# Job name:
#SBATCH --job-name=mkxml
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1

## Set up the job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.9-foss-2018b
module load Ruby/2.6.1-GCCcore-7.3.0
module list

## Feedback.
echo
echo "make_snapp_xmls.slurm was called with arguments ${@}."
echo

## Get command-line arguments.
gzvcf=`readlink -f ${1}`
res_dir=`readlink -f ${2}`
window_size=${3}
min_n_sites=${4}
max_n_sites=${5}
min_site_dist=${6}
log=`readlink -f ${7}`

# Get the snapp_prep.rb script.
if [ ! -f snapp_prep.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb
fi
cp snapp_prep.rb ${SCRATCH}

# Change to the temporary directory.
cd ${SCRATCH}

# Get the list of samples in the vcf.
lg_id=`basename ${gzvcf} | cut -d "." -f 1`
if [ ${lg_id} == "LG12" ]
then
    samples=`bcftools view -h ${gzvcf} -s ^Gadmor_lfc1,Gadmor_lfc2 | tail -n 1 | cut -f 10-`
else
    samples=`bcftools view -h ${gzvcf} | tail -n 1 | cut -f 10-`
fi

# Write a temporary samples file.                                        
echo -e "species\tspecimen" > tmp.spc.txt
last_sample=""
for sample in ${samples}
do
    echo -e "${sample:0:10}_spc\t${sample}" >> tmp.spc.txt
    last_sample=${sample}
done

# Write a temporary constraints file.
echo -ne "lognormal(0,3.83,0.093)\tcrown\t" > tmp.con.txt # The constraint is from the largest tree topology subset of the aim analysis with a prior distribution on the root age. 
samples_unique=`for sample in ${samples}; do echo ${sample:0:10}; done | sort | uniq`
for sample in ${samples_unique}
do
    echo -n "${sample:0:10}_spc" >> tmp.con.txt
    if [ ! ${sample} == ${last_sample} ]
    then
        echo -n "," >> tmp.con.txt
    fi
done

# Extract window vcf files and produce xml files for each of these.
window_start=1
window_end=$(( window_start + window_size - 1 ))
last_vcf_line=`bcftools view -H ${gzvcf} | tail -n 1`
lg_id=`echo ${last_vcf_line} | tr -s " " | cut -d " " -f 1`
last_pos=`echo ${last_vcf_line} | tr -s " " | cut -d " " -f 2`
while [ ${window_start} -lt ${last_pos} ]
do
    # Get the xml name.
    window_start_padded=${window_start}
    while [ ${#window_start_padded} -ne 9 ];
    do
	    window_start_padded="0"${window_start_padded}
    done
    window_end_padded=${window_end}
    while [ ${#window_end_padded} -ne 9 ];
    do
	    window_end_padded="0"${window_end_padded}
    done
    xml=${res_dir}/${lg_id}_${window_start_padded}_${window_end_padded}.xml
    xml_base=`basename ${xml}`

    # Extract the window data.
    if [ ${lg_id} == "LG12" ]
    then
        bcftools view -s ^Gadmor_lfc1,Gadmor_lfc2 -O v -o tmp.vcf ${gzvcf} ${lg_id}:${window_start}-${window_end}
    else
        bcftools view -O v -o tmp.vcf ${gzvcf} ${lg_id}:${window_start}-${window_end}
    fi

    # Make the xml.
    echo -n "Generating xml file for window ${lg_id}:${window_start}-${window_end}..."
    ruby snapp_prep.rb -v tmp.vcf -t tmp.spc.txt -c tmp.con.txt -l 100000 -m ${max_n_sites} -q ${min_site_dist} -x ${xml} -o ${xml_base%.xml}
    echo " done."

    # Advance the window.
    window_start=$(( window_end + 1 ))
    window_end=$(( window_start + window_size - 1 ))
done
