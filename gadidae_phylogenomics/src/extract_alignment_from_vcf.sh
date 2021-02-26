# m_matschiner Sun Jul 22 20:47:55 CEST 2018

## Load the bcftools and ruby modules.
module load bcftools/1.6
module load ruby/2.1.5

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
bcftools view -r "${region}" ${gzvcf_with_relative_path} > region.vcf

# Get the number of sites in the window vcf.
n_sites_calculated=`cat ${info_with_absolute_path} | grep "n_sites:" | wc -l`
if [[ ${n_sites_calculated} == 0 ]]
then
    n_sites=`bcftools view -H region.vcf | wc -l`
    echo "n_sites:${n_sites}" >> ${info_with_absolute_path}
else
    n_sites=`cat ${info_with_absolute_path} | grep "n_sites:" | cut -d ":" -f 2`
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
	mv alignment.phy ${align_with_absolute_path}
	echo "alignment_extracted:yes" >> ${info_with_absolute_path}
    elif [ ! -f ${align_with_absolute_path} ]
    then
	ruby ${script_with_absolute_path} region.vcf ${region} ${min_likelihood} alignment.phy
	if [ -f alignment.phy ]
	then
	    mv alignment.phy ${align_with_absolute_path}
	fi
    fi
fi

# Clean up.
rm -f region.vcf
rm -f alignment.phy