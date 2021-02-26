# m_matschiner Mon Feb 22 11:40:25 CET 2021

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0

# Quantify reference bias for all genome-wide vcfs.
for ref_id in gadMor2 melAeg_in_gadMor2_coords
do
    for gadoga_setting in original_gadoga alternative_gadoga
    do
        # Set the compressed vcf.
        gzvcf=../res/gatk/${ref_id}/${gadoga_setting}/colinear.vcf.gz

        # Uncompress the vcf.
        gunzip -c ${gzvcf} > tmp.vcf

        # Set the output file.
        table=../res/tables/ref_bias.${ref_id}_${gadoga_setting}.txt
        
        # Calculate the proportion of reads at heterozygous sites that have the reference allele, for all samples.
        ruby quantify_reference_bias.rb tmp.vcf > ${table}
        
        # Clean up.
        rm -f tmp.vcf
    done
done
