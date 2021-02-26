# michaelm Mon Feb 8 00:42:47 CET 2021

# Repeat gatk1 for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/melAeg_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Repeat first with the original gadoga specimen and then with the alternative one.
    for gadoga_setting in original_gadoga alternative_gadoga
    do

        # Generate filtered versions of all vcf files.
        for vcfgz in ../res/gatk/${ref_id}/${gadoga_setting}/LG??.vcf.gz
        do
            lg=`basename ${vcfgz%.vcf.gz}`
            filtered_vcfgz=${vcfgz%.vcf.gz}.filtered.vcf.gz
            if [ ! -f ${filtered_vcfgz} ]
            then
                out=../log/gatk/filter.${ref_id}_${gadoga_setting}_${lg}.out
                log=../log/gatk/filter.${ref_id}_${gadoga_setting}_${lg}.log
                # Remove log files if they exist already.
                rm -f ${out}
                rm -f ${log}
                sbatch -o ${out} apply_filter.slurm ${vcfgz} ${ref} ${filtered_vcfgz} ${log}
            fi
        done
    done
done
