# michaelm Mon Feb 8 00:40:08 CET 2021

# Make the prepare log directory.
mkdir -p ../log/prepare

# Define a function to sleep if any unc_idx jobs are queued or running.
function sleep_while_busy {
    n_unc_idx_jobs=`squeue -u michaelm | grep unc_idx | wc -l`
    while [ $n_unc_idx_jobs -gt 0 ]
    do
        sleep 10
        n_unc_idx_jobs=`squeue -u michaelm | grep unc_idx | wc -l`
    done
}

# Repeat gatk2 for each reference.
for ref in ../data/assemblies/gadMor3.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Repeat only with the original gadoga specimen.
    for gadoga_setting in original_gadoga
    do

        # Make the result directories.
        mkdir -p ../res/gatk/${ref_id}/${gadoga_setting}

        # Run run_gatk2.slurm with vcf files of all chromosomes.
        for i in `seq 23`
        do
            chromosome_id=chr${i}
            gzvcf="../res/gatk/${ref_id}/${gadoga_setting}/${chromosome_id}.vcf.gz"
            log="../log/gatk/run_gatk2.${ref_id}_${gadoga_setting}_${chromosome_id}.log"
            out="../log/gatk/run_gatk2.${ref_id}_${gadoga_setting}_${chromosome_id}.out"
            if [ ${gadoga_setting} == "original_gadoga" ]
            then
                gzgvcf_string=`for i in ../res/gatk/${ref_id}/??????.${chromosome_id}.g.vcf.gz ../res/gatk/${ref_id}/??????_????.${chromosome_id}.g.vcf.gz; do echo -n "${i}" | grep -v Gadog2; done`
            elif [ ${gadoga_setting} == "alternative_gadoga" ]
            then
                gzgvcf_string=`for i in ../res/gatk/${ref_id}/??????.${chromosome_id}.g.vcf.gz ../res/gatk/${ref_id}/??????_????.${chromosome_id}.g.vcf.gz; do echo -n "${i}" | grep -v Gadoga; done`
            else
                echo "ERROR: Unknown gadoga_setting: ${gadoga_setting}!"
                exit 1
            fi
            if [ ! -f ${gzvcf} ]
            then
	        
		# Run gatks's genotypegvcf program.
	        rm -f ${out}
	        rm -f ${log}
	        sbatch --account=nn9244k -o ${out} run_gatk2.slurm ${chromosome_id} ${ref} ${gzvcf} ${log} ${gzgvcf_string}
            fi
        done
    done
    
done
