# m_matschiner Mon Sep 24 13:29:05 CEST 2018

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

# Run run_gatk2.slurm with all vcf files of chromosomes.
for i in `seq -w 23`
do
    chromosome_id=LG${i}
    gzvcf="../res/gatk/${chromosome_id}.vcf.gz"
    log="../log/gatk/run_gatk2.${chromosome_id}.log"
    out="../log/gatk/run_gatk2.${chromosome_id}.out"
    gvcf_string=`for i in ../res/gatk/??????.${chromosome_id}.g.vcf.gz; do echo -n "${i%.gz} "; done`
    if [ ! -f ${gzvcf} ]
    then
        # Make sure that vcf files have been uncompressed.
        for i in ${gvcf_string}
        do
        if [ ! -f ${i} ]
        then
            file_id=`basename ${i%.g.vcf}`
            echo "Starting uncompression of file ${i}.gz."
            sbatch -o "../log/prepare/${file_id}.out" uncompress_and_index.slurm ${i}.gz ${chromosome_id} ../data/assemblies/gadMor2.fasta
        fi
    done
    sleep_while_busy

    # Make sure index files have been created for each vcf file.
    n_vcf_files=`ls ../res/gatk/??????.${chromosome_id}.g.vcf | wc -l`
    n_idx_files=`ls ../res/gatk/??????.${chromosome_id}.g.vcf.idx | wc -l`
    if [ $n_vcf_files -gt $n_idx_files ]
    then
        echo "ERROR: Not all index files have been created (currently present: ${n_idx_files})!"
        exit 1
    fi

    # Run gatks's genotypegvcf program.
    sbatch -o ${out} run_gatk2.slurm ${chromosome_id} ../data/assemblies/gadMor2.fasta ${gzvcf} ${log} ${gvcf_string}
    fi
    
done
