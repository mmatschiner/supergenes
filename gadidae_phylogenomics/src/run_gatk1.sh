# m_matschiner Wed Sep 19 23:09:52 CEST 2018

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    while [ $n_jobs -gt 350 ]
    do
        sleep 10
        n_jobs=`squeue -u michaelm | wc -l`
    done
}

# Make the output directories.
mkdir -p ../res/gatk
mkdir -p ../log/gatk

# Run gatk's haplotypecaller for each specimen and chromosome separately.
for specimen_id in `ls ../res/mapping/??????.bam | sed 's/.bam//g' | cut -d "/" -f 4`
do
    bam_with_relative_path="../res/mapping/${specimen_id}.bam"
    bam=`basename ${bam_with_relative_path}`
    truncated_bam=${bam%.bam}
    for chromosome_id in LG0{1..9} LG{10..23}
    do
        out_with_relative_path="../log/gatk/run_gatk1.${truncated_bam}.${chromosome_id}.out"
        gvcf_with_relative_path="../res/gatk/${truncated_bam}.${chromosome_id}.g.vcf.gz"
        log_with_relative_path="../log/gatk/run_gatk1.${truncated_bam}.${chromosome_id}.log"
        if [ ! -f ${gvcf_with_relative_path} ]
        then
            sbatch -o ${out_with_relative_path} run_gatk1.slurm ../data/assemblies/gadMor2.fasta ${bam_with_relative_path} ${gvcf_with_relative_path} ${log_with_relative_path} ${chromosome_id}
            sleep_while_too_busy
        fi
    done
done
