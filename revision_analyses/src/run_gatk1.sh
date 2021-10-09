# michaelm Mon Feb 8 00:40:11 CET 2021

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    while [ $n_jobs -gt 400 ]
    do
        sleep 10
        n_jobs=`squeue -u michaelm | wc -l`
    done
}

# Set the account.
acct=nn9244k

# Repeat gatk1 for each of the alternative references.
for ref in ../data/assemblies/gadMor3.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Make the output directories.
    mkdir -p ../res/gatk/${ref_id}
    mkdir -p ../log/gatk/

    # Run gatk's haplotypecaller.
    for bam in ../res/mapping/${ref_id}/??????_????.bam ../res/mapping/${ref_id}/??????.bam
    do
        for lg_id in chr{1..23}
        do
            bam_id=`basename ${bam%.bam}`
            out=`readlink -f ../log/gatk/run_gatk1.${ref_id}_${bam_id}.${lg_id}.out`
            gzgvcf=`readlink -f ../res/gatk/${ref_id}/${bam_id}.${lg_id}.g.vcf.gz`
            log=`readlink -f ../log/gatk/run_gatk1.${ref_id}_${bam_id}.${lg_id}.log`
            if [ ! -f ${gzgvcf} ]
            then
                rm -f ${out}
                rm -f ${log}
                sbatch --account=${acct} -o ${out} run_gatk1.slurm ${ref} ${bam} ${gzgvcf} ${log} ${lg_id}
                sleep_while_too_busy
            fi
        done
    done
done
