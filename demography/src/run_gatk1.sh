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

# Run gatk's haplotypecaller.
for bam in ../data/bam/??????_????.bam ../data/bam/Gadcha.bam
do
    for lg_id in LG0{1..9} LG{10..23}
    do
        bam_id=`basename ${bam%.bam}`
        out=`readlink -f ../log/gatk/run_gatk1.${bam_id}.${lg_id}.out`
        gzgvcf=`readlink -f ../res/gatk/${bam_id}.${lg_id}.g.vcf.gz`
        log=`readlink -f ../log/gatk/run_gatk1.${bam_id}.${lg_id}.log`
        if [ ! -f ${gzgvcf} ]
        then
            rm -f ${out}
            rm -f ${log}
            sbatch -o ${out} run_gatk1.slurm ../data/assemblies/gadMor2.fasta ${bam} ${gzgvcf} ${log} ${lg_id}
            sleep_while_too_busy
        fi
    done
done
