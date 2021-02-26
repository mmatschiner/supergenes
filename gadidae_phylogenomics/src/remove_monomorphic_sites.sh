# m_matschiner Tue Oct 29 14:22:04 CET 2019

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    while [ $n_jobs -gt 300 ]
    do
        sleep 10
        n_jobs=`squeue -u michaelm | wc -l`
    done
}

# Generate filtered versions of all vcf files.
for vcfgz in ../res/gatk/LG??.masked.vcf.gz
do
    lg=`basename ${vcfgz%.vcf.gz}`
    variable_vcfgz=${vcfgz%.vcf.gz}.variable.vcf.gz
    if [ ! -f ${variable_vcfgz} ]
    then
        out=../log/gatk/remove_monomorphic.${lg}.out
        log=../log/gatk/remove_monomorphic.${lg}.log
        # Remove log files if they exist already.
        rm -f ${out}
        rm -f ${log}
        sbatch -o ${out} remove_monomorphic_sites.slurm ${vcfgz} ../data/assemblies/gadMor2.fasta ${variable_vcfgz} ${log}
        sleep_while_too_busy
    fi
done