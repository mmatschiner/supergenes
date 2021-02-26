# m_matschiner Mon Sep 24 13:29:05 CEST 2018

# Run run_gatk2.slurm with vcf files of all chromosomes.
for i in `seq -w 23`
do
    lg_id=LG${i}
    gzvcf="../res/gatk/${lg_id}.vcf.gz"
    log="../log/gatk/run_gatk2.${lg_id}.log"
    out="../log/gatk/run_gatk2.${lg_id}.out"
    gzgvcf_string=`for i in ../res/gatk/Gadcha.${lg_id}.g.vcf.gz ../res/gatk/??????_????.${lg_id}.g.vcf.gz; do echo -n "${i} "; done`
    if [ ! -f ${gzvcf} ]
    then
        # Run gatks's genotypegvcf program.
        rm -f ${out}
        rm -f ${log}
        sbatch -o ${out} run_gatk2.slurm ${lg_id} ../data/assemblies/gadMor2.fasta ${gzvcf} ${log} ${gzgvcf_string}
    fi
done
