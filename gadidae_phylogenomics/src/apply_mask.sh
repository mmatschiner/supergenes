# m_matschiner Thu Oct 4 12:42:41 CEST 2018

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    while [ $n_jobs -gt 300 ]
    do
        sleep 10
        n_jobs=`squeue -u michaelm | wc -l`
    done
}

# Copy masks.
mkdir -p ../data/masks
cp ../../three_taxon_alignments/results/masks/*.bed ../data/masks

# Generate filtered versions of all vcf files.
for vcfgz in ../res/gatk/LG??.vcf.gz
do
    lg=`basename ${vcfgz%.vcf.gz}`
    mask=../data/masks/lg${lg:2:3}.bed
    masked_vcfgz=${vcfgz%.vcf.gz}.masked.vcf.gz
    if [ ! -f ${masked_vcfgz} ]
    then
        out=../log/gatk/mask.${lg}.out
        log=../log/gatk/mask.${lg}.log
        # Remove log files if they exist already.
        rm -f ${out}
        rm -f ${log}
        sbatch -o ${out} apply_mask.slurm ${vcfgz} ../data/assemblies/gadMor2.fasta ${mask} ${masked_vcfgz} ${log}
        sleep_while_too_busy
    fi
done