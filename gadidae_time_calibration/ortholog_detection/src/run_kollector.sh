# m_matschiner Wed Jul 18 15:48:49 CEST 2018

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/kollector

# Make the log directory if it doesn't exist yet.
mkdir -p ../log/kollector

# Set the query file.
query=../data/queries/gadmor_exons_nucl.fasta

# Start a slurm script per specimen.
for sample in gadmac gadoga
do
    out=../log/kollector/${sample}.out
    log=../log/kollector/${sample}.log
    res=../res/kollector/${sample}.fasta
    if [ ! -d ${res} ]
    then
        rm -f ${out}
        rm -f ${log}
        echo -n "Submitting job for sample ${sample}..."
        sbatch -o ${out} run_kollector.slurm ${sample} ${query} ${res} ${log} &> /dev/null
        echo " done."
    fi
done