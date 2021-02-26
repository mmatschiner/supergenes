# m_matschiner Thu Aug 30 13:11:32 CEST 2018

# Make the output directory.
mkdir -p ../res/mapping

# Make the log directory.
mkdir -p ../log/misc

# Set the reference.
ref=../data/assemblies/gadMor2.fasta

# Map each pair of fastq files against the gadmor reference.
for i in ../data/fastq/*/*_R1.fastq.gz
do
    specimen_id=`basename ${i} | cut -d "_" -f 1`
    library_id=`basename ${i} | cut -d "_" -f 2`
    fastq1=${i}
    fastq2=${i%_R1.fastq.gz}_R2.fastq.gz
    bam=../res/bam/${specimen_id}_${library_id}.bam
    out=../log/misc/map.${specimen_id}_${library_id}.out
    log=../log/misc/map.${specimen_id}_${library_id}.log
    rm -f ${out}
    sbatch -o ${out} map.slurm ${ref} ${bam} ${log} ${fastq1} ${fastq2}
done