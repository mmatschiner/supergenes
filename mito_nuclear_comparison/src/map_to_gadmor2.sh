# m_matschiner Wed Oct 23 15:18:28 CEST 2019

# Load modules.
module load bwa/0.7.10
module load samtools/1.9
module load picard/2.6.0

# Make the log directory.
mkdir -p ../log/mapping

# Make the results directory.
mkdir -p ../res/mapping/bam

# Set the reference.
ref=../data/ref/gadMor2.fasta

# Prepare the reference.
if [ ! -f ${ref}.fai ]
then
    bwa index ${ref} &> ../log/mapping/prepare.log
    samtools faidx ${ref}
    java -jar picard/2.6.0/picard.jar CreateSequenceDictionary \
        R=${ref} \
        O=${ref%.fasta}.dict \
        URI=./${ref} &> ../log/mapping/prepare.log
fi
if [ ! -f ${ref%.fasta}.dict ]
then
    java -jar picard/2.6.0/picard.jar CreateSequenceDictionary \
        R=${ref} \
        O=${ref%.fasta}.dict \
        URI=./${ref} &> ../log/mapping/prepare.log
fi

# Map all paired fastq files to the reference.
for fastq1 in ../data/reads/*R1.fastq.gz
do
    species_id=`basename ${fastq1} | cut -d "_" -f 1`
    specimen_id=`basename ${fastq1} | cut -d "_" -f 2`
    fastq2=${fastq1%R1.fastq.gz}R2.fastq.gz
    bam=../res/mapping/bam/${species_id}_${specimen_id}.bam
    log=../log/mapping/${species_id}_${specimen_id}.log
    out=../log/mapping/${species_id}_${specimen_id}.out
    if [ ! -f ${bam} ]
    then
        rm -f ${out}
        sbatch -o ${out} --account=nn9244k map.slurm ${fastq1} ${fastq2} ${ref} ${bam} ${log}
    fi
done