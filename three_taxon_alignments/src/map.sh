# m_matschiner Sun Jul 15 21:01:19 CEST 2018

# Make the output directory.
mkdir -p ../res/bam

# Set the reference.
ref=../data/assemblies/gadMor2/gadMor2.fasta

# Map fastq files.
for assembly_id in gadMor2 gadMor_Stat melAeg
do
	bam=../res/bam/${assembly_id}.bam
	out=../log/misc/map.${assembly_id}.out
	log=../log/misc/map.${assembly_id}.log
	fastq1=../data/fastq/${assembly_id}_R1.fastq.gz
	fastq2=../data/fastq/${assembly_id}_R2.fastq.gz
	sbatch -o ${out} map.slurm ${ref} ${bam} ${log} ${fastq1} ${fastq2}
done