# m_matschiner Sun Jul 15 23:56:44 CEST 2018

# Make the output directory.
mkdir -p ../res/fasta/converted_from_bam

# Set the reference.
ref=../data/assemblies/gadMor2/gadMor2.fasta

# Convert bam files to fasta consensus sequences.
for assembly_id in gadMor2 gadMor_Stat melAeg
do
	bam=../res/bam/${assembly_id}.bam
	out=../log/misc/bam2fasta.${assembly_id}.out
	fasta=../res/fasta/converted_from_bam/${assembly_id}.fasta
	sbatch -o ${out} convert_bam_to_fasta.slurm ${bam} ${ref} ${fasta}
done
