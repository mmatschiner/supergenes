# m_matschiner Sat Mar 18 15:38:05 CET 2017

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/orthologs/03

# Specify the input and output directories.
indir=../res/alignments/orthologs/02/
outdir=../res/alignments/orthologs/03/
threshold=0.3

# Filter sequences by dNdS.
out=../log/misc/filter_sequences_by_dNdS.out
rm -f ${out}
sbatch -o ${out} filter_sequences_by_dNdS.slurm ${indir} ${outdir} ${threshold}
