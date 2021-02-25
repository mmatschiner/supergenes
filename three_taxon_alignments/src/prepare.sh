# m_matschiner Sun Jul 15 20:35:56 CEST 2018

# Generate dictionary and index files for the reference.
out=../log/misc/prepare.out
log=../log/misc/prepare.log
sbatch -o ${out} prepare.slurm ../data/assemblies/gadMor2/gadMor2.fasta ${log}