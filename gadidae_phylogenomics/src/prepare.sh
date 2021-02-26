# m_matschiner Thu Aug 30 13:13:31 CEST 2018

# Generate dictionary and index files for the reference.
out=../log/misc/prepare.out
log=../log/misc/prepare.log
sbatch -o ${out} prepare.slurm ../data/assemblies/gadMor2.fasta ${log}