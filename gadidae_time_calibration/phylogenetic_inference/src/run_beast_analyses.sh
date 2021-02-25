# m_matschiner Thu Feb 25 20:25:15 CET 2021

# Launch beast analyses in each replicate directory.
for rep_dir in ../res/beast/replicates/r??
do
	cd ${rep_dir}
	sbatch run_beast.slurm
	cd -
done