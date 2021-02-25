# m_matschiner Thu Feb 25 20:22:42 CET 2021

# Prepare five replicate directories for beast analyses.
for n in `seq 5`
do
	rep_dir=../res/replicates/r0${n}
	mkdir -p ${rep_dir}
	cp run_beast.slurm ${rep_dir}
	cp ../data/xml/91genes.xml ${rep_dir}
done