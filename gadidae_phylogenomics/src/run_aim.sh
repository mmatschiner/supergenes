# m_matschiner Tue Oct 16 01:01:23 CEST 2018

# Run aim.
rep_dir=../res/aim/age_prior/replicates
for dir in ${rep_dir}/r*
do
    cd ${dir}
    sbatch run_aim.slurm
    cd -
done
