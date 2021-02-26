# m_matschiner Thu Dec 5 13:58:57 CET 2019

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    while [ $n_jobs -gt 300 ]
    do
        sleep 10
        n_jobs=`squeue -u michaelm | wc -l`
    done
}

# Launch all replicate snapp analyses.
for snapp_dir in ../res/snapp/regions/*/replicates/r??
do
    cd ${snapp_dir}
    sbatch run_snapp.slurm
    cd -
    sleep_while_too_busy
done
