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

# Load modules.
module load R/3.5.1-foss-2018b

# Set parameters.
min_min_ess=100
max_n_samples=500000

# Launch all replicate snapp analyses.
for snapp_dir in ../res/snapp/windows/LG??_0????????_0????????/replicates/r0?
do
    window_id=`echo ${snapp_dir} | cut -d "/" -f 5`
    log=${snapp_dir}/${window_id}.log
    if [ -f ${log} ]
    then
        n_samples=`tail -n 1 ${log} | cut -f 1`
        if (( $(echo "${n_samples} < ${max_n_samples}" | bc -l ) ))
        then
            min_ess=`Rscript get_min_ess.r ${log}`
            if (( $( echo "${min_ess} < ${min_min_ess}" | bc -l ) ))
            then
                cd ${snapp_dir}
                echo "Resuming SNAPP analysis in ${snapp_dir} (min. ESS = ${min_ess})."
                sbatch run_snapp.slurm
                cd -
                sleep_while_too_busy
            else
                echo "SNAPP analysis in ${snapp_dir} has completed (min. ESS = ${min_ess})."
            fi
        else
            echo "Skipping SNAPP analysis in ${snapp_dir} (n. samples = ${n_samples})."
        fi
    else
        cd ${snapp_dir}
        echo "Starting SNAPP analysis in ${snapp_dir}."
        sbatch run_snapp.slurm
        cd - &> /dev/null
        sleep_while_too_busy
    fi
done
