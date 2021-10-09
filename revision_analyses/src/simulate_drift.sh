# m_matschiner Sat Oct 9 11:53:25 CEST 2021

# Set the account.
acct=nn9883k

# Set the number of simulations.
n_simulations=10000

# Simulate drift for different population sizes.
for pop_size in 1000 3000 10000 30000 100000
do
    # Set the output file.
    out=../res/tables/drift_${pop_size}.txt

    # Set the log file.
    log=../log/misc/drift_${pop_size}.txt

    # Simulate drift.
    sbatch --account=${acct} -o ${log} simulate_drift.slurm ${pop_size} ${n_simulations} ${out}
done
