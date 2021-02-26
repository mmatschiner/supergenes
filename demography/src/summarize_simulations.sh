# m_matschiner Fri Jul 10 13:56:56 CEST 2020

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0

# Summarize simulations with script.
ruby summarize_simulations.rb ../res/msprime/p_bottleneck.txt ../res/msprime/pi.txt > ../res/msprime/summary.txt
