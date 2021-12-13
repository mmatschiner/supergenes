# m_matschiner Tue Dec 8 11:39:12 CET 2020

# Load modules.
module load R/3.6.2-foss-2019b

# Make the plots directory.
mkdir -p ../res/plots

# Plot the repeat content for the genome-wide background and the four high-ld regions.
Rscript plot_repeat_content.r
