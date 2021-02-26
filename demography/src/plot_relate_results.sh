# m_matschiner Mon Nov 2 19:46:43 CET 2020

# Load modules.
module purge
module load R/4.0.0-foss-2020a

# Set the generation time.
gen_time=10

# Plot population sizes over time estimated with Relate.
for coal in ../res/relate/output/simulated*/*.pairwise.coal
do
	pdf=${coal%.coal}.pdf
	Rscript plot_relate.r ${coal} ${gen_time} ${pdf}
done
