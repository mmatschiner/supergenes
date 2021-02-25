# m_matschiner Wed Jul 4 09:41:04 CEST 2018

# Load the r module.
module load R/3.4.4

# Use an r script to plot the pairwise genetic distances.
for i in ../res/alignments/merged/lg??.threaded.refined.finished.txt
do
	lg_id=`basename ${i%.threaded.refined.finished.txt}`
	distance_plot=${i%.txt}.dists.pdf
	relative_distance_plot=${i%.txt}.rel_dists.pdf
	Rscript plot_window_distances.r ${i} ${lg_id} ${distance_plot} ${relative_distance_plot}
done
