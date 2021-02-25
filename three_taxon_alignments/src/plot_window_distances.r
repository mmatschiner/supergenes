# m_matschiner Wed Jul 4 09:31:47 CEST 2018

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
table_file_name <- args[1]
lg_id <- args[2]
distance_plot_file_name <- args[3]
relative_distance_plot_file_name <- args[4]

# Read the table.
t <- read.table(table_file_name, header=T)

# Plot the absolute distances.
pdf(distance_plot_file_name, height=7, width=7)
plot(t$window_center, t$distance_between_melAeg_and_first/100000, type="l", xlab="Chromosome position", ylab="Distance", ylim=c(0,0.09), main=lg_id, col="#657b83")
lines(t$window_center, t$distance_between_gadMor_Stat_and_first/100000, col="#073642")
dev.off()

# Plot the relative distances.
pdf(relative_distance_plot_file_name, height=7, width=7)
plot(t$window_center, t$distance_between_gadMor_Stat_and_first/t$distance_between_melAeg_and_first, type="l", xlab="Chromosome position", ylab="Relative distance", ylim=c(0,1), main=lg_id, col="#073642")
dev.off()