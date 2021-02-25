# m_matschiner Tue Sep 1 10:45:06 CEST 2020

# Load libraries.
library(ggplot2)

# Read the input file.
t <- read.table("../res/tables/mutational_load.txt", header=T)
proportions_of_stop_codon_plot_file_name = "../res/plots/proportions_of_stop_codons.pdf"
proportions_of_changed_aas_in_gadMor2_plot_file_name = "../res/plots/proportions_of_changed_aas_in_gadMor2.pdf"
proportions_of_changed_aas_in_gadMor2_Stat_plot_file_name = "../res/plots/proportions_of_changed_aas_in_gadMor2_Stat.pdf"

# Add column names.
colnames(t) <- c("lg","from","to","n_aas_gadMor2","n_aas_gadMor_Stat","n_aas_melAeg","n_stops_gadMor2","n_stops_gadMor_Stat","n_stops_melAeg","n_aas_unchanged_gadMor2_gadMor_Stat","n_aas_changed_gadMor2_gadMor_Stat","n_aas_unchanged_gadMor2_melAeg","n_aas_changed_gadMor2_melAeg","n_aas_unchanged_gadMor_Stat_melAeg","n_aas_changed_gadMor_Stat_melAeg")
proportion_of_stop_aas_colinear <- c()
proportion_of_stop_aas_lg01 <- c()
proportion_of_stop_aas_lg02 <- c()
proportion_of_stop_aas_lg07 <- c()
proportion_of_stop_aas_lg12 <- c()
proportion_of_changed_aas_in_gadMor2_colinear <- c()
proportion_of_changed_aas_in_gadMor2_lg01 <- c()
proportion_of_changed_aas_in_gadMor2_lg02 <- c()
proportion_of_changed_aas_in_gadMor2_lg07 <- c()
proportion_of_changed_aas_in_gadMor2_lg12 <- c()
proportion_of_changed_aas_in_gadMor2_Stat_colinear <- c()
proportion_of_changed_aas_in_gadMor2_Stat_lg01 <- c()
proportion_of_changed_aas_in_gadMor2_Stat_lg02 <- c()
proportion_of_changed_aas_in_gadMor2_Stat_lg07 <- c()
proportion_of_changed_aas_in_gadMor2_Stat_lg12 <- c()

# Add a column for the window centers.
t$center = (t$from + t$to)/2.0

# Add a column for the region.
regions <- c()

# Collect values for the proportion of stop codons, and for the proportions of amino acids changed per window.
# The first statistic is collected only for the gadMor2_Stat cod assembly because the gadMor2 should not have any stop codons given that the annotation was made for it.
# The second statistic is collected for both the gadMor2_Stat cod assembly and the gadMor2 gadMor2 assembly.
for(x in 1:length(t$lg)) {
	if( t$lg[x] == "LG01" && t$center[x] >= 9114741 &&  t$center[x] <= 26192489 ) {
		regions <- append(regions, "inversion_lg01")
		proportion_of_stop_aas_lg01 <- append(proportion_of_stop_aas_lg01, t$n_stops_gadMor_Stat[x]/t$n_aas_gadMor_Stat[x])
		proportion_of_changed_aas_in_gadMor2_lg01 <- append(proportion_of_changed_aas_in_gadMor2_lg01, t$n_aas_changed_gadMor2_melAeg[x]/(t$n_aas_changed_gadMor2_melAeg[x]+t$n_aas_unchanged_gadMor2_melAeg[x]))
		proportion_of_changed_aas_in_gadMor2_Stat_lg01 <- append(proportion_of_changed_aas_in_gadMor2_Stat_lg01, t$n_aas_changed_gadMor_Stat_melAeg[x]/(t$n_aas_changed_gadMor_Stat_melAeg[x]+t$n_aas_unchanged_gadMor_Stat_melAeg[x]))
	} else if( t$lg[x] == "LG02" && t$center[x] >= 18489307 &&  t$center[x] <= 24050282 ) {
		regions <- append(regions, "inversion_lg02")
		proportion_of_stop_aas_lg02 <- append(proportion_of_stop_aas_lg02, t$n_stops_gadMor_Stat[x]/t$n_aas_gadMor_Stat[x])
		proportion_of_changed_aas_in_gadMor2_lg02 <- append(proportion_of_changed_aas_in_gadMor2_lg02, t$n_aas_changed_gadMor2_melAeg[x]/(t$n_aas_changed_gadMor2_melAeg[x]+t$n_aas_unchanged_gadMor2_melAeg[x]))
		proportion_of_changed_aas_in_gadMor2_Stat_lg02 <- append(proportion_of_changed_aas_in_gadMor2_Stat_lg02, t$n_aas_changed_gadMor_Stat_melAeg[x]/(t$n_aas_changed_gadMor_Stat_melAeg[x]+t$n_aas_unchanged_gadMor_Stat_melAeg[x]))
	} else if( t$lg[x] == "LG07" && t$center[x] >= 13606502 &&  t$center[x] <= 23016726 ) {
		regions <- append(regions, "inversion_lg07")
		proportion_of_stop_aas_lg07 <- append(proportion_of_stop_aas_lg07, t$n_stops_gadMor_Stat[x]/t$n_aas_gadMor_Stat[x])
		proportion_of_changed_aas_in_gadMor2_lg07 <- append(proportion_of_changed_aas_in_gadMor2_lg07, t$n_aas_changed_gadMor2_melAeg[x]/(t$n_aas_changed_gadMor2_melAeg[x]+t$n_aas_unchanged_gadMor2_melAeg[x]))
		proportion_of_changed_aas_in_gadMor2_Stat_lg07 <- append(proportion_of_changed_aas_in_gadMor2_Stat_lg07, t$n_aas_changed_gadMor_Stat_melAeg[x]/(t$n_aas_changed_gadMor_Stat_melAeg[x]+t$n_aas_unchanged_gadMor_Stat_melAeg[x]))
	} else if( t$lg[x] == "LG12" && t$center[x] >= 589105 &&  t$center[x] <= 13631347 ) {
		regions <- append(regions, "inversion_lg12")
		proportion_of_stop_aas_lg12 <- append(proportion_of_stop_aas_lg12, t$n_stops_gadMor_Stat[x]/t$n_aas_gadMor_Stat[x])
		proportion_of_changed_aas_in_gadMor2_lg12 <- append(proportion_of_changed_aas_in_gadMor2_lg12, t$n_aas_changed_gadMor2_melAeg[x]/(t$n_aas_changed_gadMor2_melAeg[x]+t$n_aas_unchanged_gadMor2_melAeg[x]))
		proportion_of_changed_aas_in_gadMor2_Stat_lg12 <- append(proportion_of_changed_aas_in_gadMor2_Stat_lg12, t$n_aas_changed_gadMor_Stat_melAeg[x]/(t$n_aas_changed_gadMor_Stat_melAeg[x]+t$n_aas_unchanged_gadMor_Stat_melAeg[x]))
	} else {
		regions <- append(regions, "colinear")
		proportion_of_stop_aas_colinear <- append(proportion_of_stop_aas_colinear, t$n_stops_gadMor_Stat[x]/t$n_aas_gadMor_Stat[x])
		proportion_of_changed_aas_in_gadMor2_colinear <- append(proportion_of_changed_aas_in_gadMor2_colinear, t$n_aas_changed_gadMor2_melAeg[x]/(t$n_aas_changed_gadMor2_melAeg[x]+t$n_aas_unchanged_gadMor2_melAeg[x]))
		proportion_of_changed_aas_in_gadMor2_Stat_colinear <- append(proportion_of_changed_aas_in_gadMor2_Stat_colinear, t$n_aas_changed_gadMor_Stat_melAeg[x]/(t$n_aas_changed_gadMor_Stat_melAeg[x]+t$n_aas_unchanged_gadMor_Stat_melAeg[x]))
	}
}
t$region = regions

# Plot the proportions of stop codons.
pdf(proportions_of_stop_codon_plot_file_name, height=7, width=7)
ggplot(data = t, mapping = aes(x = region, y = t$n_stops_gadMor_Stat/t$n_aas_gadMor_Stat, color=region)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic(base_size = 12)
dev.off()

# Plot the proportions of changed amino acids in the gadMor2 genome compared to the melAeg.
pdf(proportions_of_changed_aas_in_gadMor2_plot_file_name, height=7, width=7)
ggplot(data = t, mapping = aes(x = region, y = t$n_aas_changed_gadMor2_melAeg/(t$n_aas_changed_gadMor2_melAeg+t$n_aas_unchanged_gadMor2_melAeg), color=region)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic(base_size = 12)
dev.off()

# Plot the proportions of changed amino acids in the gadMor2_Stat genome compared to the melAeg.
pdf(proportions_of_changed_aas_in_gadMor2_Stat_plot_file_name, height=7, width=7)
ggplot(data = t, mapping = aes(x = region, y = t$n_aas_changed_gadMor_Stat_melAeg/(t$n_aas_changed_gadMor_Stat_melAeg+t$n_aas_unchanged_gadMor_Stat_melAeg), color=region)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic(base_size = 12)
dev.off()


# Test if the proportions of stop codons are greater in inversion regions.
t.test(proportion_of_stop_aas_lg01, proportion_of_stop_aas_colinear, "greater")
t.test(proportion_of_stop_aas_lg02, proportion_of_stop_aas_colinear, "greater")
t.test(proportion_of_stop_aas_lg07, proportion_of_stop_aas_colinear, "greater")
t.test(proportion_of_stop_aas_lg12, proportion_of_stop_aas_colinear, "greater")

# Test if the proportions of changed amino acids are greater in inversion regions in the gadMor2 genome.
t.test(proportion_of_changed_aas_in_gadMor2_lg01, proportion_of_changed_aas_in_gadMor2_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_lg02, proportion_of_changed_aas_in_gadMor2_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_lg07, proportion_of_changed_aas_in_gadMor2_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_lg12, proportion_of_changed_aas_in_gadMor2_colinear, "greater")

# Test if the proportions of changed amino acids are greater in inversion regions in the gadMor2_Stat genome.
t.test(proportion_of_changed_aas_in_gadMor2_Stat_lg01, proportion_of_changed_aas_in_gadMor2_Stat_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_Stat_lg02, proportion_of_changed_aas_in_gadMor2_Stat_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_Stat_lg07, proportion_of_changed_aas_in_gadMor2_Stat_colinear, "greater")
t.test(proportion_of_changed_aas_in_gadMor2_Stat_lg12, proportion_of_changed_aas_in_gadMor2_Stat_colinear, "greater")
