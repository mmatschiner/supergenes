# m_matschiner Tue Sep 1 10:45:06 CEST 2020

# Load libraries.
library(ggplot2)

# Read the input file.
t <- read.table("../res/tables/repeat_content.txt", header=T)
repeat_content_plot_file_name = "../res/plots/repeat_content.pdf"

# Add column names.
colnames(t) <- c("lg","from","to","length","n_repeat_sites","repeat_content")
repeat_contents_colinear <- c()
repeat_contents_inversion_lg01 <- c()
repeat_contents_inversion_lg02 <- c()
repeat_contents_inversion_lg07 <- c()
repeat_contents_inversion_lg12 <- c()

# Add a column for the window centers.
t$center = (t$from + t$to)/2.0

# Add a column for the region.
regions <- c()
for(x in 1:length(t$lg)) {
    if( t$lg[x] == "LG01" && t$center[x] >= 9114741 &&  t$center[x] <= 26192489 ) {
        regions <- append(regions, "inversion_lg01")
	repeat_contents_inversion_lg01 <- append(repeat_contents_inversion_lg01, t$repeat_content[x])
    } else if( t$lg[x] == "LG02" && t$center[x] >= 18489307 &&  t$center[x] <= 24050282 ) {
        regions <- append(regions, "inversion_lg02")
        repeat_contents_inversion_lg02 <- append(repeat_contents_inversion_lg02, t$repeat_content[x])
    } else if( t$lg[x] == "LG07" && t$center[x] >= 13606502 &&  t$center[x] <= 23016726 ) {
        regions <- append(regions, "inversion_lg07")
        repeat_contents_inversion_lg07 <- append(repeat_contents_inversion_lg07, t$repeat_content[x])
    } else if( t$lg[x] == "LG12" && t$center[x] >= 589105 &&  t$center[x] <= 13631347 ) {
        regions <- append(regions, "inversion_lg12")
        repeat_contents_inversion_lg12 <- append(repeat_contents_inversion_lg12, t$repeat_content[x])
    } else {
        regions <- append(regions, "colinear")
	repeat_contents_colinear <- append(repeat_contents_colinear, t$repeat_content[x])
    }
}
t$region = regions

# Plot absolute measures.
pdf(repeat_content_plot_file_name, height=7, width=7)
ggplot(data = t, mapping = aes(x = region, y = repeat_content, color=region)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic(base_size = 12)
dev.off()

# Test if repeat content is greater in inversion regions.
cat("Sample sizes: ", length(repeat_contents_colinear), "/", length(repeat_contents_inversion_lg01), "\n")
t.test(repeat_contents_inversion_lg01, repeat_contents_colinear, "greater")
cat("Sample sizes: ", length(repeat_contents_colinear), "/", length(repeat_contents_inversion_lg02), "\n")
t.test(repeat_contents_inversion_lg02, repeat_contents_colinear, "greater")
cat("Sample sizes: ", length(repeat_contents_colinear), "/", length(repeat_contents_inversion_lg07), "\n")
t.test(repeat_contents_inversion_lg07, repeat_contents_colinear, "greater")
cat("Sample sizes: ", length(repeat_contents_colinear), "/", length(repeat_contents_inversion_lg12), "\n")
t.test(repeat_contents_inversion_lg12, repeat_contents_colinear, "greater")