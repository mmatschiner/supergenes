# m_matschiner Wed May 23 19:16:56 CEST 2018

# Load the dplyr library.
library(dplyr)
library(ggplot2)

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
table_file_name <- args[1]
plot_file_name <- args[2]

# Read the input files.
t <- read.table(table_file_name, header=T, stringsAsFactors=F)

# Get the number of chromosomes.
n_chr <- length(unique(t$chr))

# Add columns for per-lg and cumulative positions.
t$pos <- t$windowStart + 0.5 * (t$windowEnd-t$windowStart)
t$cumulative_pos <- NA
s <- 0
nbp <- c()
for (i in unique(t$chr)){
    nbp[i] <- max(t[t$chr == i,]$pos)
    t[t$chr == i,"cumulative_pos"] <- t[t$chr == i,"pos"] + s
    s <- s + nbp[i]
}

# Get the central position of each lg.
axis.set <- t %>% 
  	 group_by(chr) %>% 
  	 summarize(center = (max(cumulative_pos) + min(cumulative_pos)) / 2)

manhplot <- ggplot(t, aes(x=cumulative_pos, y=D)) +
	 geom_point(aes(color=as.factor(chr)), alpha = 0.75, size = 1.25) +
 	 scale_color_manual(values = rep(c("darkslateblue","cadetblue"), n_chr)) +
  	 scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  	 scale_y_continuous(expand = c(0,0), limits = c(-1,1)) +
  	 geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") + 
  	 labs(x = NULL, y = "D") + 
  	 theme_minimal() +
  	 theme( 
    	 	legend.position = "none",
    	 	panel.border = element_blank(),
    	 	panel.grid.major.x = element_blank(),
    	 	panel.grid.minor.x = element_blank(),
    		axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
 	)
ggsave(plot_file_name, manhplot, width = 30, height = 20, units = "cm")

