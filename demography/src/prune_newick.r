# m_matschiner Mon Sep 23 11:05:55 CEST 2019

# Load libraries.
library(phytools)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
intree_file_name <- args[1]
keep_list_file_name <- args[2]
outtree_file_name <- args[3]

# Read the tree.
tree <- read.newick(intree_file_name)
tip_list <- tree$tip.label

# Read the list of taxa to keep.
keep_list <- read.table(keep_list_file_name)
keep_list <- as.vector(keep_list$V1)
drop_list <- setdiff(tip_list,keep_list)

# Remove all tips not in the keep list.
pruned_tree <- drop.tip(tree, tip=drop_list)

# Write the pruned trees to file.
write.tree(pruned_tree, file=outtree_file_name)