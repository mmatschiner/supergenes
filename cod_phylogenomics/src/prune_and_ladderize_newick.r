# m_matschiner Mon Sep 23 11:05:55 CEST 2019

# Load libraries.
library(phytools)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
intrees_file_name <- args[1]
keep_list_file_name <- args[2]
outtrees_file_name <- args[3]

# Read the tree.
trees <- read.newick(intrees_file_name)
if (class(trees) == "multiPhylo") {
   tip_list <- trees[[1]]$tip.label
} else {
   tip_list <- trees$tip.label
}

# Read the list of taxa to keep.
keep_list <- read.table(keep_list_file_name)
keep_list <- as.vector(keep_list$V1)
drop_list <- setdiff(tip_list,keep_list)

# Remove all tips not in the keep list.
if (class(trees) == "multiPhylo") {
   pruned_trees <- lapply(trees, drop.tip, tip=drop_list)
   class(pruned_trees) <- "multiPhylo"
} else {
   pruned_trees <- drop.tip(trees, tip=drop_list)
}

# Ladderize the tree.
if (class(trees) == "multiPhylo") {
   ladderized_trees <- lapply(pruned_trees, ladderize)
   class(ladderized_trees) <- "multiPhylo"
} else {
   ladderized_trees <- ladderize(pruned_trees)
}

# Write the pruned trees to file.
write.tree(ladderized_trees, file=outtrees_file_name)