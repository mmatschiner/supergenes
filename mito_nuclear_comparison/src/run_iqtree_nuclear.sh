# m_matschiner Tue Nov 5 14:18:45 CET 2019

# Load modules.
module load iqtree/1.7-beta7

# Make the resul directory.
mkdir -p ../res/iqtree/nuclear

# Run iqtree for the concatenated nuclear exon alignment.
iqtree -s ../res/alignments/nuclear/concatenated.phy -B 1000

# Move the iqtree results into the result directory.
mv ../res/alignments/nuclear/concatenated.phy.* ../res/iqtree/nuclear