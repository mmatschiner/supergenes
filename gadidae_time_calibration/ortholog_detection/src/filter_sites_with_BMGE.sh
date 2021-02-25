# m_matschiner Sat Mar 18 17:28:45 CET 2017

# Load the ruby module
module load ruby/2.1.5

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/orthologs/04

# Filter sites for missing data and entropy with BMGE.
ruby filter_sites_with_BMGE.rb ../res/alignments/orthologs/03/ ../res/alignments/orthologs/04/ ../bin/BMGE.jar 0.2 0.5