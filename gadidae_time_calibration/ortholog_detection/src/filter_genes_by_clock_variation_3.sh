# m_matschiner Tue Jul 24 12:58:06 CEST 2018

# Load the ruby module.
module load Ruby/2.7.1-GCCcore-8.3.0

# Read beast log files and filter according to coefficent of variation.
ruby filter_genes_by_clock_variation_3.rb ../res/alignments/orthologs/09 ../res/alignments/orthologs/10 ../data/tables/nuclear_queries_exons.txt ../res/tables/selected_nuclear_exons.txt 0.5 1.0
