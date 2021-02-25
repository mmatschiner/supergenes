# m_matschiner Sun Jul 22 15:34:34 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory.
mkdir -p ../res/alignments/orthologs/08
rm -rf ../res/alignments/orthologs/08/*

# Remove genes with too few exon alignments.
ruby filter_genes_by_exon_number.rb ../res/alignments/orthologs/07 ../res/alignments/orthologs/08 ../res/tables/gadmor_ortholog_regions.txt ../data/tables/nuclear_queries_exons.txt 3