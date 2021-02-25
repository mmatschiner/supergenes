# m_matschiner Sun Jul 22 14:45:05 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/orthologs/05
rm -f ../res/alignments/orthologs/05/*

# Filter exon alignments by their completeness.
ruby filter_exons_by_missing_data.rb ../res/alignments/orthologs/04 ../res/alignments/orthologs/05 2 150