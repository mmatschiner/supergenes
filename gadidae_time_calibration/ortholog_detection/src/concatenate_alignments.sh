# m_matschiner Tue Jul 24 13:53:42 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Use a ruby script to concatenate all alignments.
ruby concatenate.rb -i ../res/alignments/orthologs/11/*.nex -o ../res/alignments/orthologs/concatenated.phy -f phylip