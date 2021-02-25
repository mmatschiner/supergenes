# m_matschiner Fri Jul 20 15:16:20 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/nuclear/02

# Filter sequences by their blast bitscores.
ruby filter_sequences_by_bitscore.rb ../res/alignments/nuclear/01/ ../res/alignments/nuclear/02/ 0.9