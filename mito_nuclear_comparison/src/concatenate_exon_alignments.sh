# m_matschiner Tue Jul 24 13:53:42 CEST 2018

# Load the ruby module.
module load ruby/2.1.5
module load python3/3.5.0

# Convert all alignments to nexus format and remove the outgroup sequence.
for align in ../res/alignments/nuclear/03/*.fasta
do 
    if [ ! -f ${align%.fasta}.nex ]
    then
        python3 convert.py -x danrer -f nexus ${align} ${align%.fasta}.nex
    fi
done

# Use a ruby script to concatenate all remaining alignments.
find ../res/alignments/nuclear/03 -size 0 -delete
ruby concatenate.rb -i ../res/alignments/nuclear/03/*.nex -o ../res/alignments/nuclear/concatenated.phy -f phylip