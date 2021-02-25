# m_matschiner Tue Feb 27 20:36:46 CET 2018

# Load the ruby module.
module load ruby/2.1.5

# Use a ruby script to check the alignments.
for i in ../res/alignments/masa_cudalign/pairwise_with*/*/lg??/alignment.00.txt
do
    echo -n "Checking file ${i}..."
    ruby check_alignments.rb ${i} ${i%.txt}_check.txt
    echo " done."
done
