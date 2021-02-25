# m_matschiner Tue Jul 3 00:43:39 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Get the number of complete alignment sites.
for i in ../res/alignments/merged/*finished.fasta
do
    echo -ne "${i}\t"
    ruby get_number_of_complete_alignment_sites.rb ${i}
done
