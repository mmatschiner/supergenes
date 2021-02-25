# m_matschiner Wed Jul 4 00:13:44 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Set the window size.
window_size=100000

# Run a ruby script to calculate the number of missing sites of each sequence and the distances between sequences.
for i in ../res/alignments/merged/lg??.threaded.refined.finished.fasta
do
    alignment=${i}
    table=${alignment%.fasta}.txt
    ruby get_distances_in_sliding_windows.rb ${alignment} ${window_size} > ${table}
    echo "Wrote file ${table}."
done
