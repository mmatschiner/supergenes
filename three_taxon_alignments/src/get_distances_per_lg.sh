# michaelm Wed Aug 12 19:45:16 CEST 2020

# Load the ruby module.
module load Ruby/2.7.1-GCCcore-8.3.0

# Run a ruby script to calculate the number of missing sites of each sequence and the distances between sequences.
for i in ../res/alignments/merged/lg{01,02,07,12}.threaded.refined.finished.fasta
do
    alignment=${i}
    table=${i%.fasta}.dists.txt
    ruby get_distances.rb ${alignment} > ${table}
done
