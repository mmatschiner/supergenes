# michaelm Wed Jan 18 15:27:42 CET 2017

# Load the Ruby module.
module load ruby/2.1.5

# Mark missing and repetitive regions in the target region fasta files.
for i in ../data/assemblies/gadMor2_lg??_*.fasta
do
    echo -n "Analyzing file ${i}..."
    ruby characterize_features.rb ${i} ${i%.fasta}_features.txt 100 0.5 0.25
    echo " done."
done