# m_matschiner Thu Jun 15 12:19:16 CEST 2017

# Make the output directories if they don't exist yet.
mkdir -p ../res/fasta/merged_contigs_per_lg_with_ns/gadMor_Stat
mkdir -p ../res/fasta/merged_contigs_per_lg_with_ns/melAeg
mkdir -p ../res/fasta/merged_contigs_per_lg_without_ns/gadMor_Stat
mkdir -p ../res/fasta/merged_contigs_per_lg_without_ns/melAeg

# Load the ruby module.
module load ruby/2.1.5

# Run a ruby script to merge contigs per lg.
for specimen in gadMor_Stat melAeg
do
    for i in ../res/fasta/separate_contigs_per_lg/${specimen}/lg*.fasta
    do
	    file_name=`basename $i`
	    ruby merge_contigs_per_lg.rb ${i} ../res/fasta/merged_contigs_per_lg_with_ns/${specimen}/${file_name} true
	    ruby merge_contigs_per_lg.rb ${i} ../res/fasta/merged_contigs_per_lg_without_ns/${specimen}/${file_name} false
    done
done
