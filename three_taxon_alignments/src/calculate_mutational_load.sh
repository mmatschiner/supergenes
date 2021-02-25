# m_matschiner Mon Dec 7 16:42:44 CET 2020

# Set the annotation file.
gff=../data/annotation/gadMor2_annotation_filtered_only_gene_models_extract.gff
cd ../data/annotation
gunzip gadMor2_annotation_filtered_only_gene_models_extract.gff.gz
cd -

# Make the results directory.
mkdir -p ../res/tables

# Use a ruby script to calculate the mutational load in each chromosome.
for fasta in ../res/alignments/lg??*.threaded.refined.finished.fasta
do
	lg_id=`basename ${fasta} | cut -d "." -f 1`
	out=../res/tables/mutation_load.${lg_id}.txt
	window_size=1000000
	ruby calculate_mutational_load.rb ${fasta} ${lg_id} ${gff} ${window_size} ${out}
done