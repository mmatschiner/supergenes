# m_matschiner Sun Jul 22 15:11:46 CEST 2018

# Make the result directory.
mkdir -p ../res/alignments/orthologs/06
rm -f ../res/alignments/orthologs/06/*

# Remove the reference sequence from all alignments.
for i in ../res/alignments/orthologs/05/*.fasta
do
    tail -n +3 ${i} > ../res/alignments/orthologs/06/`basename ${i}`
done