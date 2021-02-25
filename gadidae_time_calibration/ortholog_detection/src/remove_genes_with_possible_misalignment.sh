# m_matschiner Tue Jul 24 13:49:46 CEST 2018

# Make the output directory.
mkdir ../res/alignments/orthologs/11

# Copy all alignments to the new directory but remove two manually selected alignments.
for alignment in ../res/alignments/orthologs/10/*.nex
do
    cp ${alignment} ../res/alignments/orthologs/11
done
rm ../res/alignments/orthologs/11/ENSDARG00000101220.nex
rm ../res/alignments/orthologs/11/ENSDARG00000057997.nex