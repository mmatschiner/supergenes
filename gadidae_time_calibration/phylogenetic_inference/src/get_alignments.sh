# m_matschiner Thu Feb 25 20:07:10 CET 2021

# Make data directories.
mkdir -p ../data/orthologs

# Copy ortholog alignments and the concatenated alignment.
cp -r ../../ortholog_detection/res/alignments/11 ../data/orthologs
cp ../../ortholog_detection/res/alignments/concatenated.phy ../data/orthologs