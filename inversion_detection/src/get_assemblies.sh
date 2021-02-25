# m_matschiner Thu Feb 25 17:38:23 CET 2021

# Make data directories.
mkdir -p ../data/assemblies/gadMor_Stat
mkdir -p ../data/assemblies/melAeg

# Download the gadMor_Stat assembly.
cd ../data/assemblies/gadMor_Stat
wget https://www.ebi.ac.uk/ena/browser/api/fasta/XXX
mv XXX* gadMor_Stat.fasta.gz
gunzip gadMor_Stat.fasta.gz
cd -

# Download the melAeg assembly.
cd ../data/assemblies/melAeg
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OLKM01000000.1
mv OLKM01000000* melAeg.fasta.gz
gunzip melAeg.fasta.gz
cd -