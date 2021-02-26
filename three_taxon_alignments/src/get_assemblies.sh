# m_matschiner Thu Feb 25 16:36:00 CET 2021

# Make data directories.
mkdir -p ../data/assemblies/gadMor2
mkdir -p ../data/assemblies/gadMor_Stat
mkdir -p ../data/assemblies/melAeg

# Download the gadMor2 assembly.
cd ../data/assemblies/gadMor2/
wget https://ndownloader.figshare.com/files/5323414
mv 5323414 gadMor2.fasta.gz
gunzip gadMor2.fasta.gz
cd -

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