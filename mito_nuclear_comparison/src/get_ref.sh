# m_matschiner Thu Feb 25 18:13:10 CET 2021

# Make data directories.
mkdir -p ../data/ref/gadMor2

# Download the gadMor2 assembly.
cd ../data/ref/gadMor2/
wget https://ndownloader.figshare.com/files/5323414
mv 5323414 ../data/ref/gadMor2/gadMor2.fasta.gz
gunzip ../data/ref/gadMor2/gadMor2.fasta.gz
cd -
