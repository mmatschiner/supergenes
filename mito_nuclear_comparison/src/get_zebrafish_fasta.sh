# m_matschiner Thu Feb 25 18:29:04 CET 2021

# Make data directories.
mkdir -p ../data/subjects

# Download the genome assembly of Danio rerio.
if [ ! -f ../data/subjects/danrer.fasta ]
then
    wget http://ftp.ensembl.org/pub/release-87/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
    gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
    mv Danio_rerio.GRCz10.dna.toplevel.fa ../data/subjects/danrer.fasta
fi
