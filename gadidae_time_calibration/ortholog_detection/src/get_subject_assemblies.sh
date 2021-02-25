# m_matschiner Fri Jul 20 15:28:20 CEST 2018

# Make data directories.
mkdir -p ../data/subjects

# Download the genome assembly of Danio rerio.
if [ ! -f ../data/subjects/danrer.fasta ]
then
    wget http://ftp.ensembl.org/pub/release-87/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
    gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
    mv Danio_rerio.GRCz10.dna.toplevel.fa ../data/subjects/danrer.fasta
fi

# Download assemblies from the dryad repository of Malmstrom et al.
if [ ! -f ../data/subjects/arcgla.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120275/fish_1.scf.fasta.gz
    gunzip fish_1.scf.fasta.gz
    mv fish_1.scf.fasta ../data/subjects/arcgla.fasta
fi
if [ ! -f ../data/subjects/borsai.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120523/fish_2.scf.fasta.gz
    gunzip fish_2.scf.fasta.gz
    mv fish_2.scf.fasta ../data/subjects/borsai.fasta
fi
if [  ! -f ../data/subjects/trimin.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120525/fish_3.scf.fasta.gz
    gunzip fish_3.scf.fasta.gz
    mv fish_3.scf.fasta ../data/subjects/trimin.fasta
fi
if [ ! -f ../data/subjects/polvir.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120527/fish_4.scf.fasta.gz
    gunzip fish_4.scf.fasta.gz
    mv fish_4.scf.fasta ../data/subjects/polvir.fasta
fi
if [ ! -f ../data/subjects/melaeg.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120529/fish_5.scf.fasta.gz
    gunzip fish_5.scf.fasta.gz
    mv fish_5.scf.fasta ../data/subjects/melaeg.fasta
fi
if [ ! -f ../data/subjects/mermer.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120531/fish_6.scf.fasta.gz
    gunzip fish_6.scf.fasta.gz
    mv fish_6.scf.fasta ../data/subjects/mermer.fasta
fi
if [ ! -f ../data/subjects/gadcha.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120533/fish_7.scf.fasta.gz
    gunzip fish_7.scf.fasta.gz
    mv fish_7.scf.fasta ../data/subjects/gadcha.fasta
fi
if [ ! -f ../data/subjects/gadarg.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120535/fish_8.scf.fasta.gz
    gunzip fish_8.scf.fasta.gz
    mv fish_8.scf.fasta ../data/subjects/gadarg.fasta
fi
if [ ! -f ../data/subjects/brobro.fasta ]
then
    wget https://datadryad.org/bitstream/handle/10255/dryad.120545/fish_12.scf.fasta.gz
    gunzip fish_12.scf.fasta.gz
    mv fish_12.scf.fasta ../data/subjects/brobro.fasta
fi
if [ ! -f gadmor.fasta.gz ]
then
    wget https://ndownloader.figshare.com/files/5323414
    mv 5323414 gadMor2.fasta.gz
    gunzip gadMor2.fasta.gz
    mv gadMor2.fasta ../data/subjects/gadmor.fasta
fi