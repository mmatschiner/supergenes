# m_matschiner Thu Feb 25 16:51:09 CET 2021

# Make data directories.
mkdir -p ../data/fastq/gadMor2
mkdir -p ../data/fastq/gadMor_Stat
mkdir -p ../data/fastq/melAeg

# Download fastq files of a genome close to the gadMor2 genome.
cd ../data/fastq/gadMor2
wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
mv XXX* gadMor2_R1.fasta.gz
wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
mv XXX* gadMor2_R2.fasta.gz
cd -

# Download fastq files of the gadMor_Stat genoems.
cd ../data/fastq/gadMor_Stat
wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
mv XXX* gadMor_Stat_R1.fasta.gz
wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
mv XXX* gadMor_Stat_R2.fasta.gz
cd -

# Download fastq files of the melAeg genome.
cd ../data/fastq/melAeg
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/005/ERR2028455/ERR2028455_1.fastq.gz
mv ERR2028455_1.fastq.gz melAeg_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/005/ERR2028455/ERR2028455_2.fastq.gz
mv ERR2028455_2.fastq.gz melAeg_R2.fastq.gz
cd -