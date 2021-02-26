# m_matschiner Thu Aug 30 13:51:30 CEST 2018

# Download fastq files of a genome close to the gadMor2 genome.
mkdir -p ../data/fastq/gadmon
if [ ! -f ../data/fastq/gadmon/gadmon_L001_R1.fastq.gz ]
then
    wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
    mv XXX* ../data/fastq/gadmon/gadmon_L001_R1.fastq.gz
    wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
    mv XXX* ../data/fastq/gadmon/gadmon_L002_R1.fastq.gz
fi

# Download fastq files of the gadMor_Stat genoems.
mkdir -p ../data/fastq/gadmoc
if [ ! -f ../data/fastq/gadmoc/gadmoc_L001_R1.fastq.gz ]
then
    wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
    mv XXX* ../data/fastq/gadmon/gadmoc_L001_R1.fastq.gz
    wget https://www.ebi.ac.uk/ena/browser/api/fastq/XXX
    mv XXX* ../data/fastq/gadmon/gadmoc_L002_R1.fastq.gz
fi

# Download the fastq files for haddock (melaeg).
mkdir -p ../data/fastq/melaeg
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/005/ERR2028455/ERR2028455_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/005/ERR2028455/ERR2028455_2.fastq.gz
mv ERR2028455_1.fastq.gz ../data/fastq/melaeg/melaeg_L001_R1.fastq.gz
mv ERR2028455_2.fastq.gz ../data/fastq/melaeg/melaeg_L001_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/006/ERR2028456/ERR2028456_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/006/ERR2028456/ERR2028456_2.fastq.gz
mv ERR2028456_1.fastq.gz ../data/fastq/melaeg/melaeg_L002_R1.fastq.gz
mv ERR2028456_2.fastq.gz ../data/fastq/melaeg/melaeg_L002_R2.fastq.gz

# Download the fastq files for polvir, mermer, gadcha, borsai, and arcgla from ena (from malmstrom et al 2016).
mkdir -p ../data/fastq/polvir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/008/ERR1473878/ERR1473878_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/008/ERR1473878/ERR1473878_2.fastq.gz
mv ERR1473878_1.fastq.gz ../data/fastq/polvir/polvir_L001_R1.fastq.gz
mv ERR1473878_2.fastq.gz ../data/fastq/polvir/polvir_L001_R2.fastq.gz
mkdir -p ../data/fastq/mermer
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/000/ERR1473880/ERR1473880_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/000/ERR1473880/ERR1473880_2.fastq.gz
mv ERR1473880_1.fastq.gz ../data/fastq/mermer/mermer_L001_R1.fastq.gz
mv ERR1473880_2.fastq.gz ../data/fastq/mermer/mermer_L001_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/001/ERR1473881/ERR1473881_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/001/ERR1473881/ERR1473881_2.fastq.gz
mv ERR1473881_1.fastq.gz ../data/fastq/mermer/mermer_L002_R1.fastq.gz
mv ERR1473881_2.fastq.gz ../data/fastq/mermer/mermer_L002_R2.fastq.gz
mkdir -p ../data/fastq/gadcha
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/006/ERR1473886/ERR1473886_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/006/ERR1473886/ERR1473886_2.fastq.gz
mv ERR1473886_1.fastq.gz ../data/fastq/gadcha/gadcha_L001_R1.fastq.gz
mv ERR1473886_2.fastq.gz ../data/fastq/gadcha/gadcha_L002_R2.fastq.gz
mkdir -p ../data/fastq/borsai
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/004/ERR1473884/ERR1473884_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/004/ERR1473884/ERR1473884_2.fastq.gz
mv ERR1473884_1.fastq.gz ../data/fastq/borsai/borsai_L001_R1.fastq.gz
mv ERR1473884_2.fastq.gz ../data/fastq/borsai/borsai_L001_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/005/ERR1473885/ERR1473885_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/005/ERR1473885/ERR1473885_2.fastq.gz
mv ERR1473885_1.fastq.gz ../data/fastq/borsai/borsai_L002_R1.fastq.gz
mv ERR1473885_2.fastq.gz ../data/fastq/borsai/borsai_L002_R2.fastq.gz
mkdir -p ../data/fastq/arcgla
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/002/ERR1473882/ERR1473882_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/002/ERR1473882/ERR1473882_2.fastq.gz
mv ERR1473882_1.fastq.gz ../data/fastq/arcgla/arcgla_L001_R1.fastq.gz
mv ERR1473882_2.fastq.gz ../data/fastq/arcgla/arcgla_L001_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/003/ERR1473883/ERR1473883_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/003/ERR1473883/ERR1473883_2.fastq.gz
mv ERR1473883_1.fastq.gz ../data/fastq/arcgla/arcgla_L002_R1.fastq.gz
mv ERR1473883_2.fastq.gz ../data/fastq/arcgla/arcgla_L002_R2.fastq.gz

# Download the fastq files for the pacific (gadmac) and the greenland cod (gadoga) from ena (from halldorsdottir and arnason 2015).
mkdir -p ../data/fastq/gadmac
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/005/SRR2906345/SRR2906345_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/005/SRR2906345/SRR2906345_2.fastq.gz
mv SRR2906345_1.fastq.gz ../data/fastq/gadmac/gadmac_L001_R1.fastq.gz
mv SRR2906345_2.fastq.gz ../data/fastq/gadmac/gadmac_L001_R2.fastq.gz
mkdir -p ../data/fastq/gadoga
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/003/SRR2906193/SRR2906193_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/003/SRR2906193/SRR2906193_2.fastq.gz
mv SRR2906193_1.fastq.gz ../data/fastq/gadoga/gadoga_L001_R1.fastq.gz
mv SRR2906193_2.fastq.gz ../data/fastq/gadoga/gadoga_L001_R2.fastq.gz
