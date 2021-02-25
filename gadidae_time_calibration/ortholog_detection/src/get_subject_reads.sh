# michaelm Wed Mar 8 15:55:37 CET 2017

# Make the data directory if it doesn't exist yet.
mkdir -p ../data/reads

# Download the fastq files for the Pacific and the Greenland cod from ENA (from Halldorsdottir and Arnason 2015).
if [ ! -f ../data/reads/gadmac_R1.fastq.gz ]
then
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/005/SRR2906345/SRR2906345_1.fastq.gz
    mv SRR2906345_1.fastq.gz ../data/reads/gadmac_R1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/005/SRR2906345/SRR2906345_2.fastq.gz
    mv SRR2906345_2.fastq.gz ../data/reads/gadmac_R2.fastq.gz
fi
if [ ! -f ../data/reads/gadoga_R1.fastq.gz ]
then
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/003/SRR2906193/SRR2906193_1.fastq.gz
    mv SRR2906193_1.fastq.gz ../data/reads/gadoga_R1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/003/SRR2906193/SRR2906193_2.fastq.gz
    mv SRR2906193_2.fastq.gz ../data/reads/gadoga_R2.fastq.gz
fi