# michaelm Sat Feb 15 12:56:30 CET 2020

# Copy the illumina fastq files for the canadian cod to the data directory.
for fq_id in BAT1309Z03 BAT1309Z10 TWI1307Z07 TWI1307Z15 TWI1307Z12 TWI1307Z17 ICO0304Z23 ICO0304Z20 ICC0304Z11 LOF1103Z11 LOF1106Z04 LOF1106Z24 AVE1409Z09 AVE1409Z10 AVE1403Z11 AVE1403Z10 LOW1503Z06 LOW1504Z07 KIE1103Z20 KIE1102Z06 BOR1205Z07 BOR1205Z03
do
    rsync -av michaelm@login.nird.sigma2.no:/projects/NS9003K/projects/aquagenome/data/fastq_files_Z_names/fastq_files/${fq_id}_R1.fastq.gz ../data/reads/Gadmor_${fq_id}_R1.fastq.gz
    rsync -av michaelm@login.nird.sigma2.no:/projects/NS9003K/projects/aquagenome/data/fastq_files_Z_names/fastq_files/${fq_id}_R2.fastq.gz ../data/reads/Gadmor_${fq_id}_R2.fastq.gz
done
