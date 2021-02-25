# m_matschiner Wed May 3 15:32:37 CEST 2017

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/fastq

# Use picard-tools' samtofastq function to generate a fastq from each bam file.
for bam in ../res/mapping/bam/*.MT_genome.bam
do
    echo -n "Converting ${bam} to fastq..."
    bam_id=`basename ${bam%.bam}`
    java -jar picard/2.7.1/picard.jar SamToFastq \
        I=${bam} \
        FASTQ="../res/fastq/${bam_id}.R1.fastq" \
        SECOND_END_FASTQ="../res/fastq/${bam_id}.R2.fastq" \
        QUIET=TRUE \
        VALIDATION_STRINGENCY=SILENT
    echo " done."
done