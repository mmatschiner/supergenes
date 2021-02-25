# m_matschiner Thu May 4 00:50:07 CEST 2017

# Make assemblies separately for each mitochondrial bam file.
for bam in ../res/mapping/bam/*.MT_genome.bam
do
    bam_id=`basename ${bam%.MT_genome.bam}`
    if [ ! -f ../res/fasta/${bam_id}.MT_genome.fasta ]
    then
        bash make_mitochondrial_assembly.sh ${bam_id}
    fi
done