# m_matschiner Tue Apr 25 09:59:50 CEST 2017

# Define the external bam directory.
bam_directory=../res/mapping/bam

# Define the chromosome id of the mitochondrial genome.
mitogenome_id="MT_genome"

# Load the samtools module.
module load samtools/1.3.1

# Extract the mitochondrial part from each bam file.
for i in ${bam_directory}/*.bam
do
    file_id=`basename ${i%.bam}`
    # Make sure that the input bam file is not one of those produced by this script.
    if [ ${file_id} == ${file_id%${mitogenome_id}} ]
    then
        mito_bam=../res/mapping/bam/${file_id}.${mitogenome_id}.bam
        if [ ! -f ${mito_bam} ]
        then
            samtools view -b -h ${i} ${mitogenome_id} > ${mito_bam}
            samtools index ${mito_bam}
        fi
    fi
done
