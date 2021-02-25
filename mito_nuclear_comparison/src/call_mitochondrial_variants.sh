# m_matschiner Tue Apr 25 17:00:24 CEST 2017

# Load  the samtools and bcftools modules.
module load samtools/1.3.1
module load bcftools/1.3

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/vcf

# Remove results of previous interupted runs of this script.
rm -rf tmp.vcf

# Set the reference.
ref=../data/ref/gadMor2.fasta

# Use mpileup and bcftools to call variants in all sam files.
for bam in ../res/mapping/bam/*.MT_genome.bam
do
    sample_id=`basename ${bam%.MT_genome.bam}`
    if [ ! -f ../res/vcf/${sample_id}.vcf ]
    then
        samtools mpileup -uvf ${ref} ${bam} | bcftools call -mv > tmp.vcf
        mv tmp.vcf ../res/vcf/${sample_id}.vcf
    fi
done