# m_matschiner Thu Jun 18 11:48:39 CEST 2020

# Make the bam directory.
mkdir -p ../data/bam

# Make the log directory.
mkdir -p ../log/misc

# Copy mapped reads from the cod_phylogenomics directory, excluding outgroup species.
for bam in ../../cod_phylogenomics/res/mapping/Gadmor_????.bam ../../cod_phylogenomics/res/mapping/Gadcha.bam
do
    bam_id=`basename ${bam%.bam}`
    if [ ! -f ../data/bam/${bam_id}.bam ]
    then
        cp ${bam} ../data/bam
        cp ${bam}.bai ../data/bam
    fi
done
