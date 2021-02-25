# m_matschiner Fri May 5 14:25:37 CEST 2017

# Load the mafft module.
module load mafft/7.300

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments

# Concatenate all fasta files into a single file.
echo ">gadmor2" > tmp.fasta
cat ../data/ref/gadMor2.MT_genome.fasta | tail -n +2 >> tmp.fasta
for fasta in ../res/fasta/*.fasta
do
    specimen_id=`basename ${fasta%.MT_genome.fasta}`
    echo ">${specimen_id}"
    cat ${fasta} | tail -n +2
done >> tmp.fasta

# Align all mitochondrial sequences.
mafft --thread 10 --auto tmp.fasta > ../res/alignments/mt_genome.fasta

# Clean up.
rm tmp.fasta