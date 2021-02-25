# m_matschiner Sun Apr 2 00:37:24 CEST 2017

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/nuclear/01

# Make the log directory if it doesn't exist yet.
mkdir -p ../log/misc

# Generate a list of all subject sequence fasta files
# (use readlink and -d so that the absolute paths are specified).
ls -d `readlink -f ../data/subjects/`/danrer.fasta > ../data/subjects/assemblies.txt
ls -d `readlink -f ../res/kollector/nuclear/`/*.fasta >> ../data/subjects/assemblies.txt

# Split the query file into a suitable number of files.
rm -f exons_??.fasta
split -l 300 -d ../data/queries/exons.fasta exons_
for i in exons_??
do
    mv ${i} ../data/queries/${i}.fasta
done

# Search for orthologs to all queries in all subject sequences.
for i in ../data/queries/exons_??.fasta
do
    exon_set_id=`basename ${i%.fasta}`
    out=../log/misc/find_orthologs.${exon_set_id}.out
    sbatch -o ${out} find_orthologs.slurm ${i} ../data/subjects/assemblies.txt ../res/alignments/nuclear/01/
done
