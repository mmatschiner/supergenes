# michaelm Mon Jan 16 13:34:19 CET 2017

# Load the BLAST+ module.
module load blast+/2.2.29

# Make a BLAST database for the gadMor_Stat assembly.
makeblastdb -in ../data/assemblies/gadMor_Stat.fasta -dbtype nucl

# Make a BLAST database for the melAeg assembly
makeblastdb -in ../data/assemblies/melAeg.fasta -dbtype nucl

# Make a BLAST database for each fo the linkage groups of the gadMor2 assembly.
for i in `seq -w 23`
do
    ./fastagrep -t -p LG${i} ../data/assemblies/gadMor2.fasta > ../data/assemblies/gadMor2_lg${i}.fasta
    makeblastdb -in ../data/assemblies/gadMor2_lg${i}.fasta -dbtype nucl
done