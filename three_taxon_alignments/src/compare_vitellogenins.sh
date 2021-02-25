# m_matschiner Thu Feb 11 00:13:38 CET 2021

# Load modules.
module load BEDTools/2.28.0-GCC-8.2.0-2.31.1
module load MAFFT/7.429-GCC-8.2.0-2.31.1-with-extensions
module load Python/3.7.2-GCCcore-8.2.0
module load BLAST+/2.9.0-gompi-2019a

# Make the result directory.
mkdir -p ../res/gene_sequences

# Set the annotation file.
gff=../data/annotation/gadMor2_maker.putative_function.domain_added.aed_0.5_LG12_vtgs.gff

# Extract the vitellogenin exons from the gff.
cat ${gff} | grep "ID=GAMO_00027782" | grep exon > tmp.00027782.gff
cat ${gff} | grep "ID=GAMO_00027783" | grep exon > tmp.00027783.gff
cat ${gff} | grep "ID=GAMO_00027785" | grep exon > tmp.00027785.gff

# Extract lg12 from the gadmor2 assembly.
./fastagrep -t -p LG12 ../data/assemblies/gadMor2/gadMor2.fasta > tmp.gadmor2.lg12.fasta

# Extract the stationary cod version of lg12 from the threeway alignment fasta.
echo ">LG12" > tmp.gadmor_stat.lg12.fasta
cat ../res/alignments/merged/lg12.threaded.refined.finished.fasta | grep -A 1 gadMor_Stat | tail -n 1 >> tmp.gadmor_stat.lg12.fasta

# Copy the stationary cod assembly and make a blast database for it.
./fastagrep -t -p scf7180000170861 ../data/assemblies/gadMor_Stat/gadMor_Stat.fasta > tmp.gadmor_stat.scf.fasta
makeblastdb -in tmp.gadmor_stat.scf.fasta -dbtype nucl

# Get exon sequences from both versions of lg12.
for gene_id in 00027782 00027783 00027785
do
    # Add the gadmor sequence to a new combined fasta.
    echo ">gadMor2" > ../res/gene_sequences/${gene_id}.fasta
    bedtools getfasta -fi tmp.gadmor2.lg12.fasta -bed tmp.${gene_id}.gff | grep -v ">" | tr -d "\n" >> ../res/gene_sequences/${gene_id}.fasta
    echo "" >> ../res/gene_sequences/${gene_id}.fasta
    # Get the gadmor_stat sequence in gadmor2 coordinates.
    echo ">gadMor_Stat_in_gadMor2_coords" >> ../res/gene_sequences/${gene_id}.fasta
    bedtools getfasta -fi tmp.gadmor_stat.lg12.fasta -bed tmp.${gene_id}.gff | grep -v ">" | tr -d "\n" >> ../res/gene_sequences/${gene_id}.fasta
    echo "" >> ../res/gene_sequences/${gene_id}.fasta
    # Get the gadmor_stat sequence from the original gadmor_stat assembly through blast.
    bedtools getfasta -fi tmp.gadmor_stat.lg12.fasta -bed tmp.${gene_id}.gff > tmp.fasta
    python3 find_orthologs.py --overwrite --refine -s 2 -e 0.1 tmp.fasta tmp.gadmor_stat.scf.fasta
    rm -f tmp.fasta
    for fasta in LG12*.fasta
    do
        cat ${fasta} | tail -n 1 | tr -d "\n" >> tmp.fasta
    done
    echo ">gadMor_Stat" >> ../res/gene_sequences/${gene_id}.fasta
    cat tmp.fasta >> ../res/gene_sequences/${gene_id}.fasta
    echo "" >> ../res/gene_sequences/${gene_id}.fasta
    rm tmp.fasta
    rm -f LG12*.fasta
done

# Clean up.
rm -f tmp.*.lg12.fasta
rm -f tmp.*.lg12.fasta.fai
rm -f tmp.000*.gff
rm -f tmp.gadmor_stat.scf.fasta*
