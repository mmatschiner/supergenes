# m_matschiner Tue May 1 16:06:57 CEST 2018

# Make the queries directory if it doesn't exist yet.
mkdir -p ../data/queries

# Remove the query file if it should already exist.
rm -f ../data/queries/gadmor_exons.fasta
rm -f ../data/queries/gadmor_exons_nucl.fasta

# Copy the alignment directory from the pipefish analysis pipeline.
cp -r ~/priv1TB/pipefishes/ortholog_identification/res/alignments/nuclear/01 ../data/queries

# Generate two files with only gadmor sequences in nucleotide and amino-acid format.
for i in ../data/queries/01/ENS*_nucl.fasta
do
    exon_id=`basename ${i%_nucl.fasta}`
    echo ">${exon_id}" >> ../data/queries/gadmor_exons_nucl.fasta
    cat ${i} | grep -A 1 gadmor | tail -n 1 | tr -d "-" >> ../data/queries/gadmor_exons_nucl.fasta
done
for i in ../data/queries/01/ENS*[0-9].fasta
do
    exon_id=`basename ${i%.fasta}`
    echo ">${exon_id}" >> ../data/queries/gadmor_exons.fasta
    cat ${i} | grep -A 1 gadmor | tail -n 1 | tr -d "-" >> ../data/queries/gadmor_exons.fasta
done

# Clean up.
rm -rf ../data/queries/01
