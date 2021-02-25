# m_matschiner Tue Aug 7 12:26:22 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the result directories.
rm -rf ../res/alignments/12
mkdir -p ../res/alignments/12
mkdir -p ../res/tables

# Copy alignments of genes with high node support to a new directory.
while read -r line
do
    alignment_file=`echo ${line} | cut -d " " -f 1`
    gene_id=`basename ${alignment_file%.tre}`
    node_support=`echo ${line} | cut -d " " -f 2`
    if (( $(bc <<< "${node_support}>=65") > 0 ))
    then
        cp ../data/orthologs/11/${gene_id}.nex ../res/alignments/12
    fi
done < ../res/tables/node_support.txt

# Determine the chromosomal location for each gene. 
ruby remove_alignments_in_inversions.rb ../res/alignments/12 ../data/tables/selected_nuclear_exons.txt ../data/tables/gadmor_ortholog_regions.txt ../data/tables/inversion_limits.txt ../res/tables/gadmor_gene_locations.txt

# Remove genes that fall within inverted regions.
while read -r line
do
    gene_id=`echo ${line} | cut -d " " -f 1`
    remove_gene=`echo ${line} | cut -d " " -f 4`
    if [ ${remove_gene} == "true" ]
    then
        rm -f ../res/alignments/12/${gene_id}.nex
    fi
done < ../res/tables/gadmor_gene_locations.txt