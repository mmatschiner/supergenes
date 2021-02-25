# m_matschiner Tue Jul 24 10:13:39 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Generate XML input files for beast, for each gene.
for dir in ../res/alignments/orthologs/09/*
do
    ruby beauti.rb -id `basename ${dir}` -n ${dir} -u -s -o ${dir} -c ../data/constraints/constraints.xml -l 5000000
done