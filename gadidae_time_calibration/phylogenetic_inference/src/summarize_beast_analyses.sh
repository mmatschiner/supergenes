# m_matschiner Tue Aug 7 14:12:30 CEST 2018

# Load the python3 module.
module load python3/3.5.0
module load beast2/2.5.0

# Make the output directory.
mkdir -p ../res/beast/combined

# Make lists of log files.
ls ../res/beast/replicates/r??/91genes.log > ../res/beast/combined/log_files.txt
ls ../res/beast/replicates/r??/91genes_species.trees > ../res/beast/combined/tree_files.txt

# Combine the posterior distributions.
if [ ! -f ../res/beast/combined/91genes.log ]
then
    python3 logcombiner.py -b 25 ../res/beast/combined/log_files.txt > ../res/beast/combined/91genes.log
fi
if [ ! -f ../res/beast/combined/91genes_species.trees ]
then
    python3 logcombiner.py -b 25 ../res/beast/combined/tree_files.txt > ../res/beast/combined/91genes_species.trees
fi
if [ ! -f ../res/beast/combined/91genes_species_100.trees ]
then
    python3 logcombiner.py -b 25 -n 20 ../res/beast/combined/tree_files.txt > ../res/beast/combined/91genes_species_100.trees
fi

# Make a maximum-clade-credibility summary tree.
mcc_tree=../res/beast/combined/91genes_species.tre
if [ ! -f ${mcc_tree} ]
then
    treeannotator -burnin 25 -heights mean ../res/beast/combined/91genes_species.trees ${mcc_tree}
fi