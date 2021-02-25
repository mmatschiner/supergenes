# m_matschiner Tue Jul 24 13:59:15 CEST 2018

# Load the raxml and python3 modules.
module load raxml/8.2.4
module load python3/3.5.0

# Make the result directory.
mkdir -p ../res/raxml

# Run raxml for each gene alignment.
for alignment in ../data/orthologs/11/*.nex
do
    # Clean up potential leftover results from earlier runs.
    rm -f RAxML_*
    rm -f tmp.phy

    # Convert the alignment to phylip format.
    python3 convert.py -f phylip ${alignment} tmp.phy

    # Run raxml.
    raxmlHPC -s tmp.phy -n tmp -m GTRCAT -f a -p ${RANDOM} -x ${RANDOM} -N 100

    # Move the tree and info files to the result directory.
    alignment_id=`basename ${alignment%.nex}`
    mv RAxML_bipartitions.tmp ../res/raxml/${alignment_id}.tre
    mv RAxML_info.tmp ../res/raxml/${alignment_id}.info

    # Clean up.
    rm -f RAxML_*
    rm -f tmp.phy
    rm -f tmp.phy.reduced
done