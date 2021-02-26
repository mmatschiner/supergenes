# m_matschiner Wed Oct 30 10:07:50 CET 2019

# Set the outgroup.
outgroup=polvir

# Make the log directory.
mkdir -p ../log/iqtree

# Set the alignment directory.
alignment_dir=../res/iqtree/alignments

# Set and make the tree directory.
tree_dir=../res/iqtree/trees
mkdir -p ${tree_dir}

# Run iqtree.
for last_char in {0..9}
do
    for second_last_char in {0..9}
    do
        out=../log/iqtree/iqtree_${second_last_char}${last_char}.out
        rm -f ${out}
        sbatch -o ${out} run_iqtree.slurm ${alignment_dir} ${second_last_char} ${last_char} ${tree_dir} ${outgroup}
    done
done

