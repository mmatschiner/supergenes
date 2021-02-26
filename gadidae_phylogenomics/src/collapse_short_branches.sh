# m_matschiner Fri Dec 14 11:09:46 CET 2018

# Make the log directory.
mkdir -p ../log/misc

# Collapse nodes connected by short branches. 
tree_dir=../res/iqtree/trees
clean_tree_dir=../res/iqtree/clean_trees/full
mkdir -p ${clean_tree_dir}
for last_char in {0..9}
do
    out=../log/misc/collapse_nodes_${last_char}.out
    rm -f ${out}
    sbatch -o ${out} collapse_short_branches.slurm ${tree_dir} ${last_char} ${clean_tree_dir}
done

