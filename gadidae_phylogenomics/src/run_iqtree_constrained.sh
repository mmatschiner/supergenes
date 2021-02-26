# m_matschiner Wed Oct 30 22:21:59 CET 2019

# Make iqtree trees for all hypotheses.
for hypothesis in h01 h02
do
    # Set the prefix of constraint trees.
    constraint_prefix=../data/constraints/${hypothesis}

    # Make the log directory.
    mkdir -p ../log/iqtree_constrained/${hypothesis}

    # Set the alignment directory.
    alignment_dir=../res/iqtree/alignments

    # Set the tree directory.
    res_dir=../res/iqtree_constrained/${hypothesis}

    # Make the tree directory.
    mkdir -p ${res_dir}

    # Run iqtree.
    for lg in LG{01..23}
    do
        out=../log/iqtree_constrained/${hypothesis}/iqtree_lg${lg}.out
        rm -f ${out}
        sbatch -o ${out} run_iqtree_constrained.slurm ${alignment_dir} ${lg} ${constraint_prefix} ${res_dir}
    done
done