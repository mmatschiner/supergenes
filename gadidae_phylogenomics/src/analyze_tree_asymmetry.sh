# m_matschiner Wed Oct 30 14:51:43 CET 2019

# Make the output directory.
mkdir -p ../res/tables

# Analyze tree asymmetry for different sets of trees.
for dir in ../res/iqtree/clean_trees/{full,no_inversion,lg01_inversion,lg07_inversion}
do
    # Set the output table.
    dir_id=`basename ${dir}`
    table=../res/tables/tree_asymmetry_${dir_id}.txt

    # Analyze tree asymmetry.
    out=../log/misc/tree_asymmetry_${dir_id}.out
    rm -f ${out}
    sbatch -o ${out} analyze_tree_asymmetry.slurm ${dir} ${table}
done