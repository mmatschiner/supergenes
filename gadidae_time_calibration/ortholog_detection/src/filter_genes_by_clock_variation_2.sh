# m_matschiner Tue Jul 24 10:14:30 CEST 2018

# Go to each alignment directory and start beast analyses from there.
for i in ../res/alignments/orthologs/09/ENS*
do
    cd ${dir}
    trees_file=`basename ${dir}.trees`
    if [ -f ${trees_file} ]
    then
        echo "Skipping directory ${dir} as file ${trees_file} exists."
    else
        sbatch start.slurm
    fi
    cd -
done