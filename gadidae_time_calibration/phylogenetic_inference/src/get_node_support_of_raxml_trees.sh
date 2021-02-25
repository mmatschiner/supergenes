# m_matschiner Tue Jul 24 14:42:35 CEST 2018

# Make the results directory if it doesn't exist yet.
mkdir -p ../res/tables

# Load the python3 module.
module load python3

# Get the mean node support.
for tree in ../res/raxml/ENSDARG*.tre
do
    echo -ne "${tree}\t"
    python3 get_mean_node_support.py ${tree}
done | sort -k 2 > ../res/tables/node_support.txt