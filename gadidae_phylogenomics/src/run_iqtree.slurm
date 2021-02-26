#!/bin/bash

# m_matschiner Sat Dec 15 19:52:45 CET 2018

#SBATCH --job-name=iqtree
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=5

# Set up.
source /cluster/bin/jobsetup
module load iqtree/1.6.8
module load R/3.4.4

# Get the command-line arguments.
alignment_dir=`readlink -f ${1}`
second_last_char=${2}
last_char=${3}
tree_dir=`readlink -f ${4}`
outgroup=${5}
r_script=`readlink -f reroot_tree.r`

# Change to the temporary directory.
cd ${SCRATCH}

# Run iqtree for all alignments of the permissive dataset.
for align in ${alignment_dir}/*${second_last_char}${last_char}.phy
do
    marker_id=`basename ${align%.phy}`
    echo -n "Analyzing alignment ${marker_id}.phy..."
    info=${tree_dir}/${marker_id}.info.txt
    tree=${tree_dir}/${marker_id}.tre
    if [ ! -f ${tree} ]
    then
	cp -f ${align} tmp.iqtree.phy
	
	# Ensure that data for the outgroup is included.
	outgroup_included=`cat tmp.iqtree.phy | grep ${outgroup} | wc -l`
	if [[ ${outgroup_included} == 0 ]]
	then
	    echo "INFO: Skipping tree ${tree} as it doesn't include data for ${outgroup}."
	else
	
	    # Run iqtree.
	    iqtree -nt 5 -s tmp.iqtree.phy -m TEST --runs 2 -quiet

	    # Reroot the tree.
	    Rscript ${r_script} tmp.iqtree.phy.treefile tmp.rooted.tre ${outgroup}

	    # Return the result files.
	    mv -f tmp.rooted.tre ${tree}
	    mv -f tmp.iqtree.phy.iqtree ${info}

	    # Clean up.
	    rm -f tmp.iqtree.*
	fi
    fi
    echo " done."
done