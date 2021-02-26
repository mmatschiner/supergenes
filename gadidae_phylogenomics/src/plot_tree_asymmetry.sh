# m_matschiner Fri Dec 14 09:28:00 CET 2018

# Load modules
module load ruby/2.1.5

# Make the output directory.
mkdir -p ../res/plots

# Set the name order table file.
samples_table=../data/tables/samples.txt
cat ${samples_table} | cut -f 1 > tmp.samples.txt

# Plot pairwise maximum asymmetry values.
for table in ../res/tables/tree_asymmetry*.txt
do	
    # Set the name of the plot file.
    table_id=`basename ${table%.txt}`
    plot=../res/plots/${table_id}.svg

    # Plot the reduced and folded matrix as a heatmap.
    if [ ! -f ${plot} ]
    then
        ruby plot_tree_asymmetry.rb ${table} tmp.samples.txt ${plot}
        echo "Wrote file ${plot}."
    fi
done

# Clean up.
rm -f tmp.samples.txt