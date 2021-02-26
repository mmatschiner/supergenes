# m_matschiner Fri Mar 20 11:32:48 CET 2020

# Load modules.
module load Ruby/2.6.3-GCCcore-8.2.0

# Set the tree directory.
tree_dir=../res/snapp/trees

# Set the color code.
color_code=../data/tables/color_code.txt

# Make the results directory.
mkdir -p ../res/snapp/plots

# Set the dictionary file for the reference genome (with lg sizes).
dict=../data/assemblies/gadMor2.dict

# Make plots for each lg.
for lg_id in LG0{1..9} LG{10..23}
do
    plot=../res/snapp/plots/${lg_id}.svg
    length=`head -n 30 ${dict} | grep ${lg_id} | cut -f 3 | cut -d ":" -f 2`
    ruby plot_window_trees.rb ${tree_dir} ${lg_id} ${length} ${color_code} ${plot}
done
