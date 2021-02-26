# m_matschiner Wed Oct 30 10:24:23 CET 2019

# Load modules.
module load ruby/2.1.5

# Set the samples file.
samples=../data/tables/samples.txt

# Extract the first column from the samples file.
cat ${samples} | cut -f 1 > tmp.txt

# Plot the results.
for dsuite_res in ../res/dsuite/*/samples__BBAA.txt ../res/dsuite/*/samples__tree.txt
do
    ruby plot_d.rb ${dsuite_res} tmp.txt 0.3 ${dsuite_res%.txt}.svg
done

# Clean up.
rm -f tmp.txt