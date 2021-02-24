# m_matschiner Tue Jan 31 11:51:18 CET 2017

# Load the ruby module.
module load ruby/2.1.5

# Generate a ld plot for each chromosome.
for i in ../res/plink/*_strict.ld
do
	ruby plot_ld.rb ${i} ${i%.ld}.svg
done
