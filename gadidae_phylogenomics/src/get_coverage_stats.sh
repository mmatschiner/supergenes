# m_matschiner Tue Sep 25 09:38:04 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Get coverage statistics for all bam files.
for txt in ../res/mapping/*.covdist.txt
do
    ruby get_coverage_stats.rb ${txt} > ${txt%.covdist.txt}.covstats.txt
done