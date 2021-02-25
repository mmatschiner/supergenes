# m_matschiner Wed Jun 14 14:16:39 CEST 2017

# Load the ruby module.
module load ruby/2.1.5

# Use a ruby script to split and clean the list of ordered contigs.
ruby split_and_clean_contig_orders.rb ../data/tables/gadMor_Stat_contig_order.txt
ruby split_and_clean_contig_orders.rb ../data/tables/melAeg_contig_order.txt