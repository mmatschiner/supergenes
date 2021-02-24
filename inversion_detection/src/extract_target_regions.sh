# michaelm Mon Jan 16 16:40:46 CET 2017

# Load the Ruby module.
module load ruby/2.1.5

# Extract target regions from the gadMor2 assembly.
ruby extract_target_regions.rb ../data/assemblies/gadMor2.fasta ../data/tables/all_regions.txt