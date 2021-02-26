# m_matschiner Tue Oct 16 00:59:16 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory.
mkdir -p ../res/aim/xml

# Write the aim xml file.
ruby aim_prep.rb -d ../res/aim/alignments -o polvir -a "lognormal(8.56,0.08)" -l 1000000000 -x ../res/aim/xml/aim_age_prior.xml