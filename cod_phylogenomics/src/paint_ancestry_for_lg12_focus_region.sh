# m_matschiner Mon Aug 10 00:04:07 CEST 2020

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0

# Make the output directory.
mkdir -p ../res/ancestry_painting

# Download the ancestry painting scripts.
if [ ! -f get_fixed_site_gts.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rb
fi
if [ ! -f plot_fixed_site_gts.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/plot_fixed_site_gts.rb
fi

# Uncompress the vcf file.
gunzip -c ../res/gatk/focus_lg12.vcf.gz > tmp.vcf

# Get genotypes at fixed sites.
ruby get_fixed_site_gts.rb tmp.vcf ../res/ancestry_painting/fixed_sites.txt Gadmor_low1,Gadmor_low2,Gadmor_kie1,Gadmor_kie2,Gadmor_lfc3 Gadmor_bor1,Gadmor_bor2 Gadmor_avc1,Gadmor_avc2,Gadmor_avo1,Gadmor_avo2,Gadmor_bat1,Gadmor_bat2,Gadmor_icc1,Gadmor_ico1,Gadmor_ico2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2 0.9 0.9

# Plot the genotypes at fixed sites.
ruby plot_fixed_site_gts.rb ../res/ancestry_painting/fixed_sites.txt ../res/ancestry_painting/fixed_sites.svg 0.9 1

# Clean up,
rm -f tmp.vcf
