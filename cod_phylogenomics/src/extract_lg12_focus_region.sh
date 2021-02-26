# m_matschiner Sun Aug 9 23:55:27 CEST 2020

# Set the vcf file.
gzvcf=../res/gatk/gadMor2/original_gadoga/concatenated.masked.vcf.gz

# Load modules.
module load BCFtools/1.9-foss-2018b

# Extract the region in which the bornholm population clusters with cluster 1.
bcftools view -s Gadmor_avc1,Gadmor_avc2,Gadmor_avo1,Gadmor_avo2,Gadmor_bat1,Gadmor_bat2,Gadmor_bor1,Gadmor_bor2,Gadmor_icc1,Gadmor_ico1,Gadmor_ico2,Gadmor_kie1,Gadmor_kie2,Gadmor_lfc3,Gadmor_lfo1,Gadmor_lfo3,Gadmor_low1,Gadmor_low2,Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2 --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/focus_lg12.vcf.gz ${gzvcf} LG12:7000000-8000000
