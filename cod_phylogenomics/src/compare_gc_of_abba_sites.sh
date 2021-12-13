# m_matschiner Fri Aug 28 17:14:56 CEST 2020

# Load modules.
module load BCFtools/1.10.2-iccifort-2019.5.281
module load Ruby/2.7.1-GCCcore-8.3.0

# Set the compressed vcf file.
gzvcf=../res/gatk/gadMor2/original_gadoga/inversion_lg01.full.vcf.gz

# Extract a subset from the vcf for individuals from the trio supporting gene conversion or double crossover in twillinggate.
if [ ! -f tmp.vcf ]
then
    bcftools view -s 'Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_ico1,Gadmor_ico2,Gadmor_avo1,Gadmor_avo2' -m 2 -M 2 -a --min-ac=1 ${gzvcf} > tmp.vcf
fi

# Use a ruby script to compare gc contents of bbaa, abba, and baba sites.
ruby compare_gc_of_abba_sites.rb tmp.vcf 'Gadmor_lfo1,Gadmor_lfo2,Gadmor_ico1,Gadmor_ico2,Gadmor_avo1,Gadmor_avo2' 'Gadmor_two1,Gadmor_two2' 'Gadmor_twc1,Gadmor_twc2' ../res/tables/gc_bbaa.txt ../res/tables/gc_abba.txt ../res/tables/gc_baba.txt

# Clean up.
rm -f tmp.vcf
