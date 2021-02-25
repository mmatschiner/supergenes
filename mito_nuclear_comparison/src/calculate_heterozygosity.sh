# m_matschiner Tue Apr 25 17:29:45 CEST 2017

# Load the vcftools module.
module load vcftools/0.1.14.zlib.1.2.8

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/tables

# Output a header.
echo -e "sample_id\tn_hom\tn_het" > ../res/tables/heterozygosity.txt

# Calculate the number of homo- and heterozygous snps for all samples.
for i in ../res/vcf/*.vcf
do
    sample_id=`basename ${i%.vcf}`
    n_hom=`vcftools --vcf ${i} --minQ 200 --extract-FORMAT-info GT --stdout | grep "1/1" | wc -l`
    n_het=`vcftools --vcf ${i} --minQ 200 --extract-FORMAT-info GT --stdout | grep "0/1" | wc -l`
    echo -e "${sample_id}\t${n_hom}\t${n_het}" >> ../res/tables/heterozygosity.txt
done

# Clean up.
rm out.log