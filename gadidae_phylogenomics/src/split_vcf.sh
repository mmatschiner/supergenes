# m_matschiner Wed Oct 30 11:01:59 CET 2019

# Load modules.
module load vcftools/0.1.14.zlib.1.2.8
module load bcftools/1.6

# Set the vcf file.
vcfgz=../res/gatk/concatenated.masked.variable.vcf.gz

# Set the inversion table.
inversion_table=../data/tables/inversion_limits.txt

# Make temporary bed files.
echo "chrom\tchromStart\tchromEnd" > tmp.no_inversion.txt
echo "chrom\tchromStart\tchromEnd" > tmp.lg01_inversion.txt
echo "chrom\tchromStart\tchromEnd" > tmp.lg02_inversion.txt
echo "chrom\tchromStart\tchromEnd" > tmp.lg07_inversion.txt
cat ${inversion_table} | grep -v "#" >> tmp.no_inversion.txt
cat ${inversion_table} | grep LG01 >> tmp.lg01_inversion.txt
cat ${inversion_table} | grep LG02 >> tmp.lg02_inversion.txt
cat ${inversion_table} | grep LG07 >> tmp.lg07_inversion.txt

# Make a vcf excluding all inversions.
for mask in tmp.no_inversion.txt tmp.lg??_inversion.txt
do
    mask_id=`echo ${mask} | cut -d "." -f 2`
    split_vcfgz=../res/gatk/concatenated.masked.variable.${mask_id}.vcf.gz
    if [ ! -f ${split_vcfgz} ]
    then
        if [ ${mask} == tmp.no_inversion.txt ]
        then
            vcftools --gzvcf ${vcfgz} --exclude-bed ${mask} --recode --recode-INFO-all --stdout | \
            bcftools view -O z -o ${split_vcfgz}
        else
            vcftools --gzvcf ${vcfgz} --bed ${mask} --recode --recode-INFO-all --stdout | \
            bcftools view -O z -o ${split_vcfgz}
        fi
    fi
done