# michaelm Mon Aug 31 10:33:44 CEST 2020

# Load modules.
module load BCFtools/1.10.2-GCC-8.3.0
module load Ruby/2.7.1-GCCcore-8.3.0

# Set the compressed input gzvcf.
in_gzvcf=../res/relate/input/gadmor.vcf.gz

# Set the ancestor sequence.
ancestor=../res/relate/input/ancestor.fasta

# Set the mask.
mask=../res/relate/input/mask.fasta

# Make five subsets of vcf and fasta files.
# colinear.
bcftools view -O z -o ${in_gzvcf%.vcf.gz}.colinear.vcf.gz ${in_gzvcf} LG03,LG04,LG05,LG06,LG08,LG09,LG10,LG11,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23
for lg in LG03 LG04 LG05 LG06 LG08 LG09 LG10 LG11 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 LG23
do
    ./fastagrep -t -p ${lg} ${ancestor} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${ancestor%.fasta}.colinear.${lg}.fasta
done
for lg in LG03 LG04 LG05 LG06 LG08 LG09 LG10 LG11 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 LG23
do
    ./fastagrep -t -p ${lg} ${mask} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${mask%.fasta}.colinear.${lg}.fasta
done
# inversion_lg01.
bcftools view -O z -o ${in_gzvcf%.vcf.gz}.inversion_lg01.vcf.gz ${in_gzvcf} LG01:9114741-26192489
./fastagrep -t -p LG01 ${ancestor} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${ancestor%.fasta}.lg01.fasta
./fastagrep -t -p LG01 ${mask} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${mask%.fasta}.lg01.fasta
# inversion_lg02.
bcftools view -O z -o ${in_gzvcf%.vcf.gz}.inversion_lg02.vcf.gz ${in_gzvcf} LG02:18489307-24050282
./fastagrep -t -p LG02 ${ancestor} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${ancestor%.fasta}.lg02.fasta
./fastagrep -t -p LG02 ${mask} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${mask%.fasta}.lg02.fasta
# inversion_lg07.
bcftools view -O z -o ${in_gzvcf%.vcf.gz}.inversion_lg07.vcf.gz ${in_gzvcf} LG07:13606502-23016726
./fastagrep -t -p LG07 ${ancestor} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${ancestor%.fasta}.lg07.fasta
./fastagrep -t -p LG07 ${mask} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${mask%.fasta}.lg07.fasta
# inversion_lg12.
bcftools view -O z -o ${in_gzvcf%.vcf.gz}.inversion_lg12.vcf.gz ${in_gzvcf} LG12:589105-13631347
./fastagrep -t -p LG12 ${ancestor} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${ancestor%.fasta}.lg12.fasta
./fastagrep -t -p LG12 ${mask} | sed 's/>LG0/>/g' | sed 's/>LG/>/g' > ${mask%.fasta}.lg12.fasta

# Adjust snp positions in inversion vcfs and fastas.
# inversion_lg01.
gunzip -c ${in_gzvcf%.vcf.gz}.inversion_lg01.vcf.gz > tmp.vcf
ruby adjust_vcf_positions.rb tmp.vcf -9114740
gzip -c tmp.vcf > ${in_gzvcf%.vcf.gz}.inversion_lg01.adjusted.vcf.gz
rm -f tmp.vcf
ruby trim_fasta.rb ${ancestor%.fasta}.lg01.fasta 9114741 26192489 ${ancestor%.fasta}.lg01.trimmed.fasta
ruby trim_fasta.rb ${mask%.fasta}.lg01.fasta 9114741 26192489 ${mask%.fasta}.lg01.trimmed.fasta
# inversion_lg02.
gunzip -c ${in_gzvcf%.vcf.gz}.inversion_lg02.vcf.gz > tmp.vcf
ruby adjust_vcf_positions.rb tmp.vcf -18489306
gzip -c tmp.vcf > ${in_gzvcf%.vcf.gz}.inversion_lg02.adjusted.vcf.gz
rm -f tmp.vcf
ruby trim_fasta.rb ${ancestor%.fasta}.lg02.fasta 18489307 24050282 ${ancestor%.fasta}.lg02.trimmed.fasta
ruby trim_fasta.rb ${mask%.fasta}.lg02.fasta 18489307 24050282 ${mask%.fasta}.lg02.trimmed.fasta
# inversion_lg07.
gunzip -c ${in_gzvcf%.vcf.gz}.inversion_lg07.vcf.gz > tmp.vcf
ruby adjust_vcf_positions.rb tmp.vcf -13606501
gzip -c tmp.vcf > ${in_gzvcf%.vcf.gz}.inversion_lg07.adjusted.vcf.gz
rm -f tmp.vcf
ruby trim_fasta.rb ${ancestor%.fasta}.lg07.fasta 13606502 23016726 ${ancestor%.fasta}.lg07.trimmed.fasta
ruby trim_fasta.rb ${mask%.fasta}.lg07.fasta 13606502 23016726 ${mask%.fasta}.lg07.trimmed.fasta
# inversoin_lg12.
gunzip -c ${in_gzvcf%.vcf.gz}.inversion_lg12.vcf.gz > tmp.vcf
ruby adjust_vcf_positions.rb tmp.vcf -589104
gzip -c tmp.vcf > ${in_gzvcf%.vcf.gz}.inversion_lg12.adjusted.vcf.gz
rm -f tmp.vcf
ruby trim_fasta.rb ${ancestor%.fasta}.lg12.fasta 589105 13631347 ${ancestor%.fasta}.lg12.trimmed.fasta
ruby trim_fasta.rb ${mask%.fasta}.lg12.fasta 589105 13631347 ${mask%.fasta}.lg12.trimmed.fasta
