# michaelm Sun Aug 30 13:41:38 CEST 2020

# Load modules
module load BCFtools/1.10.2-iccifort-2019.5.281

# Concatenate phased vcfs per lg.
zcat ../res/beagle/LG01.filtered.vcf.gz > ../res/beagle/concatenated.filtered.vcf
for gzvcf in ../res/beagle/LG0??.filtered.vcf.gz
do
    zcat ${gzvcf} | grep -v "#" >> ../res/beagle/concatenated.filtered.vcf
done

# Compress with bgzip.
bgzip ../res/beagle/concatenated.filtered.vcf
tabix ../res/beagle/concatenated.filtered.vcf.gz
