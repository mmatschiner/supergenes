# m_matschiner Sun Aug 30 14:26:51 CEST 2020

# Load modules.
module load BCFtools/1.10.2-iccifort-2019.5.281
module load Ruby/2.7.1-GCCcore-8.3.0

# Set the compressed vcf.
in_gzvcf=../res/beagle/concatenated.filtered.vcf.gz

# Set the ancestor id.
anc_id="Gadcha"

# Set the compressed output vcf.
out_gzvcf=../res/beagle/concatenated.filtered.anc_hets_removed.vcf.gz

# Set the output mask file.
bed=../res/beagle/anc_hets.bed

# Uncompress the vcf.
bcftools view -o tmp.vcf ${in_gzvcf}

# Use a ruby script to remove sites that are heterozygous in the ancestor.
ruby remove_sites_heterozygous_in_anc.rb tmp.vcf ${anc_id} ${out_gzvcf%.gz} ${bed}

# Compress the output vcf.
bgzip ${out_gzvcf%.gz}
tabix ${out_gzvcf}

# Clean up.
rm -f tmp.vcf
