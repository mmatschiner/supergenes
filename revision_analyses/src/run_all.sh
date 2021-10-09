# m_matschiner Tue Sep 14 20:53:15 CEST 2021

# Download the gadMor3 genome assembly.
bash get_gadmor3.sh

# Prepare dictionary files for the gadMor3 genome.
bash prepare.sh

# Map reads towards the gadMor3 genome.
bash map.sh

# Merge bam files for individuals with multiple files.
bash merge.sh

# Run variant calling with gatk's haplotypecaller.
bash run_gatk1.sh

# Call genotypes with gatk's genotypegvcf.
bash run_gatk2.sh

# Filter genotypes and apply the mappability mask.
bash apply_filter.sh

# Calculate fst per window for vcf files based on both references.
bash get_window_stats.sh

# Simulate drift to estimate how soon inversions would be lost or fixed.
bash simulate_drift.sh
