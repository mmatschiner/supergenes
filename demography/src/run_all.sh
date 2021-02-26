# m_matschiner Thu Jun 11 15:10:57 CEST 2020

# Copy mapped reads from cod_phylogenomics directory.
bash get_mapped_reads.sh

# Run variant calling with gatk's haplotypecaller.
bash run_gatk1.sh

# Call genotypes with gatk's genotypegvcf.
bash run_gatk2.sh

# Filter genotypes.
bash apply_filter.sh

# Phase the split vcf files.
bash run_beagle.sh

# Concatenate the phased per-lg vcfs.
bash concatenate_phased_vcfs.sh

# Remove sites from the vcf that are heterozygous in the individual representing the ancestral sequence.
bash remove_sites_heterozygous_in_anc.sh

# Prepare input for relate.
bash prepare_relate_input.sh

# Split the relate input by region.
bash split_relate_input_by_region.sh

# Calculate a mutation rate for relate analyses.
bash get_mutation_rate.sh

# Download the program relate and place it in in the bin directory.
echo "The program relate must be downloaded manually from https://myersgroup.github.io/relate. The downloaded directory must then be placed in ../bin and called \"relate\"."
exit

# Run relate for colinear regions.
bash run_relate_colinear1.sh

# Run relate again for colinear regions to infer the demography.
bash run_relate_colinear2.sh

# Run relate for inversion regions regions.
bash run_relate_inversion_regions1.sh

# Run relate again for inversion regions to infer the demography.
bash run_relate_inversion_regions2.sh

# Use msprime to simulate under a scenario of population divergence followed by an extreme bottleneck in one pop, checking whether that bottleneck can still be detected today.
bash simulate.sh

# Summarize the simulations.
bash summarize_simulations.sh

# Run relate with the simulated data.
bash run_relate_simulated_wrapper.sh

# Plot relate results based on empirical and simulated data.
bash plot_relate_results.sh

# Check the repeat content (this is done here because it uses the mask that is also used for relate).
bash check_repetitive_content.sh
