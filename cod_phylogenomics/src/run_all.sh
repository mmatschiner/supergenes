# m_matschiner Fri Nov 1 20:43:48 CET 2019

## divergence-time estimation with snapp.
# Get published illumina fastq files for selected specimens.
bash get_published_reads.sh

# Get illumina fastq files for the specimen used for the coastal assembly.
bash get_gadMor_Stat_reads.sh

# Get fastq files for newly sequenced specimens.
bash get_new_reads.sh

# Download the gadMor2 assembly and prepare an alternative reference for melAeg to use for tests of reference bias in d-statistics.
bash get_assemblies.sh

# Map reads of all taxa against the reference to generate bam files.
bash map.sh

# Merge bam files from several libraries.
bash merge.sh

# Run variant calling with gatk's haplotypecaller.
bash run_gatk1.sh

# Call genotypes with gatk's genotypegvcf.
bash run_gatk2.sh

# Filter genotypes and apply the mappability mask.
bash apply_filter.sh

# Concatenate the vcf files of all lgs.
bash concatenate_vcfs.sh

# Split the vcf file into separate vcfs per inversion and colinear regions.
bash split_vcf.sh

# Make xml files for snapp analyses of regions (colinear, inversions).
bash make_snapp_xmls_regions.sh

# Prepare directories for snapp analyses of regions.
bash prepare_snapp_regions.sh

# Run snapp for regions.
bash run_snapp_regions.sh

# Combine snapp results per region.
bash combine_snapp_results_regions.sh

# Make xml files for snapp analyses of sliding windows.
bash make_snapp_xmls_windows.sh

# Prepare directories for snapp analyses of sliding windows.
bash prepare_snapp_analyses_windows.sh

# Run snapp for sliding windows (resume repeatedly).
bash run_snapp_windows.sh

# Combine snapp results from two replicates per sliding window.
bash combine_snapp_results_windows.sh

# Select sliding window snapp results for plotting, excluding those with poor convergence.
bash select_combined_snapp_windows.sh

# Produce plots of the sliding-window snapp phylogenies.
bash plot_window_trees.sh

## tests for introgression with dsuite.
# Run dsuite separately for each region, and for each vcf, to test for introgression and possible gene conversion in the inversions.
bash run_dsuite.sh

# Plot all dsuite results.
bash plot_dsuite.sh

## Additional tests.
# Calculate population parameters within and outside of inversion regions.
bash get_population_parameters.sh

# Compare the gc-content of abba sites with that of other sites as a test for gene conversion.
bash compare_gc_of_abba_sites.sh

# Test if reference bias in the gadus ogac sample could be responsible for introgression signals.
bash quantify_reference_bias.sh

# Extract the region on lg12 in which bornholm clusters with cluster1.
bash extract_lg12_focus_region.sh

# Do an ancestry painting for the region on lg12.
bash paint_ancestry_for_lg12_focus_region.sh
