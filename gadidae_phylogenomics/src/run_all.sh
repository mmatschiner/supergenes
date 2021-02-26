# m_matschiner Thu Aug 30 13:14:29 CEST 2018

## Analyses of divergence times and introgression with AIM.
# Download the gadMor2 assembly.
bash get_assemblies.sh

# Generate dictionary and index files for the reference.
bash prepare.sh

# Download fastq files.
bash get_fastqs.sh

# Map reads of all taxa against the reference to generate bam files.
bash map.sh

# Merge bam files from several libraries.
bash merge.sh

# Calculate the coverage distribution for each bam file.
bash get_coverage_distribution.sh

# Get coverage statistics for each bam file.
bash get_coverage_stats.sh

# Run variant calling with gatk's haplotypecaller.
bash run_gatk1.sh

# Combine g.vcf files produced by haplotypecaller.
bash run_gatk2.sh

# Apply masks to the vcfs.
bash apply_mask.sh

# Extract alignments in phylip format from the masked vcf.
bash extract_alignments_from_vcfs.sh

# Remove info files no longer needed.
bash move_info_files_to_nobackup.sh

# Get stats for all alignments, these are written to the info files.
bash get_alignment_stats.sh

# Make a table summarizing the alignment statistics from the info files.
bash make_alignment_stats_table.sh

# Download software.
bash get_software.sh

# Select alignments for aim analyses based on position (excluding inversion regions), number of variant sites and number of hemiplasies.
bash select_alignments_for_aim.sh

# Make an xml for aim.
bash make_aim_xml.sh

# Prepare the directories for aim analyses.
bash prepare_aim_analyses.sh

# Run replicate aim analyses.
bash run_aim.sh

# Analyze aim results.
bash analyze_aim_results.sh


## Introgression tests with Dsuite.
# Remove monomorphic sites from the per-lg vcf files.
bash remove_monomorphic_sites.sh

# Concatenate the vcf files of all lgs.
bash concatenate_vcfs.sh

# Generate four further vcf files, one for each inversion region and one excluding all inversion regions.
bash split_vcf.sh

# Run dsuite for d-statistics.
bash run_dsuite_dtrios.sh

# Plot the results of the dsuite dtrios analysis as a heatmap.
bash plot_dsuite_dtrios_results.sh

# Run the dinvestigate function of dsuite to explore arcgla and gadoga introgression in more detail.
bash run_dsuite_dinvestigate.sh

# Plot the results of the dinvestigate function in the form of manhattan plots.
bash plot_dsuite_dinvestigate_results.sh


## Tree-based introgression tests.
# Select alignments for iqtree analyses based on position (excluding inversion regions), number of variant sites and number of hemiplasies.
bash select_alignments_for_iqtree.sh

# Make maximum-likelihood trees with iqtree and root them.
bash run_iqtree.sh

# Collapse nodes connected by extremely short branches into polytomies.
bash collapse_short_branches.sh

# Analyse the sets of exon and gene trees for asymmetry of alternative topologies.
bash analyze_tree_asymmetry.sh

# Plot the asymmetries in alternative topologies as a heatmap.
bash plot_tree_asymmetry.sh

# Run iqtree once again to compare the likelihoods of three different hypotheses.
bash run_iqtree_constrained.sh

# Summarize the likelihoods of the three different hypotheses.
bash summarize_constrained_iqtree_analyses.sh

# Plot the differences between the best and second-best likelihood, per best-supported hypothesis.
bash plot_constrained_likelihoods.sh
