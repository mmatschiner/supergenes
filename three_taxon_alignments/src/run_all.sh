# m_matschiner Wed Jun 14 14:47:16 CEST 2017

## 1. The blast-based contig-placement approach.
# Split the manually prepared files specifying the order of contigs per chromosome.
bash split_and_clean_contig_orders.sh

# Get contigs listed in the order files from the assembly and write them to one file per specimen.
bash extract_and_order_contigs.sh

# Generate fasta files that contain all contigs merged into a single sequence,
# with or without Ns as spacers.
bash merge_contigs.sh

# Download assemblies.
bash get_assemblies.sh

# Generate pairwise alignments with masa-cudalign.
bash make_pairwise_alignments.sh

# Check the alignments.
bash check_alignments.sh

## 2. The lastz-based whole-genome-alignment approach.
# Download scripts by jim kent.
bash download_scripts.sh

# Generate masked version of the three assemblies, to be used for whole-genome alignment.
bash create_masked_assemblies.sh

# Split the masked version of the reference into separate files per lg.
bash split_masked_ref.sh

# Generate a whole-genome alignment.
bash make_whole_genome_alignment.sh

# Convert the alignments written by masa-cudalign and lastz to fasta format.
bash convert_to_fastas.sh

# Combine all alignments per lg and thread them to the gadMor2 reference alignment.
bash thread_alignments.sh

# Refine the threaded alignment by running mafft in sliding windows.
bash refine_threaded_alignments.sh


## 3. The mapping-based approach.
# Generate dictionary and index files for the reference.
bash prepare.sh

# Download fastq files.
bash get_fastqs.sh

# Map reads of gadMor2, gadMor_Stat, and melAeg against the reference to generate bam files.
bash map.sh

# Make a fasta consensus sequence for each bam file.
bash convert_bams_to_fastas.sh

# Remove unplaced and mtdna sequences.
bash remove_unplaced_seqs.sh


## 4. Finalizing the alignment.
# Finish the alignments by making strict-consensus sequences for gadMor_Stat and melAeg.
bash finish_alignments.sh

# Calculate pairwise genetic distances per linkage groups.
bash get_distances_per_lg.sh

# Calculate pairwise genetic distances in sliding windows.
bash get_distances_in_sliding_windows.sh

# Get number of complete alignment sites.
bash get_number_of_complete_alignment_sites.sh

# Make sliding-window plots of the pairwise genetic distances. 
bash plot_window_distances.sh

# Make a mask file.
bash make_mask.sh

# Get number of complete alignment sites.
bash get_number_of_complete_alignment_sites.sh


## Calculate and plot mutation load.
# Calculate mutation load.
bash calculate_mutational_load.sh

# Combine tables with estimates of mutation load.
bash combine_mutation_load_tables.sh

# Plot mutation load.
bash plot_mutational_load.sh


## Compare vitellogenin sequences.
bash compare_vitellogenins.sh
