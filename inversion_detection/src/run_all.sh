# michaelm Mon Jan 16 13:33:43 CET 2017

## Calculate ld from a total of 100 samples.

# Calculate ld along each chromosome with plink.
bash calculate_ld.sh

# Generate plots of ld for each chromosome.
bash plot_ld.sh


## Use PacBio data of a stationary cod and a haddock specimen to identify contiguous regions via BLAST.

# Generate a BLAST databases for the stationary assembly and for the haddock assembly.
bash make_blast_dbs.sh

# Extract the target regions from the gadMor2 assembly.
bash extract_target_regions.sh

# Run BLAST searches for all extracted target regions.
bash run_gadMor2_blast.sh

# Filter the BLAST output.
bash filter_gadMor2_blast_out.sh

# Mark missing and repetitive regions in the target region fasta files.
bash characterize_features.sh

# Plot the output as SVG.
bash plot_filtered_gadMor2_blast_out.sh

# Manually generate file ../res/manual/tables/candidate_contigs.txt summarizing contigs that potentially bridge breakpoints.
echo "File ../res/manual/tables/candidate_contigs.txt must be created manually. Continue afterwards with the remaining commands."
exit

# Run BLAST searches with the candidate scaffolds as queries and gadMor2 linkage groups as databases.
bash run_contig_blast.sh

# Filter the BLAST output.
bash filter_contig_blast_out.sh
