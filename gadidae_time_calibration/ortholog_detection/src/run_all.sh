# m_matschiner Wed Jul 18 15:43:36 CEST 2018

# Get the assemblies of 10 gadidae species from the repository of malmstrøm et al.
bash get_subject_assemblies.sh

# Get read data for gadmac and gadoga from the repository of kirubakaran et al.
bash get_subject_reads.sh

# Get queries for kollector. XXX This uses the resources from the pipefishes analysis. Update when that analysis is published.
bash get_gadmor_queries.sh

# Download bmge.
bash get_software.sh

# Use kollector to assemble exon sequences of gadmac and gadoga.
bash run_kollector.sh

# Make databases for blast.
bash make_blast_dbs.sh

# Generate alignments of putative orthologs.
bash find_orthologs.sh

# Filter sequences by their bitscores.
bash filter_sequences_by_bitscore.sh

# Filter sequences by signatures of positive selection when comparing them to the outgroup.
bash filter_sequences_by_dNdS.sh

# Filter alignment sites by gap rate and entropy with bmge.
bash filter_sites_with_BMGE.sh

# Filter exons by missing data.
bash filter_exons_by_missing_data.sh

# Remove the reference sequence from all alignments.
bash remove_reference_sequence.sh

# Filter exons by gc-content variation.
bash filter_exons_by_GC_content_variation.sh

# Determine the locations of exons on the gadmor genome reference.
bash get_exon_regions.sh

# Filter genes by exon number and sort exons into directories according to genes.
bash filter_genes_by_exon_number.sh

# Filter genes by exon tree congruence, determined with concaterpillar.
bash filter_genes_by_exon_tree_congruence.sh

# Filter genes by clock variation.
bash filter_genes_by_clock_variation_1.sh
bash filter_genes_by_clock_variation_2.sh
bash filter_genes_by_clock_variation_3.sh

# Remove selected genes after visual inspection of alignments.
bash remove_genes_with_possible_misalignment.sh

# Concatenate all alignments into a single alignment in phylip format.
bash concatenate_alignments.sh
