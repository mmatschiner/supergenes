# m_matschiner Wed Oct 23 15:09:02 CEST 2019

# Get the sequencing reads for all samples except atlantic cod from arnason and halldorsdottir.
bash get_reads.sh

# Run kollector for local alignment of mitochondrial and nuclear query genes.
bash run_kollector.sh

# Get the gadMor2 assembly.
bash get_ref.sh

# Map all samples to the gadmor2 reference to allow extracing mitochondrial reads.
bash map_to_gadmor2.sh

# Extract all reads mapping to the mitochondrial genome.
bash extract_mitogenome_from_bams.sh

# Call mitochondrial variants from each sam file individually, using mpileup.
bash call_mitochondrial_variants.sh

# As a test for possible contamination, calculate the heterozygosity of all 
# mitochondrial vcf files.
bash calculate_heterozygosity.sh

# Make a mitochondrial reference sequence by extracting the corrsponding scaffold from the gadmor2 genome.
bash make_mitochondrial_reference.sh

# Prepare fastq files from the mitochondrial bam files.
bash convert_bams_to_fastqs.sh

# Generate mitochondrial assemblies with mitobim and mira.
bash make_mitochondrial_assemblies.sh

# Align the full mitochondrial genomes produced with mitobim and mira.
bash make_mitochondrial_alignment.sh

# Filter the mitochondrial alignment.
bash filter_mitochondrial_alignment_with_bmge.sh

# Run iqtree to generate a maximum-likelihood tree of the mitochondrial data.
bash run_iqtree_mitochondrial.sh


# Nuclear phylogenetics.
# Get the zebrafish genome assembly.
bash get_zebrafish_fasta.sh

#  Make blast databases for each local assembly made with kollector.
bash make_blast_dbs.sh

# Generate alignments of putatively orthologous nuclear exons.
bash find_nuclear_orthologs.sh

# Filter nuclear exon sequences by their bitscores.
bash filter_exon_sequences_by_bitscore.sh

# Exclude all nuclear exons that contain no data for the low-coverage arcgla sample.
bash filter_exons_by_missing_data.sh

# Calculate the number of nuclear exon sequences for each specimen.
bash get_missing_exons_per_specimen.sh

# Concatenate all nuclear exon alignments.
bash concatenate_exon_alignments.sh

# Run iqtree to generate a maximum-likelihood tree of the nuclear data.
bash run_iqtree_nuclear.sh
