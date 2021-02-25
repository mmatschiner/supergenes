# m_matschiner Sun Jul 22 15:32:33 CEST 2018

# Make the output directory.
mkdir -p ../res/alignments/orthologs/09
rm -rf ../res/alignments/orthologs/09/*

# Make the log directory.
mkdir -p ../log/misc

# Set the input and output directories.
indir=../res/alignments/orthologs/08
outdir=../res/alignments/orthologs/09

# Set the name of an output table.
table=../res/tables/filtered_ortholog_exons.txt

# Set the minimum number of exons that should remain per gene.
min_n_exons_per_gene=3

# Set the paths to raxml and concaterpillar.
raxml_bin=../bin/raxmlHPC
concaterpillar_dir=../bin/concaterpillar

# Set the number of cpus to be used.
n_cpus=10

# For each gene, run concaterpillar with all exons, and remove those that do not fall into the main cluster.
out=../log/misc/concaterpillar.out
rm -f ${out}
sbatch -o ${out} filter_genes_by_exon_tree_congruence.slurm ${indir} ${outdir} ${table} ${min_n_exons_per_gene} ${raxml_bin} ${concaterpillar_dir} ${n_cpus}
