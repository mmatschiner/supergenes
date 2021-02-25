#!/bin/bash

# m_matschiner Thu Jun 14 12:14:51 CEST 2018

# Job name:
#SBATCH --job-name=cctpll
#
# Wall clock limit:
#SBATCH --time=168:00:00
#
# Processor and memory usage:
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G

## Set up the job environment.
source /cluster/bin/jobsetup
module load ruby/2.1.5
module load python2/2.7.10.gnu

## Feedback.
echo
echo "filter_genes_by_exon_tree_congruence.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
indir=${1}
outdir=${2}
table=${3}
min_n_exons_per_gene=${4}
raxml_bin=${5}
concaterpillar_dir=${6}
n_cpus=${7}

# Use a ruby script to run concaterpillar.
ruby filter_genes_by_exon_tree_congruence.rb ${indir} ${outdir} ${table} ${min_n_exons_per_gene} ${raxml_bin} ${concaterpillar_dir} ${n_cpus}
